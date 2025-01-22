from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator, Mapping
from typing import (
    DefaultDict,
    TypeAlias,
    TypeVar,
    TypeVarTuple,
    no_type_check,
)

from openff.toolkit.topology._mm_molecule import _SimpleMolecule
from openff.toolkit.topology.molecule import MoleculeLike
from openff.toolkit.utils import UndefinedStereochemistryError
from openff.units import unit
from pint import Quantity

__all__ = [
    "default_dict",
    "unwrap",
    "sort_tuple",
    "flatten",
    "with_neighbours",
    "float_or_unknown",
    "dec_hex",
    "int_or_none",
    "cryst_to_box_vectors",
    "assign_stereochemistry_from_3d",
    "__UNSET__",
]

T = TypeVar("T")
U = TypeVar("U")
V = TypeVar("V")
Ts = TypeVarTuple("Ts")

CIFValue: TypeAlias = str | float | int


class __UNSET__:
    pass


def default_dict(
    default_factory: Callable[[], T],
    map: Mapping[U, V] = {},
) -> DefaultDict[U, T | V]:
    dd: DefaultDict[U, T | V] = defaultdict(default_factory)
    dd.update(map)
    return dd


def unwrap(container: Iterable[T], msg: str = "") -> T:
    """
    Unwrap an iterable only if it has a single element; raise ValueError otherwise
    """
    if msg:
        msg += ": "

    iterator = iter(container)

    try:
        value = next(iterator)
    except StopIteration:
        raise ValueError(msg + "container has no elements")

    try:
        next(iterator)
    except StopIteration:
        return value

    raise ValueError(msg + "container has multiple elements")


def sort_tuple(tup: tuple[*Ts]) -> tuple[*Ts]:
    return tuple(sorted(tup))  # type: ignore


def flatten(container: Iterable[Iterable[T]]) -> Iterator[T]:
    for inner in container:
        yield from inner


def with_neighbours(
    iterable: Iterable[T],
    default: U = None,
) -> Iterator[tuple[T | U, T, T | U]]:
    """Return each element of the iterable with its neighbours.

        abcd -> _ab, abc, bcd, cd_

    The middle element of the tuple is the current element. Missing neighbours
    are set to ``default``. The resulting sequence has the same length as the
    original iterable.
    """
    iterator = iter(iterable)

    pred: T | U = default
    this: T
    succ: T | U

    try:
        this = next(iterator)
    except StopIteration:
        return

    for succ in iterator:
        yield (pred, this, succ)
        pred = this
        this = succ

    succ = default
    yield (pred, this, succ)


def float_or_unknown(s: str) -> float | None:
    if s == "?":
        return None
    return float(s)


def dec_hex(s: str) -> int:
    """
    Interpret a string as a decimal or hexadecimal integer.

    For a string of length n, the string is interpreted as decimal if the value
    is < 10^n. This makes the dec_hex representation identical to a decimal
    integer, except for strings that cannot be parsed as a decimal. For these
    strings, the first hexadecimal number is interpreted as 10^n, and subsequent
    numbers continue from there. For example, in PDB files, a fixed width column
    format, residue numbers for large systems follow this representation:

        "   1" -> 1
        "   2" -> 2
        ...
        "9999" -> 9999
        "A000" -> 10000
        "A001" -> 10001
        ...
        "A009" -> 10009
        "A00A" -> 10010
        "A00B" -> 10011
        ...
        "A00F" -> 10015
        "A010" -> 10016
        ...
        "FFFF" -> 34575

    Strings that can be interpreted as hex but do not have a leading hex digit
    greater than 9 are a value error, as are strings that cannot be interpreted
    as either decimal or hexadecimal integers:

        "10A2" -> ValueError
        " 2.3" -> ValueError
        "hiya" -> ValueError
    """

    try:
        return int(s, 10)
    except ValueError:
        n = len(s)
        parsed_as_hex = int(s, 16)
        smallest_hex: int = 0xA * 16 ** (n - 1)
        if parsed_as_hex < smallest_hex:
            raise ValueError("hex values must have leading digit greater than 9")
        largest_dec: int = 10**n - 1
        return parsed_as_hex - smallest_hex + largest_dec + 1


def int_or_none(s: str) -> int | None:
    if s == "":
        return None
    else:
        return int(s)


def cryst_to_box_vectors(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> Quantity:
    @no_type_check
    def inner(a, b, c, alpha, beta, gamma):
        import openmm.unit
        from openmm.app.internal.unitcell import computePeriodicBoxVectors
        from openmm.unit import nanometer as openmm_unit_nanometer

        box_vectors = computePeriodicBoxVectors(
            openmm.unit.Quantity(a, openmm.unit.angstrom),
            openmm.unit.Quantity(b, openmm.unit.angstrom),
            openmm.unit.Quantity(c, openmm.unit.angstrom),
            openmm.unit.Quantity(alpha, openmm.unit.degree),
            openmm.unit.Quantity(beta, openmm.unit.degree),
            openmm.unit.Quantity(gamma, openmm.unit.degree),
        )
        return box_vectors.value_in_unit(openmm_unit_nanometer) * unit.nanometer

    return inner(a, b, c, alpha, beta, gamma)  # type: ignore


def assign_stereochemistry_from_3d(molecule: MoleculeLike):
    @no_type_check
    def inner(molecule):
        from rdkit.Chem import AssignStereochemistryFrom3D, BondStereo

        if isinstance(molecule, _SimpleMolecule):
            # SimpleMolecules do not store stereo info
            return

        rdmol = molecule.to_rdkit()
        AssignStereochemistryFrom3D(rdmol, confId=0, replaceExistingTags=True)

        for offatom, rdatom in zip(molecule.atoms, rdmol.GetAtoms()):
            stereochemistry = None
            if rdatom.HasProp("_CIPCode"):
                stereo_code = rdatom.GetProp("_CIPCode")
                if stereo_code == "R":
                    stereochemistry = "R"
                elif stereo_code == "S":
                    stereochemistry = "S"
                else:
                    raise UndefinedStereochemistryError(
                        "In from_pdb: Expected atom stereochemistry of R or S. "
                        f"Got {stereo_code} instead.",
                    )
            offatom._stereochemistry = stereochemistry

        for offbond, rdbond in zip(molecule.bonds, rdmol.GetBonds()):
            stereochemistry = None
            tag = rdbond.GetStereo()
            if tag == BondStereo.STEREOZ:
                stereochemistry = "Z"
            elif tag == BondStereo.STEREOE:
                stereochemistry = "E"
            elif tag == BondStereo.STEREOTRANS or tag == BondStereo.STEREOCIS:
                raise ValueError(
                    f"Expected RDKit bond stereochemistry of E or Z, got {tag} instead",
                )
            offbond._stereochemistry = stereochemistry

    inner(molecule)
