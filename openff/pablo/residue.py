"""
Classes for defining custom residues.
"""

import dataclasses
from collections.abc import Collection, Iterator, Mapping
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import dataclass
from functools import cached_property
from typing import Literal, Self
from collections.abc import Iterable

from openff.toolkit import Molecule
from openff.units import elements, unit

__all__ = [
    "AtomDefinition",
    "BondDefinition",
    "ResidueDefinition",
]

_residue_definition_skip_validation = False


@contextmanager
def _defer_residue_definition_validation():  # type: ignore[deadcode]
    global _residue_definition_skip_validation
    _residue_definition_skip_validation = True
    yield
    _residue_definition_skip_validation = False


@dataclass(frozen=True)
class AtomDefinition:
    """
    Description of an atom in a residue from the Chemical Component Dictionary (CCD).
    """

    name: str
    """The canonical name of this atom"""
    synonyms: tuple[str, ...]
    """Other names this atom can have"""
    symbol: str
    """The elemental symbol for this atom"""
    leaving: bool
    """Whether this atom is absent when some bond is formed between this residue and another"""
    charge: int
    """The formal charge of this atom"""
    aromatic: bool
    """Whether this atom is aromatic"""
    stereo: Literal["S", "R"] | None
    """The chirality of this atom"""

    @classmethod
    def with_defaults(
        cls,
        name: str,
        symbol: str,
        synonyms: Iterable[str] = (),
        leaving: bool = False,
        charge: int = 0,
        aromatic: bool = False,
        stereo: Literal["S", "R"] | None = None,
    ) -> Self:
        return cls(
            name=name,
            symbol=symbol,
            synonyms=tuple(synonyms),
            leaving=leaving,
            charge=charge,
            aromatic=aromatic,
            stereo=stereo,
        )


@dataclass(frozen=True)
class BondDefinition:
    """
    Description of a bond in a residue from the Chemical Component Dictionary (CCD).
    """

    atom1: str
    """The canonical name of the first atom in this bond"""
    atom2: str
    """The canonical name of the second atom in this bond"""
    order: int
    """The bond order of this bond (1 for single bond, 2 for double bond, etc.)"""
    aromatic: bool
    """``True`` if this bond is aromatic, ``False`` otherwise"""
    stereo: Literal["E", "Z"] | None
    """The stereochemistry of this bond."""

    def flipped(self) -> Self:
        """The same bond, but with the atoms in the opposite order"""
        return dataclasses.replace(self, atom1=self.atom2, atom2=self.atom1)

    def sorted(self) -> Self:
        """The same bond, but with the atoms in sorted order"""
        sorted_atoms = sorted([self.atom1, self.atom2])
        return dataclasses.replace(self, atom1=sorted_atoms[0], atom2=sorted_atoms[1])

    @classmethod
    def with_defaults(
        cls,
        atom1: str,
        atom2: str,
        order: int = 1,
        aromatic: bool = False,
        stereo: Literal["E", "Z"] | None = None,
    ) -> Self:
        return cls(
            atom1=atom1,
            atom2=atom2,
            order=order,
            aromatic=aromatic,
            stereo=stereo,
        )


@dataclass(frozen=True)
class ResidueDefinition:
    """
    Description of a residue from the Chemical Component Dictionary (CCD).
    """

    residue_name: str
    """The 3-letter residue code used in PDB files"""
    description: str
    """A longer description of the residue"""
    linking_bond: BondDefinition | None
    """Description of how this residue may bond to its neighbours in a polymer

    If the residue is only found as a monomer, ``None``. Otherwise, a
    :py:cls:`BondDefinition`. The ``atom1`` and ``atom2`` attributes give the
    canonical atom names of the atoms that form the linking bond. ``atom1`` is
    the name of the atom in the residue preceding the bond, and ``atom2`` is in
    the residue after the bond. Any atoms that the bond between residues
    replaces should be marked as ``atom.leaving=True``.

    A linking bond will be formed to join two residues if and only if all of the
    below are true:

    - The residues' atom records are sequential in the PDB file
    - The residues have the same chain ID
    - There is no TER record between the residues in the PDB file
    - The two residues have identical ``linking_bond`` attributes
    - All leaving atoms associated with ``linking_bond.atom1`` are absent in
    the first encountered residue, and all leaving atoms associated with
    ``linking_bond.atom2`` are absent in the latter residue. Leaving atoms
    are those that have the ``AtomDefinition.leaving`` attribute set to
    ``True``. A leaving atom is associated with atom `a` if it is bonded to
    `a`, or it is bonded to an atom associated with `a`.
    - There is at least one leaving atom associated with each linking atom

    The charge of linking atoms is not modified; any change in valence is
    accounted for via the removal of leaving atoms.
    """
    crosslink: BondDefinition | None
    """Optional description of a crosslink between this residue and another.

    A crosslink is a bond between this residue and another. Crosslinks differ
    from linking bonds primarily because they may occur between residues that
    are do not appear sequentially in the PDB file. Crosslinks cannot occur
    within a residue; use a bond or residue variant instead.

    If the residue does not form cross links, ``None``. Otherwise, a
    :py:cls:`BondDefinition`. The ``atom1`` and ``atom2`` attributes give the
    canonical atom names of the atoms that form the crosslink. ``atom1`` is
    the name of the atom in this residue, and ``atom2`` is in the other residue.
    Any atoms that the bond between residues replaces should be marked as
    ``atom.leaving=True``.

    A crosslink will be formed between two residues if and only if all of the
    below are true:

    - The two residues' ``crosslink`` attributes are identical except that their
    atom names are reversed
    - All leaving atoms associated with each residues' ``crosslink.atom1``
    attribute are absent in that residue. Leaving atoms are those that have the
    ``AtomDefinition.leaving`` attribute set to ``True``. A leaving atom is
    associated with atom `a` if it is bonded to `a`, or it is bonded to an atom
    associated with `a`.
    - There is at least one leaving atom associated with each cross-linking atom

    The charge of linking atoms is not modified; any change in valence is
    accounted for via the removal of leaving atoms.
    """
    atoms: tuple[AtomDefinition, ...]
    """The atom definitions that make up this residue"""
    bonds: tuple[BondDefinition, ...]
    """The bond definitions that make up this residue"""

    def __post_init__(self):
        if _residue_definition_skip_validation:
            return

        self._validate()

    def _validate(self):
        if (
            self.linking_bond is None
            and self.crosslink is None
            and True in {atom.leaving for atom in self.atoms}
        ):
            raise ValueError(
                f"{self.residue_name}: Leaving atoms were specified, but there is no linking bond or crosslink",
                self,
            )
        if len({atom.name for atom in self.atoms}) != len(self.atoms):
            raise ValueError(
                f"{self.residue_name}: All atoms must have unique canonical names",
            )

        all_leaving_atoms = {atom.name for atom in self.atoms if atom.leaving}
        assigned_leaving_atoms = (
            self.prior_bond_leaving_atoms
            | self.posterior_bond_leaving_atoms
            | self.crosslink_leaving_atoms
        )
        unassigned_leaving_atoms = all_leaving_atoms.difference(assigned_leaving_atoms)
        if len(unassigned_leaving_atoms) != 0:
            raise ValueError(
                f"{self.residue_name}: Leaving atoms could not be assigned to a bond: {unassigned_leaving_atoms}",
            )

    @classmethod
    def from_molecule(
        cls,
        name: str,
        molecule: Molecule,
        linking_bond: BondDefinition | None = None,
        crosslink: BondDefinition | None = None,
        description: str = "",
    ) -> Self:
        atoms: list[AtomDefinition] = []
        for atom in molecule.atoms:
            atoms.append(
                AtomDefinition(
                    name=atom.name,
                    synonyms=(),
                    symbol=atom.symbol,
                    leaving=bool(atom.metadata.get("leaving_atom")),
                    charge=atom.formal_charge.m_as(unit.elementary_charge),  # type: ignore
                    stereo=atom.stereochemistry,
                    aromatic=atom.is_aromatic,
                ),
            )
        bonds: list[BondDefinition] = []
        for bond in molecule.bonds:
            bonds.append(
                BondDefinition(
                    atom1=bond.atom1.name,
                    atom2=bond.atom2.name,
                    order=bond.bond_order,
                    aromatic=bond.is_aromatic,
                    stereo=bond.stereochemistry,
                ),
            )

        return cls(
            residue_name=name,
            description=description,
            linking_bond=linking_bond,
            crosslink=crosslink,
            atoms=tuple(atoms),
            bonds=tuple(bonds),
        )

    @classmethod
    def from_capped_molecule(
        cls,
        name: str,
        molecule: Molecule,
        leaving_atom_indices: Collection[int],
        linking_bond: BondDefinition,
        crosslink: BondDefinition | None = None,
        description: str = "",
    ) -> Self:
        molecule = deepcopy(molecule)
        for i in leaving_atom_indices:
            molecule.atom(i).metadata["leaving_atom"] = True
        return cls.from_molecule(
            name=name,
            molecule=molecule,
            linking_bond=linking_bond,
            crosslink=crosslink,
            description=description,
        )

    @classmethod
    def from_smiles(
        cls,
        name: str,
        mapped_smiles: str,
        atom_names: Mapping[int, str],
        leaving_atoms: Collection[int] = (),
        linking_bond: BondDefinition | None = None,
        crosslink: BondDefinition | None = None,
        description: str = "",
    ) -> Self:
        molecule = Molecule.from_mapped_smiles(mapped_smiles)
        leaving_atom_indices = set(leaving_atoms)
        for i, atom in enumerate(molecule.atoms, start=1):
            if i in leaving_atom_indices:
                atom.metadata["leaving_atom"] = True
            atom.name = atom_names[i]

        return cls.from_molecule(
            name=name,
            molecule=molecule,
            linking_bond=linking_bond,
            description=description,
            crosslink=crosslink,
        )

    def to_openff_molecule(self) -> Molecule:
        molecule = Molecule()
        atoms: dict[str, int] = {}
        for atom in self.atoms:
            atoms[atom.name] = molecule.add_atom(
                atomic_number=elements.NUMBERS[atom.symbol],
                formal_charge=atom.charge,
                is_aromatic=atom.aromatic,
                stereochemistry=atom.stereo,
                name=atom.name,
                metadata={
                    "residue_name": self.residue_name,
                    "leaving": atom.leaving,
                },
            )

        for bond in self.bonds:
            molecule.add_bond(
                atom1=atoms[bond.atom1],
                atom2=atoms[bond.atom2],
                bond_order=bond.order,
                is_aromatic=bond.aromatic,
                stereochemistry=bond.stereo,
            )

        molecule.properties.update(
            {
                "linking_bond": self.linking_bond,
            },
        )

        return molecule

    @cached_property
    def name_to_atom(self) -> dict[str, AtomDefinition]:
        """Map from each atoms' name and synonyms to the value of a field."""
        mapping = {atom.name: atom for atom in self.atoms}
        canonical_names = set(mapping)
        for atom in self.atoms:
            for synonym in atom.synonyms:
                if synonym in mapping and mapping[synonym] != atom:
                    raise ValueError(
                        f"synonym {synonym} degenerately defined for canonical"
                        + f" names {mapping[synonym]} and {atom.name} in"
                        + f" residue {self.residue_name}",
                    )
                if synonym in canonical_names:
                    raise ValueError(
                        f"synonym {synonym} of atom {atom.name} clashes with"
                        + f" another canonical name in residue {self.residue_name}",
                    )
                mapping[synonym] = atom
        return mapping

    def atoms_bonded_to(self, atom_name: str) -> Iterator[str]:
        for bond in self.bonds:
            if bond.atom1 == atom_name:
                yield bond.atom2
            if bond.atom2 == atom_name:
                yield bond.atom1

    def _leaving_fragment_of(self, linking_atom: str) -> Iterator[str]:
        atoms_to_check = list(self.atoms_bonded_to(linking_atom))
        checked_atoms: set[str] = set()
        while atoms_to_check:
            atom_name = atoms_to_check.pop()
            if self.name_to_atom[atom_name].leaving:
                yield atom_name
                atoms_to_check.extend(
                    filter(
                        lambda x: x not in checked_atoms,
                        self.atoms_bonded_to(atom_name),
                    ),
                )
            checked_atoms.add(atom_name)

    @cached_property
    def posterior_bond_leaving_atoms(self) -> set[str]:
        return (
            set()
            if self.linking_bond is None
            else set(self._leaving_fragment_of(self.posterior_bond_linking_atom))
        )

    @cached_property
    def prior_bond_leaving_atoms(self) -> set[str]:
        return (
            set()
            if self.linking_bond is None
            else set(self._leaving_fragment_of(self.prior_bond_linking_atom))
        )

    @cached_property
    def crosslink_leaving_atoms(self) -> set[str]:
        return (
            set()
            if self.crosslink is None
            else set(self._leaving_fragment_of(self.crosslink.atom1))
        )

    @property
    def prior_bond_linking_atom(self) -> str:
        if self.linking_bond is None:
            raise ValueError("not a linking residue")
        return self.linking_bond.atom2

    @property
    def posterior_bond_linking_atom(self) -> str:
        if self.linking_bond is None:
            raise ValueError("not a linking residue")
        return self.linking_bond.atom1
