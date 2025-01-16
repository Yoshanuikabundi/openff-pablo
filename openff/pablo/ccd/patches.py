"""
Patches to add essential features to the CCD.
"""

import dataclasses
from copy import deepcopy
from itertools import combinations

from openff.pablo.chem import DISULFIDE_BOND

from .._utils import flatten, unwrap
from ..residue import (
    AtomDefinition,
    BondDefinition,
    ResidueDefinition,
)
from ._ccdcache import PEPTIDE_BOND

__all__ = [
    "ACIDIC_PROTONS",
    "BASIC_ATOMS",
    "ATOM_NAME_SYNONYMS",
    "fix_caps",
    "add_protonation_variants",
    "add_synonyms",
    "disambiguate_alt_ids",
]


ACIDIC_PROTONS: dict[str, list[str]] = {
    "ALA": ["HXT", "H2"],
    "ARG": ["HXT", "H2", "HH12"],
    "ASN": ["HXT", "H2"],
    "ASP": ["HXT", "H2", "HD2"],
    "CYS": ["HXT", "H2", "HG"],
    "GLN": ["HXT", "H2"],
    "GLU": ["HXT", "H2", "HE2"],
    "GLY": ["HXT", "H2"],
    # TODO: Special case HIS so that the neutral form is not a zwitterion
    "HIS": ["HXT", "H2", "HD1", "HE2"],
    "ILE": ["HXT", "H2"],
    "LEU": ["HXT", "H2"],
    "LYS": ["HXT", "H2", "HZ3"],
    "MET": ["HXT", "H2"],
    "PHE": ["HXT", "H2"],
    "PRO": ["HXT"],
    "SER": ["HXT", "H2", "HG"],
    "THR": ["HXT", "H2", "HG1"],
    "TRP": ["HXT", "H2", "HE1"],
    "TYR": ["HXT", "H2", "HH"],
    "VAL": ["HXT", "H2"],
}
"""Map from residue name to a list of atom names of abstractable hydrogens.

Each 3-tuple specifies an atom name to remove. This atom must have exactly one
bond. A variant residue definition is created with that atom and bond removed,
and the formal charge of the bonded atom reduced by one.

Note that all combinations of deprotonations are generated; this means a residue
with ``n`` abstractable hydrogens will have ``2**n`` variants."""


BASIC_ATOMS: dict[str, list[tuple[str, str]]] = {
    "ALA": [("N", "H3")],
    "ARG": [("N", "H3")],
    "ASN": [("N", "H3")],
    "ASP": [("N", "H3")],
    "CYS": [("N", "H3")],
    "GLN": [("N", "H3")],
    "GLU": [("N", "H3")],
    "GLY": [("N", "H3")],
    "HIS": [("N", "H3")],
    "ILE": [("N", "H3")],
    "LEU": [("N", "H3")],
    "LYS": [("N", "H3")],
    "MET": [("N", "H3")],
    "PHE": [("N", "H3")],
    "PRO": [],
    "SER": [("N", "H3")],
    "THR": [("N", "H3")],
    "TRP": [("N", "H3")],
    "TYR": [("N", "H3")],
    "VAL": [("N", "H3")],
}
"""Protonation variants that add an atom to the CCD.

Each 3-tuple specifies an atom name to protonate and the name of the added
proton. A variant residue definition is created with that atom and bond added,
and the formal charge of the atom increased by one."""

ATOM_NAME_SYNONYMS = {
    "NME": {"HN2": ["H"]},
    "NA": {"NA": ["Na"]},
    "CL": {"CL": ["Cl"]},
}
"""Map from residue name and then canonical atom name to a list of synonyms"""


def fix_caps(res: ResidueDefinition) -> list[ResidueDefinition]:
    """
    Fix ``"NON-POLYMER"`` residues so they can be used as caps for peptides.
    """

    return [
        dataclasses.replace(
            res,
            linking_bond=PEPTIDE_BOND,
            atoms=[
                dataclasses.replace(
                    atom,
                    leaving=True if atom.name == "H" else atom.leaving,
                )
                for atom in res.atoms
            ],
        ),
    ]


def add_protonation_variants(res: ResidueDefinition) -> list[ResidueDefinition]:
    """
    Add protonation variants from the ACIDIC_PROTONS and BASIC_ATOMS constants

    Note that all combinations of protonations and deprotonations are generated;
    this means a residue with ``n`` abstractable hydrogens and ``m`` acidic atoms
    will have ``2**(n+m)`` variants.
    """
    return list(
        flatten(
            add_protonated_variants(deprotonated_variant)
            for deprotonated_variant in add_deprotonated_variants(res)
        ),
    )


def add_deprotonated_variants(res: ResidueDefinition) -> list[ResidueDefinition]:
    """Add protonation variants from the ACIDIC_PROTONS constant"""
    deprotonations: list[tuple[str, str]] = []
    for hydrogen in ACIDIC_PROTONS.get(res.residue_name, []):
        bonded_atoms: list[str] = []
        for bond in res.bonds:
            if bond.atom1 == hydrogen:
                bonded_atoms.append(bond.atom2)
            elif bond.atom2 == hydrogen:
                bonded_atoms.append(bond.atom1)

        if len(bonded_atoms) != 1:
            raise ValueError("should be exactly 1 bonded atom to abstracted proton")
        deprotonations.append((hydrogen, bonded_atoms[0]))

    variants: list[ResidueDefinition] = [res]
    for i in range(len(deprotonations)):
        for combination in combinations(deprotonations, i + 1):
            hydrogens, partners = zip(*combination)

            bonds: list[BondDefinition] = [
                bond
                for bond in res.bonds
                if bond.atom1 not in hydrogens and bond.atom2 not in hydrogens
            ]

            atoms: list[AtomDefinition] = []
            for atom in res.atoms:
                if atom.name in partners:
                    atoms.append(dataclasses.replace(atom, charge=atom.charge - 1))
                elif atom.name in hydrogens:
                    if atom.symbol != "H":
                        raise ValueError(
                            "Elements of PROTONATION_VARIANTS values must be hydrogens",
                        )
                else:
                    atoms.append(atom)

            variants.append(
                dataclasses.replace(
                    res,
                    atoms=atoms,
                    bonds=bonds,
                    description=res.description + f" -{' -'.join(hydrogens)}",
                ),
            )

    return variants


def add_protonated_variants(res: ResidueDefinition) -> list[ResidueDefinition]:
    """Add protonation variants from the BASIC_ATOMS constant"""
    protonations = BASIC_ATOMS.get(res.residue_name, [])

    variants: list[ResidueDefinition] = [res]
    for i in range(len(protonations)):
        for combination in combinations(protonations, i + 1):
            bonds = [*res.bonds]
            atoms = [*res.atoms]
            for heavy_atom, hydrogen in combination:
                bonds.append(
                    BondDefinition(
                        heavy_atom,
                        hydrogen,
                        order=1,
                        aromatic=False,
                        stereo=None,
                    ),
                )

                atoms.append(
                    AtomDefinition(
                        name=hydrogen,
                        synonyms=(),
                        symbol="H",
                        leaving=False,
                        charge=0,
                        aromatic=False,
                        stereo=None,
                    ),
                )
                for i, atom in enumerate(res.atoms):
                    if atom.name == heavy_atom:
                        atoms[i] = dataclasses.replace(atom, charge=atom.charge + 1)

            hydrogens, _partners = zip(*combination)
            variants.append(
                dataclasses.replace(
                    res,
                    atoms=atoms,
                    bonds=bonds,
                    description=res.description + f" +{' +'.join(hydrogens)}",
                ),
            )

    return variants


def add_synonyms(res: ResidueDefinition) -> list[ResidueDefinition]:
    """
    Patch a residue definition to include synonyms from :py:data:`ATOM_NAME_SYNONYMS`.
    """
    return [
        dataclasses.replace(
            res,
            atoms=tuple(
                dataclasses.replace(
                    atom,
                    synonyms=tuple(
                        {
                            *atom.synonyms,
                            *ATOM_NAME_SYNONYMS.get(res.residue_name, {}).get(
                                atom.name,
                                [],
                            ),
                        },
                    ),
                )
                for atom in res.atoms
            ),
        ),
    ]


def disambiguate_alt_ids(res: ResidueDefinition) -> list[ResidueDefinition]:
    """
    CCD patch: put alt atom ids in their own residue definitions if needed

    This patch should be run before other patches that add synonyms, as it
    assumes that there is at most one synonym (from the CCD alt id
    flag).

    Some CCD residues (like GLY) have alternative atom IDs that clash with
    canonical IDs for a different atom. This breaks synonyms because the
    clashing alternate ID is never assigned; the PDB file is interpreted as
    having two copies of the canonical ID atom. To fix this, we just split
    residue definitions with this clashing problem into two definitions, one
    with the canonical IDs and the other with the alternates.

    Note that if a linking atom is part of a clash, this may interfere with
    the disambiguated residue definition's ability to link to canonically named
    residue definitions.
    """
    clashes: list[int] = []
    canonical_names = {atom.name for atom in res.atoms}
    for i, atom in enumerate(res.atoms):
        for synonym in atom.synonyms:
            if synonym in canonical_names:
                clashes.append(i)

    if clashes:
        old_to_new: dict[str, str] = {}
        for atom in res.atoms:
            if atom.synonyms:
                old_to_new[atom.name] = unwrap(atom.synonyms)
            else:
                old_to_new[atom.name] = atom.name

        res1 = dataclasses.replace(
            res,
            atoms=[
                dataclasses.replace(
                    atom,
                    synonyms=[] if i in clashes else deepcopy(atom.synonyms),
                )
                for i, atom in enumerate(res.atoms)
            ],
        )
        res2 = dataclasses.replace(
            res,
            atoms=[
                dataclasses.replace(
                    atom,
                    name=old_to_new[atom.name],
                    synonyms=[] if atom.synonyms else deepcopy(atom.synonyms),
                )
                for atom in res.atoms
            ],
            bonds=[
                dataclasses.replace(
                    bond,
                    atom1=old_to_new[bond.atom1],
                    atom2=old_to_new[bond.atom2],
                )
                for bond in res.bonds
            ],
            description=res.description + "altids",
            crosslink=(
                None
                if res.crosslink is None
                else dataclasses.replace(
                    res.crosslink,
                    atom1=old_to_new[res.crosslink.atom1],
                )
            ),
            linking_bond=(
                None
                if res.linking_bond is None
                else dataclasses.replace(
                    res.linking_bond,
                    atom1=old_to_new[res.linking_bond.atom1],
                    atom2=old_to_new.get(
                        res.linking_bond.atom2,
                        res.linking_bond.atom2,
                    ),
                )
            ),
        )
        return [res1, res2]
    else:
        return [res]


def add_disulfide_crosslink(res: ResidueDefinition) -> list[ResidueDefinition]:
    if not {"HG", "SG"}.issubset({atom.name for atom in res.atoms}):
        raise ValueError(
            "Can only add disulfide crosslink to residue with HG and SG atoms",
        )

    return [
        dataclasses.replace(
            res,
            crosslink=DISULFIDE_BOND,
            atoms=[
                atom if atom.name != "HG" else dataclasses.replace(atom, leaving=True)
                for atom in res.atoms
            ],
        ),
    ]
