"""
Patches to add essential features to the CCD.
"""

import dataclasses
from copy import deepcopy
from itertools import combinations

from .._utils import unwrap
from ..residue import (
    AtomDefinition,
    BondDefinition,
    ResidueDefinition,
)
from ._ccdcache import PEPTIDE_BOND

__all__ = [
    "PROTONATION_VARIANTS",
    "ATOM_NAME_SYNONYMS",
    "fix_caps",
    "add_protonation_variants",
    "add_synonyms",
    "disambiguate_alt_ids",
]


PROTONATION_VARIANTS: dict[str, list[str]] = {
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
"""Map from residue name to a list atom names of abstractable hydrogens.

Each 3-tuple specifies an atom name to remove. This atom must have exactly one
bond. A variant residue definition is created with that atom and bond removed,
and the formal charge of the bonded atom reduced by one.

Note that all combinations of deprotonations are generated; this means a residue
with ``n`` abstractable hydrogens will have ``2**n`` variants."""

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
    """Add protonation variants from the PROTONATION_VARIANTS constant"""
    deprotonations: list[tuple[str, str]] = []
    for hydrogen in PROTONATION_VARIANTS.get(res.residue_name, []):
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

            variants.append(dataclasses.replace(res, atoms=atoms, bonds=bonds))

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
    """
    clashes: list[int] = []
    canonical_names = {atom.name for atom in res.atoms}
    for i, atom in enumerate(res.atoms):
        for synonym in atom.synonyms:
            if synonym in canonical_names:
                clashes.append(i)

    if clashes:
        old_to_new = {}
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
        )
        return [res1, res2]
    else:
        return [res]
