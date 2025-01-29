"""
Chemical information for re-use across residue databases.
"""

from openff.pablo.residue import BondDefinition

__all__ = [
    "DISULFIDE_BOND",
    "PEPTIDE_BOND",
    "PHOSPHODIESTER_BOND",
]

DISULFIDE_BOND = BondDefinition(
    atom1="SG",
    atom2="SG",
    order=1,
    aromatic=False,
    stereo=None,
)

PEPTIDE_BOND = BondDefinition(
    atom1="C",
    atom2="N",
    order=1,
    aromatic=False,
    stereo=None,
)

PHOSPHODIESTER_BOND = BondDefinition(
    atom1="O3'",
    atom2="P",
    order=1,
    aromatic=False,
    stereo=None,
)
