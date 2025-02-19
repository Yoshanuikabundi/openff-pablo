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

# TODO: Fill in this data
_LINKING_TYPES: dict[str, BondDefinition | None] = {
    # "D-beta-peptide, C-gamma linking".upper(): [],
    # "D-gamma-peptide, C-delta linking".upper(): [],
    # "D-peptide COOH carboxy terminus".upper(): [],
    # "D-peptide NH3 amino terminus".upper(): [],
    # "D-peptide linking".upper(): [],
    # "D-saccharide".upper(): [],
    # "D-saccharide, alpha linking".upper(): [],
    # "D-saccharide, beta linking".upper(): [],
    # "DNA OH 3 prime terminus".upper(): [],
    # "DNA OH 5 prime terminus".upper(): [],
    "DNA linking".upper(): PHOSPHODIESTER_BOND,
    "L-DNA linking".upper(): PHOSPHODIESTER_BOND,
    "L-RNA linking".upper(): PHOSPHODIESTER_BOND,
    # "L-beta-peptide, C-gamma linking".upper(): [],
    # "L-gamma-peptide, C-delta linking".upper(): [],
    # "L-peptide COOH carboxy terminus".upper(): [],
    # "L-peptide NH3 amino terminus".upper(): [],
    "L-peptide linking".upper(): PEPTIDE_BOND,
    # "L-saccharide".upper(): [],
    # "L-saccharide, alpha linking".upper(): [],
    # "L-saccharide, beta linking".upper(): [],
    # "RNA OH 3 prime terminus".upper(): [],
    # "RNA OH 5 prime terminus".upper(): [],
    "RNA linking".upper(): PHOSPHODIESTER_BOND,
    "non-polymer".upper(): None,
    # "other".upper(): [],
    "peptide linking".upper(): PEPTIDE_BOND,
    "peptide-like".upper(): PEPTIDE_BOND,
    # "saccharide".upper(): [],
}
"""Map from the CCD's linking types to the bond formed between two such monomers"""
