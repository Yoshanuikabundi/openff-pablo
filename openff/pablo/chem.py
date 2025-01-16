from openff.pablo.residue import BondDefinition

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
