from pathlib import Path

import pytest
from pkg_resources import resource_filename

from openff.pablo._pdb_data import PdbData, ResidueMatch
from openff.pablo.chem import DISULFIDE_BOND, PEPTIDE_BOND
from openff.pablo.residue import AtomDefinition, BondDefinition, ResidueDefinition


@pytest.fixture(
    params=[
        "data/5ap1_prepared.pdb",
        "data/193l_prepared.pdb",
        "data/3cu9_vicinal_disulfide.pdb",
    ],
)
def pdbfn(request: pytest.FixtureRequest) -> Path:
    return Path(resource_filename(__name__, request.param))


@pytest.fixture
def hewl_data() -> PdbData:
    return PdbData.from_file(resource_filename(__name__, "data/193l_prepared.pdb"))


@pytest.fixture
def vicinal_disulfide_data() -> PdbData:
    return PdbData.from_file(
        resource_filename(__name__, "data/3cu9_vicinal_disulfide.pdb"),
    )


@pytest.fixture
def cys_data() -> PdbData:
    return PdbData.parse_pdb(
        lines=[
            "ATOM      1  N   CYS     2      -0.664  -1.578  -0.633  1.00  0.00           N  ",
            "ATOM      2  CA  CYS     2      -0.105  -0.229  -0.633  1.00  0.00           C  ",
            "ATOM      3  C   CYS     2       1.405  -0.269  -0.633  1.00  0.00           C  ",
            "ATOM      4  O   CYS     2       2.034  -1.337  -0.633  1.00  0.00           O  ",
            "ATOM      5  CB  CYS     2      -0.674   0.537   0.577  1.00  0.00           C  ",
            "ATOM      6  SG  CYS     2      -0.163   2.270   0.535  1.00  0.00           S  ",
            "ATOM      7  OXT CYS     2       2.131   0.963  -0.633  1.00  0.00           O  ",
            "ATOM      8  H   CYS     2      -0.036  -2.457  -0.633  1.00  0.00           H  ",
            "ATOM      9  H2  CYS     2      -1.244  -1.644   0.190  1.00  0.00           H  ",
            "ATOM     10  HA  CYS     2      -0.420   0.277  -1.564  1.00  0.00           H  ",
            "ATOM     11  HB2 CYS     2      -0.349   0.082   1.535  1.00  0.00           H  ",
            "ATOM     12  HB3 CYS     2      -1.778   0.509   0.585  1.00  0.00           H  ",
            "ATOM     13  HG  CYS     2       0.749   2.191   1.501  1.00  0.00           H  ",
            "ATOM     14  HXT CYS     2       2.294   1.239   0.272  1.00  0.00           H  ",
        ],
    )


@pytest.fixture
def cys_def() -> ResidueDefinition:
    atoms = (
        AtomDefinition.with_defaults(name="N", symbol="N"),
        AtomDefinition.with_defaults(name="CA", symbol="C", stereo="R"),
        AtomDefinition.with_defaults(name="C", symbol="C"),
        AtomDefinition.with_defaults(name="O", symbol="O"),
        AtomDefinition.with_defaults(name="CB", symbol="C"),
        AtomDefinition.with_defaults(name="SG", symbol="S"),
        AtomDefinition.with_defaults(name="OXT", symbol="O", leaving=True),
        AtomDefinition.with_defaults(name="H", symbol="H"),
        AtomDefinition.with_defaults(name="H2", symbol="H", leaving=True),
        AtomDefinition.with_defaults(name="HA", symbol="H"),
        AtomDefinition.with_defaults(name="HB2", symbol="H"),
        AtomDefinition.with_defaults(name="HB3", symbol="H"),
        AtomDefinition.with_defaults(name="HG", symbol="H", leaving=True),
        AtomDefinition.with_defaults(name="HXT", symbol="H", leaving=True),
    )

    bonds = (
        BondDefinition.with_defaults("N", "CA"),
        BondDefinition.with_defaults("N", "H"),
        BondDefinition.with_defaults("N", "H2"),
        BondDefinition.with_defaults("CA", "C"),
        BondDefinition.with_defaults("CA", "CB"),
        BondDefinition.with_defaults("CA", "HA"),
        BondDefinition.with_defaults("C", "O", order=2),
        BondDefinition.with_defaults("C", "OXT"),
        BondDefinition.with_defaults("CB", "SG"),
        BondDefinition.with_defaults("CB", "HB2"),
        BondDefinition.with_defaults("CB", "HB3"),
        BondDefinition.with_defaults("SG", "HG"),
        BondDefinition.with_defaults("OXT", "HXT"),
    )

    return ResidueDefinition(
        atoms=atoms,
        bonds=bonds,
        crosslink=DISULFIDE_BOND,
        linking_bond=PEPTIDE_BOND,
        description="CYSTEINE",
        residue_name="CYS",
    )


@pytest.fixture
def cys_def_deprotonated_sidechain(
    cys_def: ResidueDefinition,
) -> ResidueDefinition:
    atoms: list[AtomDefinition] = []
    for atom in cys_def.atoms:
        if atom.name == "HG":
            pass
        elif atom.name == "SG":
            atoms.append(AtomDefinition.with_defaults(name="SG", symbol="S", charge=-1))
        else:
            atoms.append(atom)

    return ResidueDefinition(
        atoms=tuple(atoms),
        bonds=tuple(
            bond for bond in cys_def.bonds if "HG" not in [bond.atom1, bond.atom2]
        ),
        crosslink=DISULFIDE_BOND,
        linking_bond=PEPTIDE_BOND,
        description="CYSTEINE",
        residue_name="CYS",
    )


@pytest.fixture
def hoh_def() -> ResidueDefinition:
    atoms = (
        AtomDefinition.with_defaults(name="H1", symbol="H"),
        AtomDefinition.with_defaults(name="H2", symbol="H"),
        AtomDefinition.with_defaults(name="O", symbol="O"),
    )

    bonds = (
        BondDefinition.with_defaults("O", "H2"),
        BondDefinition.with_defaults("O", "H1"),
    )

    return ResidueDefinition(
        atoms=atoms,
        bonds=bonds,
        crosslink=None,
        linking_bond=None,
        description="water",
        residue_name="HOH",
    )


@pytest.fixture
def hoh_def_with_synonyms() -> ResidueDefinition:
    atoms = (
        AtomDefinition.with_defaults(name="H1", symbol="H", synonyms=["HA"]),
        AtomDefinition.with_defaults(name="H2", symbol="H", synonyms=["HB"]),
        AtomDefinition.with_defaults(name="O", symbol="O", synonyms=["O1"]),
    )

    bonds = (
        BondDefinition.with_defaults("O", "H2"),
        BondDefinition.with_defaults("O", "H1"),
    )

    return ResidueDefinition(
        atoms=atoms,
        bonds=bonds,
        crosslink=None,
        linking_bond=None,
        description="water",
        residue_name="HOH",
    )


@pytest.fixture
def cys_match(cys_def: ResidueDefinition) -> ResidueMatch:
    return ResidueMatch(
        residue_definition=cys_def,
        crosslink=None,
        index_to_atomdef={i: atom for i, atom in enumerate(cys_def.atoms)},
    )


@pytest.fixture
def cys_match_no_leaving(cys_def: ResidueDefinition) -> ResidueMatch:
    counter = iter(range(len(cys_def.atoms)))
    return ResidueMatch(
        residue_definition=cys_def,
        crosslink=None,
        index_to_atomdef={
            next(counter): atom for atom in cys_def.atoms if not atom.leaving
        },
    )


@pytest.fixture
def cys_match_no_leaving_deprotonated_sidechain(
    cys_def_deprotonated_sidechain: ResidueDefinition,
) -> ResidueMatch:
    counter = iter(range(len(cys_def_deprotonated_sidechain.atoms)))
    return ResidueMatch(
        residue_definition=cys_def_deprotonated_sidechain,
        crosslink=None,
        index_to_atomdef={
            next(counter): atom
            for atom in cys_def_deprotonated_sidechain.atoms
            if not atom.leaving
        },
    )


@pytest.fixture
def hoh_match(hoh_def: ResidueDefinition) -> ResidueMatch:
    return ResidueMatch(
        residue_definition=hoh_def,
        crosslink=None,
        index_to_atomdef={i: atom for i, atom in enumerate(hoh_def.atoms)},
    )


@pytest.fixture(
    params=[
        "hoh_match",
        "cys_match",
        "cys_match_no_leaving",
        "cys_match_no_leaving_deprotonated_sidechain",
    ],
)
def any_match(request: pytest.FixtureRequest) -> ResidueMatch:
    return request.getfixturevalue(request.param)
