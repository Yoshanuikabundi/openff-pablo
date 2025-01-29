import pytest
from openff.toolkit import Molecule, Topology
from pkg_resources import resource_filename

from openff.pablo._pdb import topology_from_pdb
from openff.pablo._utils import sort_tuple
from openff.pablo.residue import ResidueDefinition

FAST_PDBS: list[tuple[str, list[Molecule], list[ResidueDefinition]]] = [
    (
        "data/193l_prepared.pdb",  # Filename
        [],  # unknown_molecules
        [],  # additional_substructures
    ),
    (
        "data/3cu9_vicinal_disulfide.pdb",
        [],
        [],
    ),
    (
        "data/e2_7nel.pdb",
        [
            Molecule.from_smiles(
                r"C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O",
            ),
        ],
        [],
    ),
]
SLOW_PDBS: list[tuple[str, list[Molecule], list[ResidueDefinition]]] = [
    (
        "data/5ap1_prepared.pdb",
        [
            Molecule.from_smiles(
                "O=C([O-])Cn1cc(cn1)c2ccc(cc2OCC#N)Nc3ccc(c(n3)NC4CCCCC4)C#N",
            ),
        ],
        [],
    ),
]


@pytest.mark.slow
@pytest.mark.parametrize(
    ("pdbfile", "unknown_molecules", "additional_substructures"),
    SLOW_PDBS,
)
def test_connectivity_and_atom_order_and_net_residue_charge_and_metadata_matches_legacy_slow(
    pdbfile: str,
    unknown_molecules: list[Molecule],
    additional_substructures: list[ResidueDefinition],
):
    connectivity_and_atom_order_and_net_residue_charge_and_metadata_matches_legacy(
        pdbfile,
        unknown_molecules,
        additional_substructures,
    )


@pytest.mark.parametrize(
    ("pdbfile", "unknown_molecules", "additional_substructures"),
    FAST_PDBS,
)
def test_connectivity_and_atom_order_and_net_residue_charge_and_metadata_matches_legacy_fast(
    pdbfile: str,
    unknown_molecules: list[Molecule],
    additional_substructures: list[ResidueDefinition],
):
    connectivity_and_atom_order_and_net_residue_charge_and_metadata_matches_legacy(
        pdbfile,
        unknown_molecules,
        additional_substructures,
    )


def connectivity_and_atom_order_and_net_residue_charge_and_metadata_matches_legacy(
    pdbfile: str,
    unknown_molecules: list[Molecule],
    additional_substructures: list[ResidueDefinition],
):
    filename = resource_filename(__name__, pdbfile)

    pablo_top = topology_from_pdb(
        filename,
        unknown_molecules=unknown_molecules,
        additional_substructures=additional_substructures,
        verbose_errors=True,
    )
    legacy_top: Topology = Topology.from_pdb(
        filename,
        unique_molecules=unknown_molecules,
        _additional_substructures=[
            resdef.to_openff_molecule() for resdef in additional_substructures
        ],
    )

    assert pablo_top.n_molecules == legacy_top.n_molecules
    for pablo_mol, legacy_mol in zip(pablo_top.molecules, legacy_top.molecules):
        assert pablo_mol.n_atoms == legacy_mol.n_atoms
        for pablo_atom, legacy_atom in zip(pablo_mol.atoms, legacy_mol.atoms):
            assert pablo_atom.name == legacy_atom.name
            assert pablo_atom.symbol == legacy_atom.symbol
            assert pablo_atom.formal_charge == legacy_atom.formal_charge
            for key in [
                "residue_name",
                "chain_id",
                "residue_number",
                "insertion_code",
            ]:
                assert pablo_atom.metadata[key] == legacy_atom.metadata[key]

        pablo_bonds = {
            sort_tuple((bond.atom1_index, bond.atom2_index)) for bond in pablo_mol.bonds
        }
        legacy_bonds = {
            sort_tuple((bond.atom1_index, bond.atom2_index))
            for bond in legacy_mol.bonds
        }
        assert pablo_bonds == legacy_bonds
