import pytest
from openff.toolkit import Molecule, Topology
from pkg_resources import resource_filename

from openff.pablo._pdb import topology_from_pdb
from openff.pablo._utils import sort_tuple


@pytest.mark.parametrize(
    ("pdbfile", "unknown_smiles"),
    [
        (
            "data/5ap1_prepared.pdb",
            ["O=C([O-])Cn1cc(cn1)c2ccc(cc2OCC#N)Nc3ccc(c(n3)NC4CCCCC4)C#N"],
        ),
    ],
)
def test_connectivity_and_atom_order_and_net_residue_charge_and_metadata_matches_legacy(
    pdbfile: str,
    unknown_smiles: list[str],
):
    unknown_molecules = [Molecule.from_smiles(s) for s in unknown_smiles]
    filename = resource_filename(__name__, pdbfile)

    pablo_top = topology_from_pdb(
        filename,
        unknown_molecules=unknown_molecules,
    )
    legacy_top: Topology = Topology.from_pdb(
        filename,
        unique_molecules=unknown_molecules,
    )

    assert pablo_top.n_molecules == legacy_top.n_molecules
    for pablo_mol, legacy_mol in zip(pablo_top.molecules, legacy_top.molecules):
        assert pablo_mol.n_atoms == legacy_mol.n_atoms
        for pablo_atom, legacy_atom in zip(pablo_mol.atoms, legacy_mol.atoms):
            assert pablo_atom.name == legacy_atom.name
            assert pablo_atom.symbol == legacy_atom.symbol
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

    for pablo_res, legacy_res in zip(pablo_top.residues, legacy_top.residues):
        pablo_res_charge, legacy_res_charge = 0, 0
        for pablo_atom, legacy_atom in zip(pablo_res.atoms, legacy_res.atoms):
            pablo_res_charge += pablo_atom.formal_charge
            legacy_res_charge += legacy_atom.formal_charge
        assert pablo_res_charge == legacy_res_charge
