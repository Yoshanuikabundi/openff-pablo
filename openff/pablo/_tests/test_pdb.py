import random

import pytest
from openff.toolkit import Molecule

from openff.pablo._pdb import _add_to_molecule, _match_unknown_molecules
from openff.pablo._pdb_data import PdbData, ResidueMatch


class TestMatchUnknownMolecules:
    def test_returned_molecule_in_pdb_order(self, e2_data: PdbData):
        e2_mol = Molecule.from_smiles(
            r"C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O",
        )

        atom_indices = list(range(e2_mol.n_atoms))
        random.shuffle(atom_indices)
        e2_mol.remap({i: j for (i, j) in enumerate(atom_indices)})

        match_mol = _match_unknown_molecules(
            data=e2_data,
            indices=tuple(range(e2_mol.n_atoms)),
            unknown_molecules=[e2_mol],
        )
        assert match_mol is not None
        assert match_mol.is_isomorphic_with(e2_mol)

        for i, atom in enumerate(match_mol.atoms):
            assert atom.metadata["pdb_index"] == i
            assert atom.name == e2_data.name[i]
            assert atom.formal_charge == e2_data.charge[i]
            assert atom.symbol == e2_data.element[i]

    def test_returned_molecule_has_metadata(
        self,
        e2_data: PdbData,
    ):
        e2_mol = Molecule.from_smiles(
            r"C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O",
        )
        match_mol = _match_unknown_molecules(
            data=e2_data,
            indices=tuple(range(e2_mol.n_atoms)),
            unknown_molecules=[e2_mol],
        )
        assert match_mol is not None

        pdb_idx_to_mol_atom_idx: dict[int, int] = match_mol.properties[
            "pdb_idx_to_mol_atom_idx"
        ]
        assert pdb_idx_to_mol_atom_idx == {i: i for i in range(match_mol.n_atoms)}

        for i, atom in enumerate(match_mol.atoms):
            assert atom.metadata["residue_name"] == e2_data.res_name[i]
            assert atom.metadata["residue_number"] == str(e2_data.res_seq[i])
            assert atom.metadata["res_seq"] == e2_data.res_seq[i]
            assert atom.metadata["insertion_code"] == e2_data.i_code[i]
            assert atom.metadata["chain_id"] == e2_data.chain_id[i]
            assert atom.metadata["pdb_index"] == i
            assert atom.metadata["atom_serial"] == e2_data.serial[i]
            assert float(atom.metadata["b_factor"]) == e2_data.temp_factor[i]
            assert float(atom.metadata["occupancy"]) == e2_data.occupancy[i]
            assert "used_synonym" not in atom.metadata
            assert "canonical_name" not in atom.metadata
            assert "matched_residue_description" not in atom.metadata

    def test_returned_molecule_isomorphic_with_unknown_molecules(
        self,
        e2_data: PdbData,
    ):
        e2_mol = Molecule.from_smiles(
            r"C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O",
        )
        hoh_mol = Molecule.from_smiles("O")

        match = _match_unknown_molecules(
            data=e2_data,
            indices=tuple(range(e2_mol.n_atoms)),
            unknown_molecules=[e2_mol, hoh_mol],
        )
        assert match is not None
        assert match.is_isomorphic_with(e2_mol)

    def test_conect_to_other_residue_raises(
        self,
        vicinal_disulfide_data: PdbData,
    ):
        cys_mol = Molecule.from_smiles("N[C@H](CS)C(O)=O")
        with pytest.raises(ValueError):
            _match_unknown_molecules(
                data=vicinal_disulfide_data,
                indices=tuple(range(cys_mol.n_atoms)),
                unknown_molecules=[cys_mol],
            )

    def test_returns_none_on_failure_to_match(
        self,
        e2_data: PdbData,
    ):
        hoh_mol = Molecule.from_smiles("O")

        match = _match_unknown_molecules(
            data=e2_data,
            indices=tuple(range(44)),
            unknown_molecules=[hoh_mol],
        )
        assert match is None


def test_add_to_molecule_when_empty(cys_data: PdbData, cys_match: ResidueMatch):
    this_molecule = Molecule()
    molecules = [this_molecule]
    mol_out = _add_to_molecule(
        data=cys_data,
        molecules=molecules,
        this_molecule=this_molecule,
        residue_match=cys_match,
        res_atom_idcs=tuple(cys_match.res_atom_idcs),
        use_canonical_names=False,
    )
    assert mol_out is this_molecule
    matched, _ = Molecule.are_isomorphic(
        mol_out,
        cys_match.residue_definition.to_openff_molecule(),
        atom_stereochemistry_matching=False,
    )
    assert matched
