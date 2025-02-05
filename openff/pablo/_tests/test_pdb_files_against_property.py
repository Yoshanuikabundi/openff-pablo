import pytest

from openff.pablo._pdb import topology_from_pdb
from openff.pablo._tests.utils import get_test_data_path


@pytest.mark.slow
def test_2zuq_cross_chain_disulfide_discontinuous():
    pdbfn = get_test_data_path("prepared_pdbs/2zuq_prepared.pdb")
    with pytest.warns():
        topology = topology_from_pdb(pdbfn)

    # Correct number of molecules
    assert topology.n_molecules == 3
    # Correct number of chains loaded in correct order
    assert [elem.identifier[0] for elem in topology.chains] == ["A", "C", "B", "D"]
    # Chains belong to correct molecules
    assert all(
        elem.scheme.parent is topology.molecule(mol_idx)
        for elem, mol_idx in zip(topology.chains, [0, 0, 1, 2])
    )
    # All molecules represent a single molecule
    assert not topology.molecule(0)._has_multiple_molecules()  # type: ignore
    assert not topology.molecule(1)._has_multiple_molecules()  # type: ignore
    assert not topology.molecule(2)._has_multiple_molecules()  # type: ignore
    # First molecule's pdb_index values are contiguous except for one discontinuity
    index_offset = 0
    for i, atom in enumerate(topology.molecule(0).atoms):
        pdb_index: int = atom.metadata["pdb_index"]  # type: ignore
        if pdb_index != i and index_offset == 0:
            index_offset = pdb_index - i

        assert pdb_index == i + index_offset
    # Other molecule's pdb_index values are contiguous
    assert all(
        i == atom.metadata["pdb_index"]
        for i, atom in enumerate(
            topology.molecule(1).atoms,
            start=topology.molecule(1).atom(0).metadata["pdb_index"],  # type: ignore
        )
    )
    assert all(
        i == atom.metadata["pdb_index"]
        for i, atom in enumerate(
            topology.molecule(2).atoms,
            start=topology.molecule(2).atom(0).metadata["pdb_index"],  # type: ignore
        )
    )
