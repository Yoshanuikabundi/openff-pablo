import pytest
from pkg_resources import resource_filename

from openff.pablo._pdb_data import PdbData


@pytest.mark.parametrize("pdbfile", ["data/5ap1_prepared.pdb"])
def test_can_load_pdb_file_as_data(pdbfile: str):
    data = PdbData.from_file(resource_filename(__name__, pdbfile))
    assert len(data.name) == 53520
