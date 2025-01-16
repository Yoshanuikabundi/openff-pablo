import pytest
from pkg_resources import resource_filename

from openff.pablo._pdb_data import PdbData
from openff.pablo.exceptions import UnknownOrAmbiguousSerialInConectError


@pytest.mark.parametrize("pdbfile", ["data/5ap1_prepared.pdb"])
def test_can_load_pdb_file_as_data(pdbfile: str):
    data = PdbData.from_file(resource_filename(__name__, pdbfile))
    assert len(data.name) == 53520


def test_process_conects_produces_indices():
    serial_to_index: dict[int, list[int]] = {3: [0], 4: [1]}
    lines: list[str] = ["CONECT    3    4"]

    conects = PdbData._process_conects(
        lines,
        serial_to_index,
        conects=[set(), set()],
    )

    assert conects == [{1}, {0}]


def test_process_conects_raises_when_ambiguous():
    serial_to_index: dict[int, list[int]] = {3: [0], 4: [1, 2]}
    lines: list[str] = ["CONECT    3    4"]

    with pytest.raises(UnknownOrAmbiguousSerialInConectError):
        PdbData._process_conects(
            lines,
            serial_to_index,
            conects=[set(), set(), set()],
        )


def test_process_conects_raises_when_serial_missing():
    serial_to_index: dict[int, list[int]] = {3: [0], 4: [1, 2]}
    lines: list[str] = ["CONECT    3    5"]

    with pytest.raises(UnknownOrAmbiguousSerialInConectError):
        PdbData._process_conects(
            lines,
            serial_to_index,
            conects=[set(), set(), set()],
        )
