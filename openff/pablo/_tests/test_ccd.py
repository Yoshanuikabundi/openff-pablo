import pytest

from openff.pablo.ccd import CCD_RESIDUE_DEFINITION_CACHE


@pytest.mark.parametrize(
    "resname",
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "DG",
        "DA",
        "DT",
        "DC",
        "G",
        "A",
        "U",
        "C",
        "ACE",
        "NME",
    ],
)
def test_ccdcache_can_load_common_residues(resname: str):
    CCD_RESIDUE_DEFINITION_CACHE[resname]
