import gzip
from typing import no_type_check
from urllib.request import urlopen

from openmm.app import PDBFile
from pdbfixer import PDBFixer

PDBID = "3h34"


def fixer_from_pdbid(pdbid: str) -> PDBFixer:
    """
    Load a PDBFixer object from a PDBID.

    Differs from ``PDBFixer(pdbid=id)`` because this function downloads the PDB
    file in gzipped format, which reduces network usage, whereas the
    ``PDBFixer`` constructor downloads in uncompressed .PDB format.
    """
    with urlopen(f"https://files.rcsb.org/download/{pdbid}.pdb.gz") as gzfile:
        pdbfile = gzip.GzipFile(fileobj=gzfile)
        fixer = PDBFixer(pdbfile=pdbfile)
    return fixer


@no_type_check
def main():
    fixer = fixer_from_pdbid(PDBID)

    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    fnout = f"{PDBID}_prepared.pdb"

    with open(fnout, mode="w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


if __name__ == "__main__":
    main()
