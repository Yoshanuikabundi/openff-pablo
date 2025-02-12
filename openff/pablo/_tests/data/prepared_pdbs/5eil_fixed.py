import gzip
from io import StringIO
from typing import no_type_check
from urllib.request import urlopen

from openmm.app import PDBFile
from pdbfixer import PDBFixer

PDBID = "5eil"
"""2zuq's model includes two biological assemblies, the first in chains A, B and
C and the second in chains D, E, and F. Cross-chain disulfide bonds exist
between chains B and C and between chains E and F."""


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

    fnout = f"{PDBID}_fixed.pdb"
    with StringIO() as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
        pdblines = f.getvalue().splitlines()

    with open(fnout, mode="w") as f:
        for line in pdblines:
            fe_serial = 7179
            if line.startswith("CONECT") and str(fe_serial) in line:
                a, *bs = (
                    int(line[start:stop].strip())
                    for start, stop in [(6, 11), (11, 16), (16, 21), (21, 26), (26, 31)]
                    if len(line[start:stop].strip()) > 0
                )
                print(a, bs)
                if a == fe_serial:
                    print("skipping")
                    continue
                if fe_serial in bs:
                    line = f"CONECT{a: >5}" + "".join(
                        f"{b: >5}" for b in bs if b != fe_serial
                    )
                    print(line)
            print(line, file=f)


if __name__ == "__main__":
    main()
