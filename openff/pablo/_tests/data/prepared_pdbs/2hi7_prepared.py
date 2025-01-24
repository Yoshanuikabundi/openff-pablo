import gzip
from typing import no_type_check
from urllib.request import urlopen

import openmm
from openmm.app import PDBFile
from pdbfixer import PDBFixer

PDBID = "2hi7"


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
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    forcefield = openmm.app.ForceField("amber99sbildn.xml", "amber99_obc.xml")
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=openmm.app.CutoffPeriodic,
        nonbondedCutoff=1.6 * openmm.unit.nanometer,
        constraints=None,
        rigidWater=True,
    )
    integrator = openmm.LangevinMiddleIntegrator(
        300 * openmm.unit.kelvin,
        5.0 / openmm.unit.picosecond,
        1.0 * openmm.unit.femtosecond,
    )
    # integrator.setConstraintTolerance(0.00001)
    simulation = openmm.app.Simulation(
        fixer.topology,
        system,
        integrator,
    )

    simulation.context.setPositions(fixer.positions)

    print("Minimizing...")
    simulation.minimizeEnergy(maxIterations=10_000)

    fnout = f"{PDBID}_prepared.pdb"

    with open(fnout, mode="w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


if __name__ == "__main__":
    main()
