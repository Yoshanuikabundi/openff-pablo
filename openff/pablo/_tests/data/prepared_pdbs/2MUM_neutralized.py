import gzip
from io import StringIO
from pathlib import Path
from typing import no_type_check
from urllib.request import urlopen

import openmm
from openmm.app import PDBFile
from pdbfixer import PDBFixer

from openff.pablo._utils import unwrap

PDBID = "2MUM"


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
def energy_minimise_fixer(
    fixer: PDBFixer,
    forcefield: openmm.app.ForceField = openmm.app.ForceField(
        "charmm36.xml",
    ),
):
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=openmm.app.CutoffNonPeriodic,
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
    simulation.minimizeEnergy(maxIterations=100_000)

    fixer.positions = simulation.context.getState(positions=True).getPositions()


def main():
    fixer = fixer_from_pdbid(PDBID)
    fixer.removeHeterogens()
    fixer.addMissingHydrogens(pH=0)

    last_res = list(fixer.topology.residues())[-1]
    last_res_oxt = unwrap(atom for atom in last_res.atoms() if atom.name == "OXT")
    last_res_hxt = fixer.topology.addAtom(
        name="HXT",
        element=openmm.app.element.hydrogen,
        residue=last_res,
    )
    fixer.topology.addBond(atom1=last_res_oxt, atom2=last_res_hxt, order=1)
    fixer.positions.append(fixer.positions[last_res_oxt.index])

    energy_minimise_fixer(fixer)

    fnout = Path(f"{PDBID}_neutralized.pdb")

    with StringIO() as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
        lines = f.getvalue().splitlines()

    with open(fnout, "w") as f:
        for line in lines:
            resname = line[17:20].strip()
            atomname = line[12:16].strip()
            if (line.startswith("ATOM  ") or line.startswith("HETATM")) and (
                (resname == "ARG" and atomname == "HH22")
                or (resname == "LYS" and atomname == "HZ3")
                or (resname == "HIS" and atomname == "HD1")
                or (atomname == "H3")
            ):
                continue
            f.write(line)
            f.write("\n")


if __name__ == "__main__":
    main()
