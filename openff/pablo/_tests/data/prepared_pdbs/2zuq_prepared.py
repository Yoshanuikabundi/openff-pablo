import gzip
from typing import no_type_check
from urllib.request import urlopen

import numpy as np
import openmm
from openmm.app import PDBFile, Topology
from pdbfixer import PDBFixer

PDBID = "2zuq"
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
def energy_minimise_fixer(
    fixer: PDBFixer,
    forcefield: openmm.app.ForceField = openmm.app.ForceField(
        "amber99sbildn.xml",
        "amber99_obc.xml",
    ),
):
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
    simulation.minimizeEnergy(maxIterations=100_000)

    fixer.positions = simulation.context.getState(positions=True).getPositions()


def reorder_chains(fixer: PDBFixer, chain_ids: list[str]):
    topology: Topology = fixer.topology
    new_top: Topology = Topology()
    positions = np.asarray(fixer.positions.value_in_unit(openmm.unit.nanometer))  # type:ignore
    old_to_new = {}
    new_positions = []
    for chain_to_place in chain_ids:
        for chain in topology.chains():
            if chain.id == chain_to_place:
                new_chain = new_top.addChain(chain.id)
                for residue in chain.residues():
                    new_residue = new_top.addResidue(
                        name=residue.name,
                        chain=new_chain,
                        id=residue.id,
                        insertionCode=residue.insertionCode,
                    )
                    for atom in residue.atoms():
                        new_atom = new_top.addAtom(
                            name=atom.name,
                            element=atom.element,
                            residue=new_residue,
                            id=atom.id,
                            formalCharge=atom.formalCharge,
                        )
                        old_to_new[atom] = new_atom
                        new_positions.append(positions[atom.index])
    for bond in topology.bonds():
        if bond.atom1 not in old_to_new or bond.atom2 not in old_to_new:
            continue
        new_top.addBond(
            old_to_new[bond.atom1],
            old_to_new[bond.atom2],
            type=bond.type,
            order=bond.order,
        )

    fixer.topology = new_top
    fixer.positions = new_positions * openmm.unit.nanometer  # type:ignore


@no_type_check
def main():
    fixer = fixer_from_pdbid(PDBID)

    # Remove most of additional biological assembly, but preserve a chain for
    # testing chains following broken molecule
    fixer.removeChains(chainIds=["E", "F"])

    fixer.missingResidues = {}
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    energy_minimise_fixer(fixer)

    reorder_chains(fixer, ["B", "A", "C", "D"])

    fnout = f"{PDBID}_prepared.pdb"

    with open(fnout, mode="w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


if __name__ == "__main__":
    main()
