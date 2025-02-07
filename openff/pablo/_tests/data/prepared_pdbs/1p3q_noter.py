from io import StringIO
from pathlib import Path
from typing import no_type_check

import openmm
from openmm.app import PDBFile, Topology, element
from pdbfixer import PDBFixer

from openff.pablo._utils import unwrap, with_neighbours

PDBID = "1p3q"


@no_type_check
def energy_minimise(
    fixer: PDBFixer | PDBFile,
    forcefield: openmm.app.ForceField = openmm.app.ForceField(
        "amber14-all.xml",
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


def cap_to_charged_ter(fixer: PDBFixer) -> PDBFixer:
    old_pos = fixer.positions.value_in_unit(openmm.unit.angstrom)
    new_top = Topology()
    new_pos = []
    old_to_new = {}
    for chain in fixer.topology.chains():
        new_chain = new_top.addChain(chain.id)
        for prev_res, this_res, next_res in with_neighbours(chain.residues()):
            if this_res.name in {"ACE", "NME"}:
                continue
            new_res = new_top.addResidue(
                name=this_res.name,
                chain=new_chain,
                id=this_res.id,
                insertionCode=this_res.insertionCode,
            )
            for atom in this_res.atoms():
                old_to_new[atom] = new_top.addAtom(
                    name=atom.name,
                    element=atom.element,
                    residue=new_res,
                    id=atom.id,
                    formalCharge=atom.formalCharge,
                )
                new_pos.append(old_pos[atom.index])
            if next_res is not None and next_res.name == "NME":
                atom = unwrap(atom for atom in next_res.atoms() if atom.name == "N")
                old_to_new[atom] = new_top.addAtom(
                    "OXT",
                    element=element.oxygen,
                    residue=new_res,
                    id=atom.id,
                )
                new_pos.append(old_pos[atom.index])
            if prev_res is not None and prev_res.name == "ACE":
                atom = unwrap(atom for atom in prev_res.atoms() if atom.name == "C")
                old_to_new[atom] = new_top.addAtom(
                    "H2",
                    element=element.hydrogen,
                    residue=new_res,
                    id=atom.id,
                )
                new_pos.append(old_pos[atom.index])
                h3 = new_top.addAtom(
                    "H3",
                    element=element.hydrogen,
                    residue=new_res,
                )
                new_pos.append(old_pos[atom.index])
                new_top.addBond(
                    atom1=unwrap(atom for atom in new_res.atoms() if atom.name == "N"),
                    atom2=h3,
                )
                h1 = unwrap(atom for atom in new_res.atoms() if atom.name == "H")
                h1.name = "H1"

    for bond in fixer.topology.bonds():
        try:
            new_top.addBond(
                atom1=old_to_new[bond.atom1],
                atom2=old_to_new[bond.atom2],
                type=bond.type,
                order=bond.order,
            )
        except KeyError:
            pass

    fixer.topology = new_top

    fixer.positions = new_pos * openmm.unit.angstrom
    return fixer


def main():
    fixer = PDBFixer(pdbid=PDBID)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    last_resi = {
        chain.index: len([*chain.residues()]) for chain in fixer.topology.chains()
    }
    for (chain, res), missing in list(fixer.missingResidues.items()):
        if res == 0 or res == last_resi[chain]:
            del fixer.missingResidues[chain, res]
        elif len(missing) > 1:
            fixer.missingResidues[chain, res] = ["NME", "ACE"]
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    fixer = cap_to_charged_ter(fixer)

    energy_minimise(fixer)

    fnout = Path(f"{PDBID}_noter.pdb")

    with StringIO() as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
        lines = f.getvalue().splitlines()

    with open(fnout, "w") as f:
        for line in lines:
            if line.startswith("TER   "):
                continue
            f.write(line)
            f.write("\n")


if __name__ == "__main__":
    main()
