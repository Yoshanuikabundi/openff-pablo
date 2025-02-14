import math
from copy import deepcopy
from pathlib import Path
from collections.abc import Sequence

import openmm
import openmm.app
import openmm.unit
from openff.toolkit import ForceField, Quantity, Topology, unit

from openff.interchange import Interchange
from openff.interchange.components.potentials import Potential
from openff.interchange.exceptions import NonIntegralMoleculeChargeError
from openff.interchange.models import (
    LibraryChargeTopologyKey,
    PotentialKey,
    SingleAtomChargeTopologyKey,
)
from openff.nagl import GNNModel
from openff.nagl_models import validate_nagl_model_path
from openff.pablo._utils import draw_molecule

__all__ = [
    "draw_molecule",
]


def nglview_show_openmm(
    topology, positions: str | Path | Quantity, image_molecules=False
):
    import mdtraj
    import nglview
    import numpy as np
    from openff.units import ensure_quantity

    top = mdtraj.Topology.from_openmm(topology)

    if isinstance(positions, str) or isinstance(positions, Path):
        traj = mdtraj.load(positions, top=top)
        if image_molecules:
            traj.image_molecules(inplace=True)
    else:
        positions = ensure_quantity(positions, "openmm").value_in_unit(
            openmm.unit.nanometer
        )
        xyz = np.asarray([positions])
        box_vectors = topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            l1, l2, l3, alpha, beta, gamma = (
                mdtraj.utils.box_vectors_to_lengths_and_angles(
                    *np.asarray(box_vectors.value_in_unit(openmm.unit.nanometer))
                )
            )
            unitcell_angles, unitcell_lengths = [alpha, beta, gamma], [l1, l2, l3]
        else:
            unitcell_angles, unitcell_lengths = None, None
        traj = mdtraj.Trajectory(
            xyz, top, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles
        )
    widget = nglview.show_mdtraj(traj)
    widget.clear_representations()
    widget.add_cartoon()
    widget.add_line(opacity=0.5, crossSize=1.0)
    return widget


def get_charge_sum(
    interchange: Interchange,
    topology_indices: Sequence[int],
) -> Quantity:
    charges = {
        key.atom_indices[0]: value.m
        for key, value in interchange["Electrostatics"].charges.items()
    }

    return Quantity(
        sum([charges[index] for index in topology_indices]),
        "elementary_charge",
    )


def smear_charges(
    interchange: Interchange,
    topology_indices: Sequence[int],
) -> Interchange:
    total_formal_charge_of_topology_indices = sum(
        [interchange.topology.atom(i).formal_charge for i in topology_indices]
    )

    charge_to_smear = (
        get_charge_sum(interchange, topology_indices)
        - total_formal_charge_of_topology_indices
    )

    initial_charge_sum = get_charge_sum(interchange, topology_indices)

    per_atom_difference = charge_to_smear / len(topology_indices)

    interchange["Electrostatics"]._charges_cached = False

    for index in topology_indices:
        topology_key = SingleAtomChargeTopologyKey(this_atom_index=index)
        potential_key = interchange["Electrostatics"].key_map[topology_key]

        interchange["Electrostatics"].potentials[potential_key].parameters[
            "charge"
        ] -= per_atom_difference

    new_charge_sum = get_charge_sum(interchange, topology_indices)

    assert math.isclose((initial_charge_sum - new_charge_sum).m, charge_to_smear.m)

    assert math.isclose(
        get_charge_sum(interchange, topology_indices).m,
        total_formal_charge_of_topology_indices.m,
        abs_tol=1e-10,
        rel_tol=0,
    )

    return interchange


def get_total_charge(system: openmm.System) -> float:
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            return sum(
                [
                    force.getParticleParameters(index)[0]._value
                    for index in range(force.getNumParticles())
                ]
            )


def parametrize_with_nagl(
    force_field: ForceField,
    topology: Topology,
    nagl_method: str = "openff-gnn-am1bcc-0.1.0-rc.3.pt",
    allow_nonintegral_charges: bool = False,
) -> Interchange:
    """The protein has to be the first molecule"""
    print("adding dummy charges to force field ...")
    ff = deepcopy(force_field)
    # Add a dummy 0.0 library charge at the _top_ so it's only used as a last resort
    ff["LibraryCharges"].add_parameter(
        parameter_kwargs={
            "smirks": "[*:1]",
            "charge1": Quantity(0.0, "elementary_charge"),
            "name": "dummy",
        },
        before=0,  # "[#3+1:1]",
    )

    protein = topology.molecule(0)

    print("assigning graph charges ...")
    # protein.assign_partial_charges(
    #     partial_charge_method=nagl_method,
    #     toolkit_registry=NAGLToolkitWrapper(),
    # )

    nagl_path = validate_nagl_model_path(model=nagl_method)
    model = GNNModel.load(nagl_path, eval_mode=True, weights_only=False)
    charges = model.compute_property(
        protein,
        as_numpy=True,
        error_if_unsupported=True,
    )
    protein.partial_charges = Quantity(
        charges.astype(float),
        unit.elementary_charge,
    )

    print("making Interchange ...")
    interchange: Interchange = ff.create_interchange(
        topology,
        allow_nonintegral_charges=True,
    )

    potential_keys_to_remove = list()
    topology_keys_to_remove = list()

    nagl_indices = tuple(
        key.atom_indices[0]
        for key, val in interchange["Electrostatics"].key_map.items()
        if val.id == "[*:1]"
    )

    print("replacing dummy charges with NAGL charges ... ")
    for key, charge in interchange["Electrostatics"].charges.items():
        if key.atom_indices[0] in nagl_indices:
            # only modify charges where dummy placeholder of 0.0 was assigned from "[*:1]" parameter
            index = key.atom_indices[0]

            # must make new "single atom" topology key for each atom since the current
            # 1:many representation from the dummy charge is no longer valid
            new_potential_key = PotentialKey(
                id="inserted_graph_charges",
                associated_handler="molecules_with_preset_charges",
                mult=index,
            )
            # TODO: Compute graph charges as-needed here
            new_potential = Potential(
                parameters={"charge": protein.partial_charges[index]}
            )

            interchange["Electrostatics"].key_map[
                SingleAtomChargeTopologyKey(this_atom_index=index)
            ] = new_potential_key
            interchange["Electrostatics"].potentials.update(
                {new_potential_key: new_potential}
            )

            # remove the keys associated with the dummy library charge
            potential_keys_to_remove.append(
                interchange["Electrostatics"].key_map[
                    LibraryChargeTopologyKey(this_atom_index=index)
                ]
            )

            topology_keys_to_remove.append(
                LibraryChargeTopologyKey(this_atom_index=index)
            )

    for key_to_remove in topology_keys_to_remove:
        interchange["Electrostatics"].key_map.pop(key_to_remove)

    interchange["Electrostatics"]._charges_cached = False

    interchange = smear_charges(
        interchange=interchange,
        topology_indices=nagl_indices,
    )

    if not allow_nonintegral_charges:
        total_formal_charge = sum(
            atom.formal_charge for atom in interchange.topology.atoms
        ).m_as("elementary_charge")

        net_charge = get_charge_sum(
            interchange=interchange,
            topology_indices=range(interchange.topology.n_atoms),
        ).m_as("elementary_charge")

        if abs(total_formal_charge - net_charge) > 0.01:
            raise NonIntegralMoleculeChargeError(
                f"Interchange has a net charge of {net_charge} compared to a"
                + f" total formal charge  of {total_formal_charge}.",
            )

    return interchange
