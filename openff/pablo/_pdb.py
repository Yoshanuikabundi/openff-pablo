import itertools
from collections.abc import Mapping
from os import PathLike
from pathlib import Path
from typing import Iterable

import numpy as np
from openff.toolkit import Molecule, Topology
from openff.units import elements, unit
from typing_extensions import assert_never

from ._pdb_data import PdbData, ResidueMatch
from ._utils import (
    assign_stereochemistry_from_3d,
    cryst_to_box_vectors,
    sort_tuple,
)
from .ccd import CCD_RESIDUE_DEFINITION_CACHE
from .exceptions import (
    MultipleMatchingResidueDefinitionsError,
    NoMatchingResidueDefinitionError,
)
from .residue import ResidueDefinition

__all__ = [
    "topology_from_pdb",
]


def _load_unknown_residue(
    data: PdbData,
    indices: tuple[int, ...],
    unknown_molecules: Iterable[Molecule],
) -> Molecule | None:
    conects: set[tuple[int, int]] = set()
    serial_to_mol_index: dict[int, int] = {}
    pdbmol = Molecule()
    for pdb_index in indices:
        serial_to_mol_index[data.serial[pdb_index]] = pdbmol.add_atom(
            atomic_number=elements.NUMBERS[data.element[pdb_index]],
            formal_charge=data.charge[pdb_index],
            is_aromatic=False,
            stereochemistry=None,
            name=data.name[pdb_index],
            metadata={
                "residue_name": data.res_name[pdb_index],
                "leaving": False,
                "pdb_index": pdb_index,
                "residue_number": str(data.res_seq[pdb_index]),
                "res_seq": data.res_seq[pdb_index],
                "insertion_code": data.i_code[pdb_index],
                "chain_id": data.chain_id[pdb_index],
                "atom_serial": data.serial[pdb_index],
            },
        )
        for conect_serial in data.conects[pdb_index]:
            conects.add(sort_tuple((data.serial[pdb_index], conect_serial)))

    for serial1, serial2 in conects:
        pdbmol.add_bond(
            atom1=serial_to_mol_index[serial1],
            atom2=serial_to_mol_index[serial2],
            bond_order=1,
            is_aromatic=False,
        )

    for molecule in unknown_molecules:
        (match_found, mapping) = Molecule.are_isomorphic(
            molecule,
            pdbmol,
            return_atom_map=True,
            aromatic_matching=False,
            formal_charge_matching=False,
            bond_order_matching=False,
            atom_stereochemistry_matching=False,
            bond_stereochemistry_matching=False,
            strip_pyrimidal_n_atom_stereo=True,
        )
        if match_found:
            assert mapping is not None
            molecule = molecule.remap(mapping)
            for atom, pdbatom in zip(molecule.atoms, pdbmol.atoms):
                atom.metadata.update(pdbatom.metadata)
                atom.name = pdbatom.name
            molecule.generate_conformers(n_conformers=0, clear_existing=True)

            return molecule
    else:
        return None


def topology_from_pdb(
    path: PathLike[str],
    unknown_molecules: Iterable[Molecule] = [],
    residue_database: Mapping[
        str, Iterable[ResidueDefinition]
    ] = CCD_RESIDUE_DEFINITION_CACHE,
    additional_substructures: Iterable[ResidueDefinition] = [],
    use_canonical_names: bool = False,
    ignore_unknown_CONECT_records: bool = False,
    set_stereochemistry: bool = True,
) -> Topology:
    """
    Load a PDB file into an OpenFF ``Topology``.

    This function requires all hydrogens (and all other atoms) to be present in
    the PDB file, and that atom and residue names are consistent with the
    ``residue_database``. In return, it provides full chemical information on
    the entire PDB file.

    To load a PDB file with molecules including any residue not found in the
    CCD, or with residues that differ from that specified under a particular
    residue name, provide your own ``residue_database``. Any mapping from a
    residue name to a list of :py:data:`ResidueDefinition
    <openff.pdbscan.pdb.residue.ResidueDefinition>` objects may be used,
    but the :py:mod:`ccd <openff.pdbscan.pdb.ccd>` module  provides tools for
    augmenting the CCD.

    Alternatively, to load a single-residue molecule that is not present in the
    CCD, name that molecule ``"UNL"`` (or any name not present in the
    ``residue_database``), specify its CONECT records, and provide the
    appropriate molecule to the ``unknown_molecules`` argument.

    Parameters
    ----------
    path
        The path to the PDB file.
    unknown_molecules
        A list of molecules to match residues not found in the
        ``residue_database`` against. Unlike ``residue_database``, this requires
        that CONECT records be present and performs a match between the chemical
        graphs rather than using residue and atom names to detect chemistry.
    residue_database
        The database of residues to identify the atoms in the PDB file by. By
        default, a patched version of the CCD. Chemistry is identified by atom
        and residue names. If multiple residue definitions match a particular
        residue, the first one encountered is applied.
    additional_substructures
        Additional residue definitions to match against all residues that found
        no matches in the ``residue_database``. These definitions can match
        whether or not the residue name matches. To use this argument with
        OpenFF ``Molecule`` objects or SMILES strings, see the
        ``ResidueDefinition.from_*`` class methods.
    use_canonical_names
        If ``True``, atom names in the PDB file will be replaced by the
        canonical name for the same atom from the residue database.
    ignore_unknown_CONECT_records
        CONECT records do not include chemical information such as bond order
        and cannot be used on their own to add bonds beyond those specified
        through the residue database and unknown molecules. By default, any
        CONECT records not reflected in the final topology raise an error.
        If this argument is ``True``, this error is suppressed.
    set_stereochemistry
        If ``True``, stereochemistry will be set according to the structure of
        the PDB file. This takes considerable time. If ``False``, leave stereo
        unset.

    Notes
    -----

    This function uses a residue database to load a PDB file from its atom and
    residue names without guessing bonds. Bonds will be added by comparing atom
    and residue names to the residues defined in the ``residue_database``
    argument, which by default uses a patched version of the RCSB Chemical
    Component Dictionary (CCD). This is the dictionary of residue and atom names
    that the RCSB PDB is referenced against. The CCD is very large and cannot be
    distributed with this software, so by default internet access is required to
    use it.

    The following metadata are specified for all atoms produced by this function
    and can be accessed via ``topology.atom(i).metadata[key]``:

    ``"residue_name"``
        The residue name
    ``"residue_number"``
        The residue number as a string
    ``"res_seq"``
        The residue number as an integer
    ``"insertion_code"``
        The icode for the atom's residue. Used to align residue numbers between
        proteins with indels.
    ``"chain_id"``
        The letter identifier for the atom's chain.
    ``"pdb_index"``
        The atom's index in the PDB file. Sometimes called rank. Not to be
        confused with ``"atom_serial"``, which is the number given to the atom
        in the second column of the PDB file. Guaranteed to be unique and to
        match the index of the atom within the topology.
    ``"used_synonym"``
        The name of the atom that was found in the PDB file. By default,
        `atom.name` is set to this.
    ``"canonical_name"``
        The canonical name of the atom in the residue database. `atom.name` can
        be set to this with the `use_canonical_names` argument.
    ``"atom_serial"``
        The serial number of the atom, found in the second column of the PDB
        file. Not guaranteed to be unique.
    ``"matched_residue_description"``
        The residue description found in the residue database.
    ``"b_factor"``
        The temperature b-factor for the atom.
    ``"occupancy"``
        The occupancy for the atom.

    """
    # TODO: support streams and gzipped files
    path = Path(path)
    data = PdbData.parse_pdb(path.read_text().splitlines())

    molecules: list[Molecule] = []
    this_molecule = Molecule()
    prev_chain_id = data.chain_id[0]
    prev_model = data.model[0]
    for res_atom_idcs, matches in zip(
        data.residue_indices,
        data.get_residue_matches(residue_database, additional_substructures),
    ):
        # Check that we have a unique match, and error out or consult
        # unique_molecules as appropriate
        chemical_data: Molecule | ResidueMatch
        if len(matches) == 0:
            unknown_molecule = _load_unknown_residue(
                data, res_atom_idcs, unknown_molecules
            )
            if unknown_molecule is None:
                raise NoMatchingResidueDefinitionError(res_atom_idcs, data)
            else:
                chemical_data = unknown_molecule
        # If all matches would assign the same chemistry, accept it
        # assert all([]) == True
        elif all(a.agrees_with(b) for a, b in itertools.pairwise(matches)):
            chemical_data = matches[0]

            # this is a debug assert, if it triggers there's a bug
            assert set(res_atom_idcs) == chemical_data.res_atom_idcs
        else:
            raise MultipleMatchingResidueDefinitionsError(matches, res_atom_idcs, data)

        prototype_index = res_atom_idcs[0]

        # Terminate the previous molecule and start a new one if we can see that
        # this is the start of a new molecule
        if this_molecule.n_atoms > 0 and (
            data.chain_id[prototype_index] != prev_chain_id
            or data.model[prototype_index] != prev_model
            or (
                isinstance(chemical_data, ResidueMatch)
                and not chemical_data.expect_prior_bond
            )
            or (isinstance(chemical_data, Molecule))
        ):
            molecules.append(this_molecule)
            this_molecule = Molecule()

        # Apply the chemical data we've collected
        if isinstance(chemical_data, Molecule):
            this_molecule = chemical_data
        elif isinstance(chemical_data, ResidueMatch):
            add_to_molecule(
                this_molecule,
                res_atom_idcs,
                chemical_data,
                data,
                use_canonical_names,
            )
        else:
            assert_never(chemical_data)

        # Terminate the current molecule if we can see that this is the last residue
        if (
            data.terminated[prototype_index]
            or isinstance(chemical_data, Molecule)
            or (
                isinstance(chemical_data, ResidueMatch)
                and not chemical_data.expect_posterior_bond
            )
        ):
            molecules.append(this_molecule)
            this_molecule = Molecule()

        # TODO: Load other data from PDB file
        # TODO: Incorporate CONECT records
        # TODO: Deal with multi-model files

        prev_chain_id = data.chain_id[prototype_index]
        prev_model = data.model[prototype_index]
    molecules.append(this_molecule)

    for offmol in molecules:
        offmol._invalidate_cached_properties()
        offmol.add_default_hierarchy_schemes()

    topology = Topology.from_molecules(molecules)
    topology.set_positions(np.stack([data.x, data.y, data.z], axis=-1) * unit.angstrom)  # type: ignore

    if set_stereochemistry:
        for molecule in topology.molecules:
            # TODO: Speed this up
            #   - Build up molecules in RDMol form to skip conversion step?
            # This accounts for nearly half of the time to load 5ap1_prepared.pdb
            assign_stereochemistry_from_3d(molecule)

    if not ignore_unknown_CONECT_records:
        check_all_conects(topology, data)

    set_box_vectors(topology, data)

    return topology


def check_all_conects(topology: Topology, data: PdbData):
    all_bonds: set[tuple[int, int]] = {
        tuple(
            sorted([topology.atom_index(bond.atom1), topology.atom_index(bond.atom2)])
        )
        for bond in topology.bonds
    }  # type:ignore[assignment]
    index_of = {serial: i for i, serial in enumerate(data.serial)}

    conect_bonds: set[tuple[int, int]] = set()
    for i_idx, js in enumerate(data.conects):
        for j in js:
            j_idx = index_of[j]
            conect_bonds.add((i_idx, j_idx) if i_idx < j_idx else (j_idx, i_idx))
    if not conect_bonds.issubset(all_bonds):
        raise ValueError(
            "CONECT records without chemical information not supported",
            sorted(
                sorted([data.serial[a], data.serial[b]])
                for a, b in conect_bonds.difference(all_bonds)
            ),
        )


def set_box_vectors(topology: Topology, data: PdbData):
    if (
        data.cryst1_a is not None
        and data.cryst1_b is not None
        and data.cryst1_c is not None
        and data.cryst1_alpha is not None
        and data.cryst1_beta is not None
        and data.cryst1_gamma is not None
    ):
        topology.box_vectors = cryst_to_box_vectors(
            data.cryst1_a,
            data.cryst1_b,
            data.cryst1_c,
            data.cryst1_alpha,
            data.cryst1_beta,
            data.cryst1_gamma,
        )


def add_to_molecule(
    this_molecule: Molecule,
    res_atom_idcs: tuple[int, ...],
    residue_match: ResidueMatch,
    data: PdbData,
    use_canonical_names: bool,
) -> None:
    # Identify the previous linking atom
    linking_atom_idx: None | int = None
    if residue_match.expect_prior_bond:
        assert residue_match.residue_definition.linking_bond is not None
        linking_atom_name = residue_match.residue_definition.linking_bond.atom1
        for i in reversed(range(this_molecule.n_atoms)):
            if this_molecule.atom(i).metadata["canonical_name"] == linking_atom_name:
                linking_atom_idx = i
                break
        assert (
            linking_atom_idx is not None
        ), "Expecting a prior bond, but no linking atom found"

    # Add the residue to the current molecule
    atom_name_to_mol_idx: dict[str, int] = {}
    for pdb_index in res_atom_idcs:
        atom_def = residue_match.atom(pdb_index)

        if data.alt_loc[pdb_index] != "":
            # TODO: Support altlocs (probably in PdbData, maybe PdbData.residues()?)
            raise ValueError("altloc not yet supported")

        atom_name_to_mol_idx[atom_def.name] = this_molecule._add_atom(
            atomic_number=elements.NUMBERS[atom_def.symbol],
            formal_charge=atom_def.charge,
            is_aromatic=atom_def.aromatic,
            stereochemistry=None,
            name=atom_def.name if use_canonical_names else data.name[pdb_index],
            metadata={
                "residue_name": data.res_name[pdb_index],
                "res_seq": data.res_seq[pdb_index],
                "residue_number": str(data.res_seq[pdb_index]),
                "insertion_code": data.i_code[pdb_index],
                "chain_id": data.chain_id[pdb_index],
                "pdb_index": pdb_index,
                "used_synonym": data.name[pdb_index],
                "canonical_name": atom_def.name,
                "atom_serial": data.serial[pdb_index],
                "matched_residue_description": residue_match.residue_definition.description,
                "b_factor": str(data.temp_factor[pdb_index]),
                "occupancy": str(data.occupancy[pdb_index]),
            },
            invalidate_cache=False,
        )

    for bond in residue_match.residue_definition.bonds:
        if bond.atom1 in atom_name_to_mol_idx and bond.atom2 in atom_name_to_mol_idx:
            this_molecule._add_bond(
                atom1=atom_name_to_mol_idx[bond.atom1],
                atom2=atom_name_to_mol_idx[bond.atom2],
                bond_order=bond.order,
                is_aromatic=bond.aromatic,
                stereochemistry=None,
                invalidate_cache=False,
            )

    if linking_atom_idx is not None:
        linking_bond = residue_match.residue_definition.linking_bond
        assert (
            linking_bond is not None
        ), "linking_atom_idx is only set when linking_atom_idx is None"
        this_molecule._add_bond(
            atom1=linking_atom_idx,
            atom2=atom_name_to_mol_idx[linking_bond.atom2],
            bond_order=linking_bond.order,
            is_aromatic=linking_bond.aromatic,
            stereochemistry=None,
            invalidate_cache=False,
        )
