import dataclasses
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from functools import cached_property
from os import PathLike
from pathlib import Path
from typing import Any, DefaultDict, Self

from ._utils import __UNSET__, dec_hex, int_or_none, with_neighbours
from .exceptions import UnknownOrAmbiguousSerialInConectError
from .residue import AtomDefinition, ResidueDefinition


@dataclass(frozen=True)
class ResidueMatch:
    residue_definition: ResidueDefinition
    index_to_atomdef: dict[int, AtomDefinition]
    crosslink: tuple[int, int] | None = None
    """PDB indices of each bonded atom"""

    def atom(self, identifier: int | str) -> AtomDefinition:
        if isinstance(identifier, int):
            return self.index_to_atomdef[identifier]
        elif isinstance(identifier, str):
            return self.residue_definition.name_to_atom[identifier]
        else:
            raise TypeError(f"unknown identifier type {type(identifier)}")

    def set_crosslink(self, atom1_idx: int, atom2_idx: int) -> None:
        if (
            self.residue_definition.crosslink is None
            or atom1_idx not in self.res_atom_idcs
            or self.atom(atom1_idx).name != self.residue_definition.crosslink.atom1
            or atom2_idx in self.res_atom_idcs
        ):
            raise ValueError("bad crosslink index(es)")
        object.__setattr__(self, "crosslink", (atom1_idx, atom2_idx))

    @cached_property
    def res_atom_idcs(self) -> set[int]:
        return set(self.index_to_atomdef)

    @cached_property
    def prototype_index(self) -> int:
        return next(iter(self.index_to_atomdef))

    @cached_property
    def missing_atoms(self) -> set[str]:
        return {
            atom.name
            for atom in self.residue_definition.atoms
            if atom.name not in self.canonical_atom_name_to_index
        }

    @cached_property
    def missing_leaving_atoms(self) -> set[str]:
        return {
            atom_name
            for atom_name in self.missing_atoms
            if self.atom(atom_name).leaving
        }

    @cached_property
    def canonical_atom_name_to_index(self) -> dict[str, int]:
        return {atom.name: i for i, atom in self.index_to_atomdef.items()}

    @cached_property
    def expect_prior_bond(self) -> bool:
        if self.residue_definition.linking_bond is None:
            return False

        linking_atom = self.residue_definition.prior_bond_linking_atom
        expected_leaving_atoms = self.residue_definition.prior_bond_leaving_atoms

        return (
            linking_atom in self.canonical_atom_name_to_index
            and len(expected_leaving_atoms) > 0
            and expected_leaving_atoms.issubset(self.missing_leaving_atoms)
        )

    @cached_property
    def expect_posterior_bond(self) -> bool:
        if self.residue_definition.linking_bond is None:
            return False

        linking_atom = self.residue_definition.posterior_bond_linking_atom
        expected_leaving_atoms = self.residue_definition.posterior_bond_leaving_atoms

        return (
            linking_atom in self.canonical_atom_name_to_index
            and len(expected_leaving_atoms) > 0
            and expected_leaving_atoms.issubset(self.missing_leaving_atoms)
        )

    @cached_property
    def expect_crosslink(self) -> bool:
        if self.residue_definition.crosslink is None:
            return False

        linking_atom = self.residue_definition.crosslink.atom1
        expected_leaving_atoms = self.residue_definition.crosslink_leaving_atoms

        return (
            linking_atom in self.canonical_atom_name_to_index
            and len(expected_leaving_atoms) > 0
            and expected_leaving_atoms.issubset(self.missing_leaving_atoms)
        )

    def agrees_with(self, other: Self) -> bool:
        """True if both matches would assign the same chemistry, False otherwise"""
        if set(self.index_to_atomdef.keys()) != set(other.index_to_atomdef.keys()):
            return False

        name_map: dict[str, str] = {}
        for i, self_atom in self.index_to_atomdef.items():
            other_atom = other.index_to_atomdef[i]
            if not (
                self_atom.aromatic == other_atom.aromatic
                and self_atom.charge == other_atom.charge
                and self_atom.symbol == other_atom.symbol
                and self_atom.stereo == other_atom.stereo
            ):
                return False
            name_map[self_atom.name] = other_atom.name

        self_bonds = {
            (
                *sorted([name_map[bond.atom1], name_map[bond.atom2]]),
                bond.aromatic,
                bond.order,
                bond.stereo,
            )
            for bond in self.residue_definition.bonds
            if bond.atom1 in self.canonical_atom_name_to_index
            and bond.atom2 in self.canonical_atom_name_to_index
        }
        other_bonds = {
            (
                *sorted([bond.atom1, bond.atom2]),
                bond.aromatic,
                bond.order,
                bond.stereo,
            )
            for bond in other.residue_definition.bonds
            if bond.atom1 in self.canonical_atom_name_to_index
            and bond.atom2 in self.canonical_atom_name_to_index
        }
        if self_bonds != other_bonds:
            return False

        if self.expect_crosslink and (
            self.crosslink != other.crosslink
            or self.residue_definition.crosslink != other.residue_definition.crosslink
        ):
            return False

        if (self.expect_prior_bond or self.expect_posterior_bond) and (
            self.residue_definition.linking_bond
            != other.residue_definition.linking_bond
        ):
            return False

        return (
            self.expect_crosslink == other.expect_crosslink
            and self.expect_prior_bond == other.expect_prior_bond
            and self.expect_posterior_bond == other.expect_posterior_bond
        )


@dataclass
class PdbData:
    model: list[int | None] = field(default_factory=list)
    serial: list[int] = field(default_factory=list)
    name: list[str] = field(default_factory=list)
    alt_loc: list[str] = field(default_factory=list)
    res_name: list[str] = field(default_factory=list)
    chain_id: list[str] = field(default_factory=list)
    res_seq: list[int] = field(default_factory=list)
    i_code: list[str] = field(default_factory=list)
    x: list[float] = field(default_factory=list)
    y: list[float] = field(default_factory=list)
    z: list[float] = field(default_factory=list)
    occupancy: list[float] = field(default_factory=list)
    temp_factor: list[float] = field(default_factory=list)
    element: list[str] = field(default_factory=list)
    charge: list[int | None] = field(default_factory=list)
    terminated: list[bool] = field(default_factory=list)
    serial_to_index: DefaultDict[int, list[int]] = field(
        default_factory=lambda: defaultdict(list),
    )
    conects: list[set[int]] = field(default_factory=list)
    """The ith set contains atom indices CONECTed to atom index i"""
    cryst1_a: float | None = None
    cryst1_b: float | None = None
    cryst1_c: float | None = None
    cryst1_alpha: float | None = None
    cryst1_beta: float | None = None
    cryst1_gamma: float | None = None

    @classmethod
    def from_file(cls, path: str | PathLike[str]) -> Self:
        return cls.parse_pdb(Path(path).read_text().splitlines())

    def _append_coord_line(self, line: str):
        for field_ in dataclasses.fields(self):
            value = getattr(self, field_.name)
            if hasattr(value, "append"):
                value.append(__UNSET__)
                assert value[-1] is __UNSET__

        self.model[-1] = None
        self.serial[-1] = int(line[6:11])
        self.serial_to_index[self.serial[-1]].append(len(self.serial) - 1)
        self.name[-1] = line[12:16].strip()
        self.alt_loc[-1] = line[16].strip() or ""
        self.res_name[-1] = line[17:20].strip()
        self.chain_id[-1] = line[21].strip()
        self.res_seq[-1] = dec_hex(line[22:26])
        self.i_code[-1] = line[26].strip() or " "
        self.x[-1] = float(line[30:38])
        self.y[-1] = float(line[38:46])
        self.z[-1] = float(line[46:54])
        self.occupancy[-1] = float(line[54:60])
        self.temp_factor[-1] = float(line[60:66])
        self.element[-1] = line[76:78].strip()
        self.charge[-1] = int_or_none(line[78:80].strip())
        self.terminated[-1] = False
        self.conects[-1] = set()

        # Ensure we've assigned a value to every field
        for field_ in dataclasses.fields(self):
            value = getattr(self, field_.name)
            if hasattr(value, "append"):
                assert value[-1] is not __UNSET__

    @classmethod
    def parse_pdb(cls, lines: Iterable[str]) -> Self:
        model_n = None
        data = cls()
        for line in lines:
            if line.startswith("MODEL "):
                model_n = int(line[10:14])
            if line.startswith("ENDMDL "):
                model_n = None
            if line.startswith("HETATM") or line.startswith("ATOM  "):
                data._append_coord_line(line)
                data.model[-1] = model_n
            if line.startswith("TER   "):
                terminated_resname = line[17:20].strip() or data.res_name[-1]
                terminated_chainid = line[21].strip() or data.chain_id[-1]
                terminated_resseq = dec_hex(line[22:26]) or data.res_seq[-1]
                for i in range(-1, -99999, -1):
                    if (
                        data.res_name[i] == terminated_resname
                        and data.chain_id[i] == terminated_chainid
                        and data.res_seq[i] == terminated_resseq
                    ):
                        data.terminated[i] = True
                    else:
                        break
                else:
                    assert False, "last residue too big"
            if line.startswith("CRYST1"):
                data.cryst1_a = float(line[6:15])
                data.cryst1_b = float(line[15:24])
                data.cryst1_c = float(line[24:33])
                data.cryst1_alpha = float(line[33:40])
                data.cryst1_beta = float(line[40:47])
                data.cryst1_gamma = float(line[47:54])

        # Read all CONECT records
        data.conects = cls._process_conects(lines, data.serial_to_index, data.conects)

        return data

    @staticmethod
    def _process_conects(
        lines: Iterable[str],
        serial_to_index: dict[int, list[int]],
        conects: list[set[int]],
    ) -> list[set[int]]:
        for line in lines:
            if line.startswith("CONECT "):
                a = int(line[6:11])
                a_idcs = serial_to_index.get(a, [])
                if len(a_idcs) != 1:
                    raise UnknownOrAmbiguousSerialInConectError(a, a_idcs)
                a_idx = a_idcs[0]

                for start, stop in [(11, 16), (16, 21), (21, 26), (26, 31)]:
                    try:
                        b = int(line[start:stop])
                    except (ValueError, IndexError):
                        continue

                    b_idcs = serial_to_index.get(b, [])
                    if len(b_idcs) != 1:
                        raise UnknownOrAmbiguousSerialInConectError(b, b_idcs)
                    b_idx = b_idcs[0]

                    conects[a_idx].add(b_idx)
                    conects[b_idx].add(a_idx)
        return conects

    @property
    def residue_indices(self) -> Iterator[tuple[int, ...]]:
        indices = []
        prev = None
        for atom_idx, residue_info in enumerate(
            zip(
                self.model,
                self.res_name,
                self.chain_id,
                self.res_seq,
                self.i_code,
            ),
        ):
            if prev == residue_info or prev is None:
                indices.append(atom_idx)
            else:
                yield tuple(indices)
                indices = [atom_idx]
            prev = residue_info

        yield tuple(indices)

    def subset_matches_residue(
        self,
        res_atom_idcs: Sequence[int],
        residue_definition: ResidueDefinition,
    ) -> ResidueMatch | None:
        # Raise an error if the match would be empty - this way the
        # return value's truthiness always reflects whether there was a match
        if len(res_atom_idcs) == 0:
            raise ValueError("cannot match empty res_atom_idcs")

        # Skip definitions with too few atoms
        if len(residue_definition.atoms) < len(res_atom_idcs):
            return None

        # Skip non-(cross)linking definitions with the wrong number of atoms
        if (
            residue_definition.linking_bond is None
            and residue_definition.crosslink is None
            and len(
                residue_definition.atoms,
            )
            != len(res_atom_idcs)
        ):
            return None

        # Get the map from the canonical names to the indices
        try:
            index_to_atomdef = {
                i: residue_definition.name_to_atom[self.name[i]] for i in res_atom_idcs
            }
        except KeyError:
            return None

        matched_atoms = {atom.name for atom in index_to_atomdef.values()}

        # Fail to match if any atoms in PDB file got matched to more than one name
        if len(matched_atoms) != len(res_atom_idcs):
            return None

        # This assert should be guaranteed by the above
        assert set(index_to_atomdef.keys()) == set(res_atom_idcs)

        # Check that elements match, but tolerate missing columns and wrong case
        if any(
            self.element[i] != "" and self.element[i].lower() != atom.symbol.lower()
            for i, atom in index_to_atomdef.items()
        ):
            return None

        # Check that charges match, but tolerate missing columns
        if any(
            self.charge[i] is not None and self.charge[i] != atom.charge
            for i, atom in index_to_atomdef.items()
        ):
            return None

        missing_atoms = [
            atom for atom in residue_definition.atoms if atom.name not in matched_atoms
        ]

        # Match only if all the leaving atoms associated with each linking atom
        # is either entirely present or entirely absent
        missing_atom_names = {atom.name for atom in missing_atoms}
        if any(not atom.leaving for atom in missing_atoms):
            return None
        elif (
            (
                missing_atom_names.issuperset(
                    residue_definition.prior_bond_leaving_atoms,
                )
                or missing_atom_names.isdisjoint(
                    residue_definition.prior_bond_leaving_atoms,
                )
            )
            and (
                missing_atom_names.issuperset(
                    residue_definition.posterior_bond_leaving_atoms,
                )
                or missing_atom_names.isdisjoint(
                    residue_definition.posterior_bond_leaving_atoms,
                )
            )
            and (
                missing_atom_names.issuperset(
                    residue_definition.crosslink_leaving_atoms,
                )
                or missing_atom_names.isdisjoint(
                    residue_definition.crosslink_leaving_atoms,
                )
            )
        ):
            return ResidueMatch(
                index_to_atomdef=index_to_atomdef,
                residue_definition=residue_definition,
            )
        else:
            return None

    def get_residue_matches(
        self,
        residue_database: Mapping[str, Iterable[ResidueDefinition]],
        additional_substructures: Iterable[ResidueDefinition],
    ) -> Iterator[list[ResidueMatch]]:
        # Identify possible matches based on atom and residue names
        name_matches: list[list[ResidueMatch]] = []
        atom_idx_to_res_idx: dict[int, int] = {}
        for res_idx, res_atom_idcs in enumerate(self.residue_indices):
            prototype_index = res_atom_idcs[0]
            res_name = self.res_name[prototype_index]
            atom_idx_to_res_idx.update({i: res_idx for i in res_atom_idcs})

            residue_matches: list[ResidueMatch] = []
            for residue_definition in residue_database.get(res_name, []):
                match = self.subset_matches_residue(
                    res_atom_idcs,
                    residue_definition,
                )

                if match is not None:
                    residue_matches.append(match)

            if len(residue_matches) == 0:
                for residue_definition in additional_substructures:
                    match = self.subset_matches_residue(
                        res_atom_idcs,
                        residue_definition,
                    )
                    if match is not None:
                        residue_matches.append(match)

            name_matches.append(residue_matches)

        # Check for polymer bonds
        prev_filtered_matches: list[ResidueMatch] = []
        linkage_matches: list[list[ResidueMatch]] = []
        for _, this_matches, next_matches in with_neighbours(
            name_matches,
            default=(),
        ):
            neighbours_support_posterior_bond = any(
                next_match.expect_prior_bond for next_match in next_matches
            )
            neighbours_support_prior_bond = any(
                prev_match.expect_posterior_bond for prev_match in prev_filtered_matches
            )
            neighbours_support_molecule_end = (
                any(not next_match.expect_prior_bond for next_match in next_matches)
                or len(next_matches) == 0
            )
            neighbours_support_molecule_start = (
                any(
                    not prev_match.expect_posterior_bond
                    for prev_match in prev_filtered_matches
                )
                or len(prev_filtered_matches) == 0
            )
            this_filtered_matches: list[ResidueMatch] = []
            for match in this_matches:
                if len(match.missing_atoms) != 0:
                    prior_bond_mismatched = (
                        match.expect_prior_bond != neighbours_support_prior_bond
                    )
                    if prior_bond_mismatched:
                        continue

                    posterior_bond_mismatched = (
                        match.expect_posterior_bond != neighbours_support_posterior_bond
                    )
                    if posterior_bond_mismatched:
                        continue

                    this_filtered_matches.append(match)
                elif (
                    neighbours_support_molecule_end
                    and neighbours_support_molecule_start
                ):
                    this_filtered_matches.append(match)

            linkage_matches.append(this_filtered_matches)
            prev_filtered_matches = this_filtered_matches

        # Check for crosslinks
        # TODO: This could be simplified if we required crosslinking atoms not to have synonyms
        for residue_matches in linkage_matches:
            for match in residue_matches:
                if match.crosslink is not None:
                    # This match's crosslink has already been assigned
                    continue
                if not match.expect_crosslink:
                    continue
                this_crosslink_def = match.residue_definition.crosslink
                if this_crosslink_def is None:
                    # No crosslink defined for this match
                    continue
                this_crosslink_atom_idx = match.canonical_atom_name_to_index[
                    this_crosslink_def.atom1
                ]

                this_crosslink_conects = self.conects[this_crosslink_atom_idx]
                for other_crosslink_atom_idx in this_crosslink_conects:
                    other_crosslink_res_idx = atom_idx_to_res_idx[
                        other_crosslink_atom_idx
                    ]
                    other_matches = linkage_matches[other_crosslink_res_idx]
                    for other_match in other_matches:
                        other_crosslink_def = other_match.residue_definition.crosslink
                        other_crosslink_atom_canonical_name = other_match.atom(
                            other_crosslink_atom_idx,
                        ).name
                        if (
                            other_match.expect_crosslink
                            and other_crosslink_def is not None
                            and other_crosslink_def.flipped() == this_crosslink_def
                            and other_crosslink_def.atom2
                            == other_crosslink_atom_canonical_name
                        ):
                            # We've found a crosslink!
                            match.set_crosslink(
                                this_crosslink_atom_idx,
                                other_crosslink_atom_idx,
                            )
                            other_match.set_crosslink(
                                other_crosslink_atom_idx,
                                this_crosslink_atom_idx,
                            )
        # If there is a crosslink, we know we want it because it's in CONECTs
        # so filter out anything else
        for residue_matches in linkage_matches:
            if any(map(lambda x: x.crosslink is not None, residue_matches)):
                yielded_matches: list[ResidueMatch] = []
                for match in residue_matches:
                    if match.crosslink is not None:
                        yielded_matches.append(match)
                yield yielded_matches
            else:
                yield [match for match in residue_matches if not match.expect_crosslink]

    def __getitem__(self, index: int) -> dict[str, Any]:
        return {
            field.name: getattr(self, field.name)[index]
            for field in dataclasses.fields(self)
        }
