import dataclasses
from dataclasses import dataclass, field
from functools import cached_property
from symtable import Symbol
from typing import Any, Iterable, Iterator, Mapping, Self, Sequence

from ._utils import __UNSET__, dec_hex, with_neighbours
from .residue import AtomDefinition, ResidueDefinition


@dataclass(frozen=True)
class ResidueMatch:
    index_to_atomdef: dict[int, AtomDefinition]
    residue_definition: ResidueDefinition
    missing_atoms: set[str]

    def atom(self, identifier: int | str) -> AtomDefinition:
        if isinstance(identifier, int):
            return self.index_to_atomdef[identifier]
        elif isinstance(identifier, str):
            return self.residue_definition.name_to_atom[identifier]
        else:
            raise TypeError(f"unknown identifier type {type(identifier)}")

    @cached_property
    def res_atom_idcs(self) -> set[int]:
        return set(self.index_to_atomdef)

    @cached_property
    def prototype_index(self) -> int:
        return next(iter(self.index_to_atomdef))

    @cached_property
    def missing_leaving_atoms(self) -> set[str]:
        return {
            atom_name
            for atom_name in self.missing_atoms
            if self.atom(atom_name).leaving
        }

    @cached_property
    def matched_canonical_atom_names(self) -> set[str]:
        return {atom.name for atom in self.index_to_atomdef.values()}

    @cached_property
    def expect_prior_bond(self) -> bool:
        if self.residue_definition.linking_bond is None:
            return False

        linking_atom = self.residue_definition.prior_bond_linking_atom
        expected_leaving_atoms = self.residue_definition.prior_bond_leaving_atoms

        return (
            linking_atom in self.matched_canonical_atom_names
            and len(expected_leaving_atoms) > 0
            and self.missing_leaving_atoms.intersection(expected_leaving_atoms)
            == expected_leaving_atoms
        )

    @cached_property
    def expect_posterior_bond(self) -> bool:
        if self.residue_definition.linking_bond is None:
            return False

        linking_atom = self.residue_definition.posterior_bond_linking_atom
        expected_leaving_atoms = self.residue_definition.posterior_bond_leaving_atoms

        return (
            linking_atom in self.matched_canonical_atom_names
            and len(expected_leaving_atoms) > 0
            and self.missing_leaving_atoms.intersection(expected_leaving_atoms)
            == expected_leaving_atoms
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
            if bond.atom1 in self.matched_canonical_atom_names
            and bond.atom2 in self.matched_canonical_atom_names
        }
        other_bonds = {
            (
                *sorted([bond.atom1, bond.atom2]),
                bond.aromatic,
                bond.order,
                bond.stereo,
            )
            for bond in other.residue_definition.bonds
            if bond.atom1 in self.matched_canonical_atom_names
            and bond.atom2 in self.matched_canonical_atom_names
        }
        return self_bonds == other_bonds


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
    charge: list[int] = field(default_factory=list)
    terminated: list[bool] = field(default_factory=list)
    conects: list[set[int]] = field(default_factory=list)
    cryst1_a: float | None = None
    cryst1_b: float | None = None
    cryst1_c: float | None = None
    cryst1_alpha: float | None = None
    cryst1_beta: float | None = None
    cryst1_gamma: float | None = None

    def _append_coord_line(self, line: str):
        for field_ in dataclasses.fields(self):
            value = getattr(self, field_.name)
            if hasattr(value, "append"):
                value.append(__UNSET__)
                assert value[-1] is __UNSET__

        self.model[-1] = None
        self.serial[-1] = int(line[6:11])
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
        self.charge[-1] = int(line[78:80].strip() or 0)
        self.terminated[-1] = False
        self.conects[-1] = set()

        # Ensure we've assigned a value to every field
        for field_ in dataclasses.fields(self):
            value = getattr(self, field_.name)
            if hasattr(value, "append"):
                assert value[-1] is not __UNSET__

    @classmethod
    def parse_pdb(cls, lines: Iterable[str]) -> Self:
        conects = {}
        # Read all CONECT records
        for line in lines:
            if line.startswith("CONECT "):
                a = int(line[6:11])
                bs = []
                for start, stop in [(11, 16), (16, 21), (21, 26), (26, 31)]:
                    try:
                        b = int(line[start:stop])
                    except (ValueError, IndexError):
                        continue
                    bs.append(b)
                    conects.setdefault(b, set()).add(a)
                conects.setdefault(a, set()).update(bs)

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
                data.conects[-1].update(conects.get(data.serial[-1], []))
            if line.startswith("TER   "):
                terminated_resname = line[17:20].strip() or data.res_name[-1]
                terminated_chainid = line[21].strip() or data.chain_id[-1]
                terminated_resseq = dec_hex(line[22:26]) or data.res_seq[-1]
                for i in range(-1, -999, -1):
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

        return data

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
            )
        ):
            if prev == residue_info or prev is None:
                indices.append(atom_idx)
            else:
                yield tuple(indices)
                indices = [atom_idx]
            prev = residue_info

        yield tuple(indices)

    def get_residue_matches(
        self,
        residue_database: Mapping[str, Iterable[ResidueDefinition]],
        additional_substructures: Iterable[ResidueDefinition],
    ) -> Iterator[list[ResidueMatch]]:
        all_matches: list[tuple[ResidueMatch, ...]] = []
        for res_atom_idcs in self.residue_indices:
            prototype_index = res_atom_idcs[0]
            res_name = self.res_name[prototype_index]

            matches: list[ResidueMatch] = []
            for residue_definition in residue_database.get(res_name, []):
                match = self.subset_matches_residue(
                    res_atom_idcs,
                    residue_definition,
                )

                if match is not None:
                    matches.append(match)

            if len(matches) == 0:
                for residue_definition in additional_substructures:
                    match = self.subset_matches_residue(
                        res_atom_idcs, residue_definition
                    )
                    if match is not None:
                        matches.append(match)

            all_matches.append(tuple(matches))

        prev_filtered_matches: list[ResidueMatch] = []
        for _, this_matches, next_matches in with_neighbours(
            all_matches,
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

                    # assert any([]) == False
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

            yield this_filtered_matches

            prev_filtered_matches = this_filtered_matches

    def __getitem__(self, index: int) -> dict[str, Any]:
        return {
            field.name: getattr(self, field.name)[index]
            for field in dataclasses.fields(self)
        }

    def subset_matches_residue(
        self,
        res_atom_idcs: Sequence[int],
        residue_definition: ResidueDefinition,
    ) -> ResidueMatch | None:
        # Raise an error if the returned dict would be empty - this way the
        # return value's truthiness always reflects whether there was a match
        if len(res_atom_idcs) == 0:
            raise ValueError("cannot match empty res_atom_idcs")

        # Skip definitions with too few atoms
        if len(residue_definition.atoms) < len(res_atom_idcs):
            return None

        # Skip non-linking definitions with the wrong number of atoms
        if residue_definition.linking_bond is None and len(
            residue_definition.atoms
        ) != len(res_atom_idcs):
            return None

        # Get the map from the canonical names to the indices
        try:
            index_to_atomdef = {
                i: residue_definition.name_to_atom[self.name[i]] for i in res_atom_idcs
            }
        except KeyError as e:
            return None

        matched_atoms = set(atom.name for atom in index_to_atomdef.values())

        # Fail to match if any atoms in PDB file got matched to more than one name
        if len(matched_atoms) != len(res_atom_idcs):
            return None

        # This assert should be guaranteed by the above
        assert set(index_to_atomdef.keys()) == set(res_atom_idcs)

        missing_atoms = [
            atom for atom in residue_definition.atoms if atom.name not in matched_atoms
        ]

        # Match only if the set of all missing atoms is one of the following:
        # - empty
        # - the prior bond leaving fragment
        # - the posterior bond leaving fragment
        # - both leaving fragments
        if any(not atom.leaving for atom in missing_atoms):
            return None
        elif set(atom.name for atom in missing_atoms) in [
            set(),
            residue_definition.prior_bond_leaving_atoms.union(
                residue_definition.posterior_bond_leaving_atoms
            ),
            residue_definition.prior_bond_leaving_atoms,
            residue_definition.posterior_bond_leaving_atoms,
        ]:
            return ResidueMatch(
                index_to_atomdef=index_to_atomdef,
                residue_definition=residue_definition,
                missing_atoms={atom.name for atom in missing_atoms},
            )
        else:
            return None

    def are_alt_locs(self, i: int, j: int) -> bool:
        if i == j:
            raise ValueError(f"i and j are the same ({i})")
        if max(i, j) - min(i, j) == 1:
            return (
                self.model[i],
                self.name[i],
                self.res_name[i],
                self.chain_id[i],
                self.res_seq[i],
                self.i_code[i],
            ) == (
                self.model[j],
                self.name[j],
                self.res_name[j],
                self.chain_id[j],
                self.res_seq[j],
                self.i_code[j],
            )
        else:
            return self.are_alt_locs(i, i + 1) and self.are_alt_locs(i + 1, j)
