"""
Exceptions for the PDB loader.
"""

__all__ = [
    "NoMatchingResidueDefinitionError",
    "MultipleMatchingResidueDefinitionsError",
]


from collections.abc import Sequence
from collections.abc import Iterable, Mapping

from openff.toolkit import Molecule

from openff.pablo._pdb_data import PdbData, ResidueMatch
from openff.pablo.residue import ResidueDefinition


class NoMatchingResidueDefinitionError(ValueError):
    """Exception raised when a residue is missing from the database"""

    def __init__(
        self,
        res_atom_idcs: Sequence[int],
        data: "PdbData",
        unknown_molecules: Iterable[Molecule],
        additional_substructures: Iterable[ResidueDefinition],
        residue_database: Mapping[
            str,
            Iterable[ResidueDefinition],
        ],
        verbose_errors: bool = False,
    ):
        i = res_atom_idcs[0]
        res_name = data.res_name[i]

        msg = [
            (
                "No residue definitions covered all atoms in residue"
                + f"{data.chain_id[i]}:{res_name}#{data.res_seq[i]}"
            ),
        ]
        if verbose_errors:
            residue_definitions = list(residue_database[res_name])
            if len(residue_definitions) == 0:
                msg.append("  No residues in database for residue name {res_name}")
            found_names = [data.name[i] for i in res_atom_idcs]
            for resdef in residue_definitions:
                resdef_names = [
                    "|".join([atom.name, *atom.synonyms]) for atom in resdef.atoms
                ]
                msg.append(f"    In {resdef.description}:")
                msg.append(f"      Expected {sorted(resdef_names)}")
                msg.append(f"      Found {sorted(found_names)}")
            # TODO: Describe residue_database and additional_substructures too

        super().__init__("\n".join(msg))


class MultipleMatchingResidueDefinitionsError(ValueError):
    def __init__(
        self,
        matches: Sequence["ResidueMatch"],
        res_atom_idcs: tuple[int, ...],
        data: "PdbData",
    ):
        i = res_atom_idcs[0]
        super().__init__(
            f"{len(matches)} residue definitions matched residue "
            + f"{data.chain_id[i]}:{data.res_name[i]}#{data.res_seq[i]}",
        )
