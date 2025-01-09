"""
Exceptions for the PDB loader.
"""

__all__ = [
    "NoMatchingResidueDefinitionError",
    "MultipleMatchingResidueDefinitionsError",
]


from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:
    from openff.pdbscan.pdb._pdb_data import PdbData, ResidueMatch


class NoMatchingResidueDefinitionError(ValueError):
    """Exception raised when a residue is missing from the database"""

    def __init__(self, res_atom_idcs: Sequence[int], data: "PdbData"):
        i = res_atom_idcs[0]
        super().__init__(
            "No residue definitions covered all atoms in residue"
            + f"{data.chain_id[i]}:{data.res_name[i]}#{data.res_seq[i]}"
        )


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
            + f"{data.chain_id[i]}:{data.res_name[i]}#{data.res_seq[i]}"
        )
