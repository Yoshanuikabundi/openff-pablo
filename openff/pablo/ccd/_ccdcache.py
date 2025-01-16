from collections.abc import Callable, Iterable, Iterator, Mapping
from io import StringIO
from pathlib import Path
from typing import no_type_check
from urllib.request import urlopen

from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader

from ..chem import PEPTIDE_BOND
from ..residue import (
    AtomDefinition,
    BondDefinition,
    ResidueDefinition,
    _defer_residue_definition_validation,
)

__all__ = [
    "CcdCache",
]


class CcdCache(Mapping[str, list[ResidueDefinition]]):
    """
    Caches, patches, and presents the CCD as a Python ``Mapping``.

    This requires internet access to work.

    Parameters
    ==========
    path
        The path to which to download CCD entries.
    preload
        A list of residue names to download when initializing the class.
    patches
        Functions to call on the given ``ResidueDefinitions`` before they are
        returned. A map from residue names to a single callable. The patch
        corresponding to key ``"*"`` will be applied to all residues before the
        more specific patches. Use :py:func:`combine_patches` to combine
        multiple patches into one.
    """

    # TODO: Methods for adding entries from mapped SMILES

    def __init__(
        self,
        path: Path,
        preload: list[str] = [],
        patches: Iterable[
            Mapping[
                str,
                Callable[[ResidueDefinition], list[ResidueDefinition]],
            ]
        ] = {},
    ):
        # TODO: Allow multiple cache paths so that some CIF files can be distributed with Pablo
        self._path = path.resolve()
        self._path.mkdir(parents=True, exist_ok=True)

        self._definitions: dict[str, list[ResidueDefinition]] = {}
        self._patches: list[
            dict[
                str,
                Callable[[ResidueDefinition], list[ResidueDefinition]],
            ]
        ] = [dict(d) for d in patches]

        for file in path.glob("*.cif"):
            self._add_definition_from_str(file.read_text())

        for resname in set(preload) - set(self._definitions):
            self[resname]

    def __repr__(self):
        return (
            f"CcdCache(path={self._path},"
            + f" preload={list(self._definitions)},"
            + f" patches={self._patches!r})"
        )

    def __getitem__(self, key: str) -> list[ResidueDefinition]:
        res_name = key.upper()
        if res_name in ["UNK", "UNL"]:
            # These residue names are reserved for unknown ligands/peptide residues
            raise KeyError(res_name)
        if res_name not in self._definitions:
            try:
                s = (self._path / f"{res_name.upper()}.cif").read_text()
            except Exception:
                s = self._download_cif(res_name)

            self._add_definition_from_str(s, res_name=res_name)
        return self._definitions[res_name]

    def _apply_patches(
        self,
        residue_definition: ResidueDefinition,
    ) -> list[ResidueDefinition]:
        with _defer_residue_definition_validation():
            definitions: list[ResidueDefinition] = [residue_definition]
            for patch_dict in self._patches:
                patched_definitions: list[ResidueDefinition] = []
                for definition in definitions:
                    star_patch = patch_dict.get("*", lambda x: [x])
                    res_patch = patch_dict.get(
                        residue_definition.residue_name.upper(),
                        lambda x: [x],
                    )
                    for star_patched_res in star_patch(definition):
                        for patched_res in res_patch(star_patched_res):
                            patched_definitions.append(patched_res)
                definitions = patched_definitions

        for definition in definitions:
            definition._validate()

        return definitions

    def _add_definition_from_str(self, s: str, res_name: str | None = None) -> None:
        definition = self._res_def_from_ccd_str(s)
        self._add_definition(definition, res_name)

    def _add_definition(
        self,
        definition: ResidueDefinition,
        res_name: str | None = None,
    ) -> None:
        if res_name is None:
            res_name = definition.residue_name.upper()

        patched_definitions = self._apply_patches(definition)

        assert all(
            res_name == definition.residue_name.upper()
            for definition in patched_definitions
        )

        self._definitions.setdefault(res_name, []).extend(patched_definitions)

    def _download_cif(self, resname: str) -> str:
        with urlopen(
            f"https://files.rcsb.org/ligands/download/{resname.upper()}.cif",
        ) as stream:
            s: str = stream.read().decode("utf-8")
        path = self._path / f"{resname.upper()}.cif"
        path.write_text(s)
        return s

    @staticmethod
    def _res_def_from_ccd_str(s: str) -> ResidueDefinition:
        @no_type_check
        def inner(s):
            # TODO: Handle residues like CL with a single atom properly (no tables)
            data = []
            with StringIO(s) as file:
                PdbxReader(file).read(data)
            block = data[0]

            parent_residue_name = (
                block.getObj("chem_comp").getValue("mon_nstd_parent_comp_id").upper()
            )
            parent_residue_name = (
                None if parent_residue_name == "?" else parent_residue_name
            )
            residueName = block.getObj("chem_comp").getValue("id").upper()
            residue_description = block.getObj("chem_comp").getValue("name")
            linking_type = block.getObj("chem_comp").getValue("type").upper()
            linking_bond = LINKING_TYPES[linking_type]

            atomData = block.getObj("chem_comp_atom")
            atomNameCol = atomData.getAttributeIndex("atom_id")
            altAtomNameCol = atomData.getAttributeIndex("alt_atom_id")
            symbolCol = atomData.getAttributeIndex("type_symbol")
            leavingCol = atomData.getAttributeIndex("pdbx_leaving_atom_flag")
            chargeCol = atomData.getAttributeIndex("charge")
            aromaticCol = atomData.getAttributeIndex("pdbx_aromatic_flag")
            stereoCol = atomData.getAttributeIndex("pdbx_stereo_config")

            atoms = [
                AtomDefinition(
                    name=row[atomNameCol],
                    synonyms=tuple(
                        [row[altAtomNameCol]]
                        if row[altAtomNameCol] != row[atomNameCol]
                        else [],
                    ),
                    symbol=row[symbolCol][0:1].upper() + row[symbolCol][1:].lower(),
                    leaving=row[leavingCol] == "Y",
                    charge=int(row[chargeCol]),
                    aromatic=row[aromaticCol] == "Y",
                    stereo=None if row[stereoCol] == "N" else row[stereoCol],
                )
                for row in atomData.getRowList()
            ]

            bondData = block.getObj("chem_comp_bond")
            if bondData is not None:
                atom1Col = bondData.getAttributeIndex("atom_id_1")
                atom2Col = bondData.getAttributeIndex("atom_id_2")
                orderCol = bondData.getAttributeIndex("value_order")
                aromaticCol = bondData.getAttributeIndex("pdbx_aromatic_flag")
                stereoCol = bondData.getAttributeIndex("pdbx_stereo_config")
                bonds = [
                    BondDefinition(
                        atom1=row[atom1Col],
                        atom2=row[atom2Col],
                        order={"SING": 1, "DOUB": 2, "TRIP": 3, "QUAD": 4}[
                            row[orderCol]
                        ],
                        aromatic=row[aromaticCol] == "Y",
                        stereo=None if row[stereoCol] == "N" else row[stereoCol],
                    )
                    for row in bondData.getRowList()
                ]
            else:
                bonds = []

            with _defer_residue_definition_validation():
                residue_definition = ResidueDefinition(
                    residue_name=residueName,
                    description=residue_description,
                    linking_bond=linking_bond,
                    crosslink=None,
                    atoms=tuple(atoms),
                    bonds=tuple(bonds),
                )

            return residue_definition

        return inner(s)

    def __contains__(self, value: object) -> bool:
        if value in self._definitions:
            return True
        if not isinstance(value, str):
            raise TypeError(
                f"CcdCache contains residue names of type str, not {type(value)}",
            )

        try:
            self[value]
        except Exception:
            # This catches KeyError, but also failures to download the residue
            return False
        else:
            return True

    def __iter__(self) -> Iterator[str]:
        return self._definitions.__iter__()

    def __len__(self) -> int:
        return self._definitions.__len__()


# TODO: Fill in this data
LINKING_TYPES: dict[str, BondDefinition | None] = {
    # "D-beta-peptide, C-gamma linking".upper(): [],
    # "D-gamma-peptide, C-delta linking".upper(): [],
    # "D-peptide COOH carboxy terminus".upper(): [],
    # "D-peptide NH3 amino terminus".upper(): [],
    # "D-peptide linking".upper(): [],
    # "D-saccharide".upper(): [],
    # "D-saccharide, alpha linking".upper(): [],
    # "D-saccharide, beta linking".upper(): [],
    # "DNA OH 3 prime terminus".upper(): [],
    # "DNA OH 5 prime terminus".upper(): [],
    # "DNA linking".upper(): [],
    # "L-DNA linking".upper(): [],
    # "L-RNA linking".upper(): [],
    # "L-beta-peptide, C-gamma linking".upper(): [],
    # "L-gamma-peptide, C-delta linking".upper(): [],
    # "L-peptide COOH carboxy terminus".upper(): [],
    # "L-peptide NH3 amino terminus".upper(): [],
    "L-peptide linking".upper(): PEPTIDE_BOND,
    # "L-saccharide".upper(): [],
    # "L-saccharide, alpha linking".upper(): [],
    # "L-saccharide, beta linking".upper(): [],
    # "RNA OH 3 prime terminus".upper(): [],
    # "RNA OH 5 prime terminus".upper(): [],
    # "RNA linking".upper(): [],
    "non-polymer".upper(): None,
    # "other".upper(): [],
    "peptide linking".upper(): PEPTIDE_BOND,
    "peptide-like".upper(): PEPTIDE_BOND,
    # "saccharide".upper(): [],
}
"""Map from the CCD's linking types to the bond formed between two such monomers"""
