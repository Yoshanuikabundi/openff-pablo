from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Self
from urllib.request import urlopen

import xdg.BaseDirectory as xdg_base_dir

from ..residue import (
    ResidueDefinition,
    _skip_residue_definition_validation,
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
    extra_definitions
        Additional residue definitions to add to the cache. Note that patches
        are not applied to these definitions.
    """

    # TODO: Methods for adding entries from mapped SMILES

    def __init__(
        self,
        library_paths: Iterable[Path],
        cache_path: Path = Path(
            xdg_base_dir.save_cache_path("openff-pablo"),
            "ccd_cache",
        ),
        preload: list[str] = [],
        patches: Iterable[
            Mapping[
                str,
                Callable[[ResidueDefinition], list[ResidueDefinition]],
            ]
        ] = {},
        extra_definitions: Mapping[str, Iterable[ResidueDefinition]] = {},
    ):
        self._cache_path = cache_path.resolve()
        self._cache_path.mkdir(parents=True, exist_ok=True)

        self._library_paths = [path.resolve() for path in library_paths]

        self._definitions: dict[str, list[ResidueDefinition]] = {}
        self._patches: list[
            dict[
                str,
                Callable[[ResidueDefinition], list[ResidueDefinition]],
            ]
        ] = [dict(d) for d in patches]

        for path in self._glob("*.cif"):
            try:
                self._add_definition_from_str(path.read_text())
            except Exception:
                # If adding a file fails, skip it - we want an error at runtime, not importtime
                pass

        for resname in set(preload) - set(self._definitions):
            try:
                self[resname]
            except Exception:
                # If a preload fails, skip it - we want an error at runtime, not importtime
                pass

        self._extra_definitions_set = False
        for resname, resdefs in extra_definitions.items():
            self._extra_definitions_set = True
            self._add_definitions(resdefs, resname)

    def __repr__(self):
        return (
            f"CcdCache(path={self._cache_path},"
            + f" preload={list(self._definitions)},"
            + f" patches={self._patches!r}"
            + (", extra_definitions={...}" if self._extra_definitions_set else "")
            + ")"
        )

    def __getitem__(self, key: str) -> list[ResidueDefinition]:
        res_name = key.upper()
        if res_name in ["UNK", "UNL"]:
            # These residue names are reserved for unknown ligands/peptide residues
            raise KeyError(res_name)
        if res_name not in self._definitions:
            try:
                s = (self._cache_path / f"{res_name.upper()}.cif").read_text()
            except Exception:
                s = self._download_cif(res_name)

            self._add_definition_from_str(s, res_name=res_name)
        return self._definitions[res_name]

    def _apply_patches(
        self,
        residue_definition: ResidueDefinition,
    ) -> list[ResidueDefinition]:
        with _skip_residue_definition_validation():
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
        definition = ResidueDefinition.from_ccd_cif_str(s)
        self._add_and_patch_definition(definition, res_name)

    def _add_and_patch_definition(
        self,
        definition: ResidueDefinition,
        res_name: str | None = None,
    ) -> None:
        if res_name is None:
            res_name = definition.residue_name.upper()

        self._add_definitions(self._apply_patches(definition), res_name)

    def _add_definitions(
        self,
        definitions: Iterable[ResidueDefinition],
        res_name: str,
    ) -> None:
        for definition in definitions:
            if res_name != definition.residue_name.upper():
                raise ValueError(
                    f"ResidueDefinition {definition.residue_name}"
                    + f" ({definition.description}) must have residue name {res_name}",
                )
        self._definitions.setdefault(res_name, []).extend(definitions)

    def _download_cif(self, resname: str) -> str:
        with urlopen(
            f"https://files.rcsb.org/ligands/download/{resname.upper()}.cif",
        ) as stream:
            s: str = stream.read().decode("utf-8")
        path = self._cache_path / f"{resname.upper()}.cif"
        path.write_text(s)
        return s

    @property
    def _paths(self) -> Iterator[Path]:
        yield self._cache_path
        yield from self._library_paths

    def _glob(self, pattern: str) -> Iterator[Path]:
        for path in self._paths:
            yield from path.glob(pattern)

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

    def with_(
        self,
        extra_definitions: Mapping[str, Sequence[ResidueDefinition]],
    ) -> Self:
        new = deepcopy(self)
        for resname, resdefs in extra_definitions.items():
            new._extra_definitions_set = True
            new._add_definitions(resdefs, resname)
        return new
