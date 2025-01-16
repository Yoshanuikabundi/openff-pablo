"""
Tools for reading and patching the PDB Chemical Component Dictionary (CCD).

The ``CCD_RESIDUE_DEFINITION_CACHE`` loads, caches, and provides residue
definitions from the CCD and is the default residue database used by
``topology_from_pdb``.
"""

from pathlib import Path

from . import patches
from ._ccdcache import CcdCache
from .patches import (
    add_disulfide_crosslink,
    add_protonation_variants,
    add_synonyms,
    disambiguate_alt_ids,
    fix_caps,
)

__all__ = [
    "CCD_RESIDUE_DEFINITION_CACHE",
    "CcdCache",
    "patches",
]

# TODO: Replace these patches with CONECT records?
CCD_RESIDUE_DEFINITION_CACHE: CcdCache = CcdCache(
    # TODO: Store the user's cache in a more appropriate location
    Path(__file__).parent / "../../../.ccd_cache",
    patches=[
        {
            "ACE": fix_caps,
            "NME": fix_caps,
            "CYS": add_disulfide_crosslink,
        },
        {"*": add_protonation_variants},
        {"*": disambiguate_alt_ids},
        {"*": add_synonyms},
    ],
)
"""The CCD, with commonly-required patches"""
