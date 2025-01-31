"""
Tools for reading and patching the PDB Chemical Component Dictionary (CCD).

The ``CCD_RESIDUE_DEFINITION_CACHE`` loads, caches, and provides residue
definitions from the CCD and is the default residue database used by
``topology_from_pdb``.
"""

from pathlib import Path

from openff.pablo.residue import ResidueDefinition

from . import patches
from ._ccdcache import CcdCache
from .patches import (
    add_dephosphorylated_5p_terminus,
    add_disulfide_crosslink,
    add_protonation_variants,
    add_synonyms,
    disambiguate_alt_ids,
    fix_caps,
    patch_his_sidechain_zwitterion,
    set_hop3_leaving,
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
        {
            "U": add_dephosphorylated_5p_terminus,
            "G": add_dephosphorylated_5p_terminus,
            "C": add_dephosphorylated_5p_terminus,
            "A": add_dephosphorylated_5p_terminus,
            "DT": add_dephosphorylated_5p_terminus,
            "DG": add_dephosphorylated_5p_terminus,
            "DC": add_dephosphorylated_5p_terminus,
            "DA": add_dephosphorylated_5p_terminus,
        },
        {
            "U": set_hop3_leaving,
            "G": set_hop3_leaving,
            "C": set_hop3_leaving,
            "A": set_hop3_leaving,
            "DT": set_hop3_leaving,
            "DG": set_hop3_leaving,
            "DC": set_hop3_leaving,
            "DA": set_hop3_leaving,
        },
        {"*": disambiguate_alt_ids},
        {"*": add_synonyms},
        {"HIS": patch_his_sidechain_zwitterion},
    ],
    extra_definitions={"I": [ResidueDefinition.from_smiles("[I-:1]", {1: "I"}, "I")]},
)
"""The CCD, with commonly-required patches"""
