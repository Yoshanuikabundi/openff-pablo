"""
PDB loader that uses the CCD to read most PDB files without guessing bonds
"""

from . import ccd, exceptions, residue
from ._pdb import topology_from_pdb
from .ccd import CCD_RESIDUE_DEFINITION_CACHE

__all__ = [
    "topology_from_pdb",
    "CCD_RESIDUE_DEFINITION_CACHE",
    "exceptions",
    "ccd",
    "residue",
]

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
