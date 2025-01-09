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

versions = get_versions()  # type: ignore[unknownType]
__version__: str = versions["version"]  # type: ignore[assignment]
__git_revision__: str = versions["full-revisionid"]  # type: ignore[assignment]
del get_versions, versions
