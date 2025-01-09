"""
Location of data files for tests
================================

Use as ::

    from openff-pablo.tests.data.files import EXAMPLE_SDF_WITH_CHARGES
    from openff.toolkit import Molecule

    molecule = Molecule.from_file(EXAMPLE_SDF_WITH_CHARGES, file_format="sdf")

"""

__all__ = [
    "EXAMPLE_SDF_WITH_CHARGES",
]

from pkg_resources import resource_filename

EXAMPLE_SDF_WITH_CHARGES = resource_filename(__name__, "C1CC1.sdf")
