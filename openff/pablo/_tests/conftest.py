"""
Global pytest fixtures
"""

# Use this file if you need to share any fixtures
# across multiple modules
# More information at
# https://docs.pytest.org/en/stable/how-to/fixtures.html#scope-sharing-fixtures-across-classes-modules-packages-or-session

import pytest

from openff.pablo._tests.data.files import EXAMPLE_SDF_WITH_CHARGES


@pytest.fixture
def example_molecule_with_charges():
    """Example fixture demonstrating how data files can be accessed"""
    from openff.toolkit import Molecule
    molecule = Molecule.from_file(EXAMPLE_SDF_WITH_CHARGES, file_format="sdf")
    return molecule
