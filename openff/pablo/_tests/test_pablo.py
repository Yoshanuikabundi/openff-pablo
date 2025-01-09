"""
Unit and regression test for the pablo package.
"""

# Import package, test suite, and other packages as needed
import openff.pablo
import pytest
import sys


def test_pablo_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pablo" in sys.modules


def test_molecule_charges(example_molecule_with_charges):
    """Example test using a fixture defined in conftest.py"""
    assert example_molecule_with_charges.partial_charges is not None
