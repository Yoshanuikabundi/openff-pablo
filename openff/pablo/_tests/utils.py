import importlib.resources
from functools import lru_cache
from pathlib import Path


@lru_cache
def get_test_data_path(path: Path | str) -> Path:
    """Get the filename of a resource"""
    # TODO: Do this non-deprecated-ly
    assert __package__ is not None
    with importlib.resources.as_file(importlib.resources.files(__package__)) as root:
        pass
    assert root.is_dir()
    return root / "data" / path
