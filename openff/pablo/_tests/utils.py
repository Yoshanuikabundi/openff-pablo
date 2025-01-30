import importlib.resources
from functools import lru_cache
from pathlib import Path


@lru_cache
def get_test_data_path(path: Path | str) -> Path:
    """Get the filename of a resource"""
    # TODO: Do this non-deprecated-ly
    with importlib.resources.as_file(importlib.resources.files()) as root:
        pass
    assert root.is_dir()
    print(root)
    return root / "data" / path
