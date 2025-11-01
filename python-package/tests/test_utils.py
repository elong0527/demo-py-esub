"""
Tests for demo001.utils module
"""

import pytest
from pathlib import Path
import polars as pl
from demo001.utils import find_project_root, load_adam_dataset


def test_find_project_root():
    """Test finding project root"""
    root = find_project_root()
    assert (root / "pyproject.toml").exists()
    assert root.name == "demo-py-esub"


def test_load_adam_dataset():
    """Test loading ADaM datasets"""
    root = find_project_root()

    # Test loading ADSL
    adsl = load_adam_dataset("adsl", root)
    assert isinstance(adsl, pl.DataFrame)
    assert adsl.height > 0

    # Test that USUBJID column exists
    assert "USUBJID" in adsl.columns


def test_load_nonexistent_dataset():
    """Test loading non-existent dataset raises error"""
    root = find_project_root()

    with pytest.raises(FileNotFoundError):
        load_adam_dataset("nonexistent", root)
