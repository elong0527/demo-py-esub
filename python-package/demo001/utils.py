"""
Utility functions for project management and common operations.
"""

from __future__ import annotations

from pathlib import Path
import polars as pl


def find_project_root(start_path: Path | None = None) -> Path:
    """
    Find the project root by looking for pyproject.toml file.

    Parameters:
    -----------
    start_path : Path, optional
        Starting path for search. Defaults to current working directory.

    Returns:
    --------
    Path
        Path to the project root directory

    Raises:
    -------
    FileNotFoundError
        If pyproject.toml is not found in any parent directory
    """
    if start_path is None:
        start_path = Path.cwd()

    current_path = start_path
    while current_path != current_path.parent:
        if (current_path / "pyproject.toml").exists():
            return current_path
        current_path = current_path.parent

    raise FileNotFoundError("Could not find pyproject.toml in any parent directory")


def load_adam_dataset(
    dataset_name: str, project_root: Path | None = None
) -> pl.DataFrame:
    """
    Load an ADaM dataset from the data directory.

    Parameters:
    -----------
    dataset_name : str
        Name of the dataset (e.g., 'adsl', 'adae', 'adlbc')
    project_root : Path, optional
        Project root path. If None, will be auto-detected.

    Returns:
    --------
    pl.DataFrame
        Loaded dataset as Polars DataFrame
    """
    if project_root is None:
        project_root = find_project_root()

    dataset_path = project_root / "data" / f"{dataset_name.lower()}.parquet"

    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset {dataset_path} not found")

    return pl.read_parquet(dataset_path)
