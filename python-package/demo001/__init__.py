"""
demo001: Python Package for DEMO-001 Study Analysis

A clinical biostatistics package for generating tables, listings, and figures
for regulatory submissions following ICH guidelines.
"""

from .utils import find_project_root, load_adam_dataset
from .population import count_by_treatment
from .baseline import summarize_continuous, summarize_categorical, get_value
from .safety import count_participants

__version__ = "0.1.0"
__author__ = "Clinical Biostatistics Team"

__all__ = [
    "find_project_root",
    "load_adam_dataset",
    "count_by_treatment",
    "summarize_continuous",
    "summarize_categorical",
    "get_value",
    "count_participants",
]
