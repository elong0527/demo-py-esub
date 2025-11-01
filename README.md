# demo-py-esub

## Overview

`demo-py-esub` is a demo project that illustrates how to organize analysis scripts using Python in a recommended folder structure for clinical study reports and regulatory submissions.

This demo project follows the concepts discussed in [Python for Clinical Study Reports and Submission](https://pycsr.org).

## Features

By using the recommended folder structure and associated development tools, this project demonstrates:

- **Consistency**: Standardized project organization
- **Automation**: Automated testing and code quality checks
- **Reproducibility**: Isolated environment management with uv
- **Compliance**: Best practices for regulatory submissions

## Folder Structure

This project follows a standardized folder structure for clinical study reports and regulatory submissions:

### Minimal Sufficient Folders and Files

- `pyproject.toml`: Python project metadata and dependencies
- `uv.lock`: Lock file for reproducible environments
- `_quarto.yml`: Quarto book configuration
- `index.qmd`: Book homepage
- `analysis/`: Analysis scripts using Quarto Markdown
- `data/`: ADaM datasets in Parquet format
- `output/`: Generated tables, listings, and figures (TLF)

### Additional Folders and Files

- `.venv/`: Virtual environment (created by uv)
- `.python-version`: Python version specification
- `python-package/demo001/`: DEMO-001 study-specific Python package
- `_book/`: Rendered Quarto book output
- `.gitignore`: Git ignore patterns

### Project Structure

```
demo-py-esub/
├── pyproject.toml           # Project configuration and dependencies
├── uv.lock                  # Dependency lock file
├── _quarto.yml             # Quarto book configuration
├── index.qmd               # Book homepage
├── analysis/               # Analysis scripts
│   ├── tlf-01-disposition.qmd
│   ├── tlf-02-population.qmd
│   ├── tlf-03-baseline.qmd
│   ├── tlf-04-efficacy-ancova.qmd
│   ├── tlf-05-ae-summary.qmd
│   └── tlf-06-specific.qmd
├── data/                   # ADaM datasets
│   ├── adsl.parquet
│   ├── adae.parquet
│   └── adlbc.parquet
├── output/                 # Generated TLF outputs
└── python-package/         # Clinical biostatistics package
    └── demo001/
        ├── __init__.py     # Package initialization
        ├── utils.py        # Utility functions
        ├── population.py   # Population analysis functions
        ├── baseline.py     # Baseline characteristics functions
        ├── safety.py       # Safety analysis functions
        └── efficacy.py     # Efficacy analysis functions
```

### Key Goals of This Structure

- **Consistency**: Standardized organization across projects
- **Automation**: Seamless integration with CI/CD pipelines
- **Reproducibility**: Isolated environments and locked dependencies
- **Compliance**: Follows regulatory submission best practices

### Benefits of This Approach

- **Code Reusability**: Functions can be used across multiple analyses
- **Standardization**: Consistent methodology across studies
- **Quality Assurance**: Centralized, tested functions reduce errors
- **Regulatory Compliance**: Standardized approaches for submissions
- **Maintainability**: Easy to update and improve functions

## Usage

### Install the Package
```bash
# Install in development mode
uv pip install -e .
```

### Run Analysis
```bash
# Render all analysis reports
quarto render

# Or render individual files
quarto render analysis/tlf-01-disposition.qmd
```

## Requirements

- Python 3.14.0
- [uv](https://docs.astral.sh/uv/) package manager

## Installation

+ Clone the repository:
   ```bash
   git https://github.com/elong0527/demo_py_esub
   cd demo-py-esub
   ```

+ Install uv (if not already installed):
   ```bash
   # Follow instructions at: https://pycsr.org/env-uv.html#installing-uv
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

+ Create and activate the project environment:
   ```bash
   uv sync
   ```

+ Activate the virtual environment:
   - **Windows**: `.venv\Scripts\activate`
   - **macOS/Linux**: `source .venv/bin/activate`