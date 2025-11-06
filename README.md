# demo-py-esub

## Overview

`demo-py-esub` is a project that illustrates how to organize analysis scripts
in a Python package for clinical study reports and regulatory submissions,
while following the packaging best practices such as the
[Python Packaging User Guide](https://packaging.python.org/).

This demo project follows the concepts discussed in
[Python for Clinical Study Reports and Submission](https://pycsr.org).

## Installation

Install uv (if not already installed):

```bash
# Follow instructions at: https://pycsr.org/env-uv.html#installing-uv
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Clone the repository:

```bash
git https://github.com/elong0527/demo-py-esub
cd demo-py-esub
```

Restore and activate the project environment:

```bash
uv sync
```

Activate the virtual environment:

- Windows: `.venv\Scripts\activate`
- macOS/Linux: `source .venv/bin/activate`

## Usage

### Install the package

```bash
# Install in development mode
uv pip install -e .
```

### Run analysis

```bash
# Render all analysis reports
quarto render

# Or render specific reports
quarto render analysis/tlf-01-disposition.qmd
```

## Folder structure

This project follows a standardized folder structure for clinical study reports
and regulatory submissions:

### Minimal sufficient folders and files

- `pyproject.toml`: Python project metadata and dependencies
- `src/demo001/`: Source code for the DEMO-001 study-specific package
- `tests/`: Unit tests for the package
- `data/`: ADaM datasets in Parquet format
- `output/`: Generated tables, listings, and figures (TLF)
- `_quarto.yml`: Quarto book configuration
- `index.qmd`: Quarto book landing page
- `analysis/`: Analysis scripts using Quarto Markdown

### Additional folders and files

- `.gitignore`: Ignore patterns for version control
- `.python-version`: Python version specification (created by uv)
- `uv.lock`: Lock file for reproducible environments (created by uv)
- `.venv/`: Virtual environment (created by uv)
- `_book/`: Rendered Quarto book output

### Project structure

```
demo-py-esub/
├── README.md               # Project overview
├── _quarto.yml             # Quarto project configuration
├── index.qmd               # Quarto book homepage
├── pyproject.toml          # Project configuration and dependencies
├── uv.lock                 # Dependency lock file
├── analysis/               # Analysis documents
│   ├── tlf-01-disposition.qmd
│   ├── tlf-02-population.qmd
│   ├── tlf-03-baseline.qmd
│   ├── tlf-04-efficacy-ancova.qmd
│   ├── tlf-05-ae-summary.qmd
│   └── tlf-06-specific.qmd
├── data/                   # ADaM datasets (parquet format)
│   ├── adae.parquet
│   ├── adlbc.parquet
│   ├── adlbh.parquet
│   ├── adlbhy.parquet
│   ├── adsl.parquet
│   ├── adtte.parquet
│   └── advs.parquet
├── output/                 # Generated TLF outputs
│   ├── tlf_disposition.rtf
│   ├── tlf_population.rtf
│   ├── tlf_baseline.rtf
│   ├── tlf_efficacy_ancova.rtf
│   ├── tlf_ae_summary.rtf
│   └── tlf_ae_specific.rtf
├── src/                    # Clinical biostatistics Python package
│   └── demo001/
│       ├── __init__.py     # Package initialization
│       ├── baseline.py     # Baseline characteristics functions
│       ├── efficacy.py     # Efficacy analysis functions
│       ├── population.py   # Population analysis functions
│       ├── safety.py       # Safety analysis functions
│       └── utils.py        # Shared utilities
└── tests/                  # Automated test suite
    ├── __init__.py
    └── test_utils.py
```

## Features

By using the recommended folder structure and associated development tools,
this project demonstrates:

- **Consistency**: Standardized project organization
- **Automation**: Automated testing and code quality checks
- **Reproducibility**: Isolated environment management with uv
- **Compliance**: Best practices for regulatory submissions

Benefits of this approach:

- **Code reusability**: Functions can be used across multiple analyses
- **Standardization**: Consistent methodology across studies
- **Quality assurance**: Centralized, tested functions reduce errors
- **Regulatory compliance**: Standardized approaches for submissions
- **Maintainability**: Easy to update and improve functions
