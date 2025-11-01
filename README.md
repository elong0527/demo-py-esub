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

## Project 

+ Run in batch `quarto render`
