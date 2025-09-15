# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

The Shedding Hub is a data repository and Python package for biomarker shedding datasets (viral RNA, drug metabolites) in human specimens. It's designed to support wastewater-based epidemiology research by providing standardized, validated datasets from peer-reviewed studies.

## Architecture

### Core Components

- **Data Storage**: Structured YAML files in `data/` following a strict schema (`data/.schema.yaml`)
- **Python Package**: Core utilities in `shedding_hub/` for data loading, processing, and analysis
- **REST API**: FastAPI-based web service in `api/` for programmatic access
- **Data Extraction**: Makefile-based system for processing extraction scripts
- **Validation**: pytest-based testing framework with schema validation

### Data Flow

1. Raw data extracted from publications via Python scripts (`*-extraction.py`) or Jupyter notebooks (`*-extraction.md`)
2. Data converted to standardized YAML format following the schema
3. Validation through pytest ensures data integrity
4. Python package provides programmatic access via `shedding_hub.load_dataset()`
5. REST API serves data over HTTP for web applications

## Common Commands

### Development and Testing
```bash
# Install dependencies
pip install -r requirements.txt

# Run all tests (validates all datasets)
pytest

# Run tests for specific functionality
pytest tests/test_util.py
pytest tests/test_data.py

# Data extraction (converts extraction scripts to YAML)
make extraction
```

### API Development
```bash
# Start API server
cd api/
python run.py
# OR
python main.py
# OR
uvicorn main:app --reload --host 0.0.0.0 --port 8000

# Test API endpoints
python test_api.py
```

### Dataset Management
```bash
# Backup existing data before extraction
make backup_data

# Verify data hasn't changed (for CI)
make assert_data_unchanged
```

## Data Schema and Validation

All datasets must conform to `data/.schema.yaml` which defines:
- Required fields (title, description, participants, measurements)
- Analyte specifications (biomarker, specimen type, units, detection limits)
- Participant attributes (demographics, clinical characteristics)
- Measurement structure (time series data with temporal offsets)

Key schema components:
- **Analytes**: Biomarkers being measured (SARS-CoV-2, PMMoV, crAssphage, etc.)
- **Specimens**: Sample types (stool, sputum, urine, swabs, etc.)
- **Reference Events**: Temporal anchors (symptom onset, enrollment, hospital admission)
- **Units**: Standardized measurement units (gc/mL, cycle threshold, etc.)

## Code Patterns

### Loading Datasets
Use `shedding_hub.load_dataset()` which supports:
- GitHub references: `load_dataset('woelfel2020virological', ref='commit_hash')`
- Pull requests: `load_dataset('dataset_name', pr=123)`
- Local files: `load_dataset('dataset_name', local='./data')`

### String Utilities
The package provides YAML string formatting utilities:
- `folded_str`: For long text with line wrapping
- `literal_str`: For preserving exact formatting
- `normalize_str`: For text cleaning with configurable options

### Shedding Analysis
Core analysis functions in `shedding_hub.shedding_duration`:
- `calc_shedding_duration()`: Individual participant analysis
- `calc_shedding_durations()`: Batch processing
- `plot_shedding_duration()`: Visualization utilities

## Important Files

- `data/.schema.yaml`: Authoritative data schema - all datasets must validate against this
- `Makefile`: Data extraction pipeline automation
- `requirements.txt`: Python dependencies for the core package
- `api/requirements.txt`: Additional dependencies for the REST API
- `pyproject.toml`: Package configuration and metadata

## Dataset Contribution Workflow

1. Create new directory: `data/study_name/`
2. Add extraction script: `study_name-extraction.py` or `study_name-extraction.md`
3. Run extraction: `make extraction` (generates `study_name.yaml`)
4. Validate: `pytest` (ensures schema compliance)
5. Use git workflow: branch → commit → pull request