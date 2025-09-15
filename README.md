# Shedding Hub [![Shedding Hub](https://github.com/shedding-hub/shedding-hub/actions/workflows/build.yaml/badge.svg)](https://github.com/shedding-hub/shedding-hub/actions/workflows/build.yaml) [![DOI](https://zenodo.org/badge/836912278.svg)](https://doi.org/10.5281/zenodo.15052772)

The Shedding Hub collates data and statistical models for biomarker shedding (such as viral RNA or drug metabolites) in different human specimen (such as stool or sputum samples). Developing wastewater-based epidemiology into a quantitative, reliable epidemiological monitoring tool motivates the project.

Datasets are extracted from appendices, figures, and supplementary materials of peer-reviewed studies. Each dataset is stored as a [`.yaml`](https://en.wikipedia.org/wiki/YAML) file and validated against our [data schema](data/.schema.yaml) to verify its integrity.

## 📊 Getting the Data

You can obtain the data by [downloading it from GitHub](https://github.com/shedding-hub/shedding-hub/tree/main/data). We also provide a [convenient Python package](http://pypi.org/project/shedding-hub/) so you can download the most recent data directly in your code or obtain a specific version of the data for reproducible analysis. Install the package by running `pip install shedding-hub` from the command line. The example below downloads the [data from Wölfel et al. (2020)](https://shedding-hub.github.io/datasets/woelfel2020virological.html) as of the commit [`259ca0d`](https://github.com/shedding-hub/shedding-hub/commit/259ca0d).

```python
>>> import shedding_hub as sh

>>> sh.load_dataset('woelfel2020virological', ref='259ca0d')
{'title': 'Virological assessment of hospitalized patients with COVID-2019',
 'doi': '10.1038/s41586-020-2196-x',
 ...}

```

## 🚀 API Access

The Shedding Hub now provides a REST API for programmatic access to datasets, built with FastAPI. The API offers comprehensive endpoints for browsing, searching, and retrieving biomarker shedding data.

### Starting the API Server

```bash
cd api
pip install -r requirements.txt
python run.py
```

The API will be available at `http://localhost:8004` with interactive documentation at `http://localhost:8004/docs`.

**Alternative methods:**
```bash
# Direct uvicorn command
uvicorn main:app --host 0.0.0.0 --port 8004

# Development mode with auto-reload
uvicorn main:app --host 0.0.0.0 --port 8004 --reload
```

### Key API Endpoints

| Endpoint | Description |
|----------|-------------|
| `GET /` | API information and available endpoints |
| `GET /datasets` | List all datasets with filtering options |
| `GET /datasets/{dataset_id}` | Get full dataset with measurements |
| `GET /datasets/{dataset_id}/metadata` | Get dataset metadata only |
| `GET /search` | Search datasets by keywords |
| `GET /stats` | Overall project statistics |
| `GET /health` | Health check endpoint |
| `GET /filters` | Available filter options |
| `POST /refresh-cache` | Refresh dataset cache |
| `GET /datasets/{dataset_id}/export` | Export dataset in CSV/JSON/Excel |

### Monitoring and Metrics

The API includes comprehensive monitoring capabilities:

| Endpoint | Description |
|----------|-------------|
| `GET /metrics` | Prometheus metrics for monitoring tools |
| `GET /admin/metrics` | JSON metrics summary for dashboards |
| `GET /admin/metrics/detailed` | Detailed performance metrics |

**Available Metrics:**
- Request count and duration by endpoint
- Cache hit/miss ratios
- System resource usage (CPU, memory)
- Error rates and rate limiting violations
- Dataset loading statistics

### Data Source Configuration

The API supports two data source modes:

**🌐 GitHub Integration (Default)**  
The API automatically fetches data directly from the live GitHub repository, ensuring you always have access to the latest datasets without manual synchronization.

**📁 Local Files**  
Traditional mode using local data files for air-gapped or offline deployments.

Toggle between modes using the `SHEDDING_HUB_GITHUB_ENABLED` setting:

```bash
# Use GitHub as data source (recommended)
SHEDDING_HUB_GITHUB_ENABLED=true

# Use local files instead
SHEDDING_HUB_GITHUB_ENABLED=false
SHEDDING_HUB_DATA_DIRECTORY=/path/to/local/data
```

**GitHub Integration Benefits:**
- ✅ Real-time access to the latest datasets
- ✅ No manual data synchronization required
- ✅ Automatic caching with configurable TTL
- ✅ Fallback to cached data if GitHub is temporarily unavailable
- ✅ Support for GitHub tokens to increase rate limits

### Configuration

The API supports extensive configuration via environment variables (prefix: `SHEDDING_HUB_`):

```bash
# Server Configuration
SHEDDING_HUB_HOST=0.0.0.0
SHEDDING_HUB_PORT=8000
SHEDDING_HUB_DEBUG=false

# Rate Limiting
SHEDDING_HUB_RATE_LIMIT_ENABLED=true
SHEDDING_HUB_RATE_LIMIT_DEFAULT=100/minute
SHEDDING_HUB_RATE_LIMIT_SEARCH=30/minute

# Cache Settings
SHEDDING_HUB_CACHE_TTL_SECONDS=3600

# GitHub Integration Settings
SHEDDING_HUB_GITHUB_ENABLED=true
SHEDDING_HUB_GITHUB_REPO_OWNER=shedding-hub
SHEDDING_HUB_GITHUB_REPO_NAME=shedding-hub
SHEDDING_HUB_GITHUB_BRANCH=main
SHEDDING_HUB_GITHUB_DATA_PATH=data
SHEDDING_HUB_GITHUB_TOKEN=your_github_token_here  # Optional, for higher rate limits
SHEDDING_HUB_GITHUB_CACHE_TTL=600  # GitHub cache TTL in seconds (10 minutes)

# Local Files (when GitHub disabled)
SHEDDING_HUB_DATA_DIRECTORY=/path/to/data
```

### Example Usage

```python
import requests

# Get all datasets
response = requests.get('http://localhost:8004/datasets?limit=10')
datasets = response.json()

# Search for COVID-related datasets
response = requests.get('http://localhost:8004/search?q=COVID&limit=5')
results = response.json()

# Get specific dataset
response = requests.get('http://localhost:8004/datasets/woelfel2020virological')
dataset = response.json()

# Export dataset as CSV
response = requests.get(
    'http://localhost:8004/datasets/woelfel2020virological/export?format=csv'
)

# Check data source information
response = requests.get('http://localhost:8004/admin/data-source')
source_info = response.json()
print(f"Data source: {source_info['source_type']}")
```

### Features

- **🌐 GitHub Integration**: Real-time data access from live repository with intelligent caching
- **🚀 High Performance**: Async processing with configurable caching strategies
- **🛡️ Rate Limiting**: Protection against abuse with configurable limits per endpoint
- **✅ Input Validation**: Comprehensive sanitization and Pydantic v2 validation
- **🔧 Error Handling**: Structured error responses with detailed logging
- **📊 Export Options**: CSV, JSON, and Excel format support
- **📈 Monitoring**: Prometheus-compatible metrics with system resource tracking
- **🌍 CORS Support**: Configurable cross-origin resource sharing
- **🗜️ Compression**: Automatic response compression for large datasets
- **🔄 Dual Data Sources**: Seamlessly switch between GitHub and local file modes

## 🤝 Contributing

Thank you for contributing your data to the Shedding Hub and supporting wastewater-based epidemiology! If you hit a bump along the road, [create a new issue](https://github.com/shedding-hub/shedding-hub/issues/new) and we'll sort it out together.

We use [pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) to add and update data, allowing for review and quality assurance. Learn more about the general workflow [here](https://docs.github.com/en/get-started/using-github/github-flow). To contribute your data, follow these easy steps (if you're already familiar with pull requests, steps 2 and 3 are for you):

1. Create a [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) of the Shedding Hub repository by clicking [here](https://github.com/shedding-hub/shedding-hub/fork) and [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the fork to your computer. You only have to do this once.
2. Create a new `my_cool_study/my_cool_study.yaml` file in the [`data`](data) directory and populate it with your data. See [here](data/woelfel2020virological/woelfel2020virological.yaml) for a comprehensive example from [Wölfel et al. (2020)](https://www.nature.com/articles/s41586-020-2196-x). A minimal example for studies with a single analyte (e.g., SARS-CoV-2 RNA concentration in stool samples) is available [here](tests/examples/valid_single_analyte.yaml), and a minimal example for studies with multiple analytes (e.g., crAssphage RNA concentration in stool samples and caffeine metabolites in urine) is available [here](tests/examples/valid_multiple_analytes.yaml).
3. Optionally, if you have a recent version of [Python](https://www.python.org) installed, you can validate your data to ensure it has the right structure before contributing it to the Shedding Hub.
    - Run `pip install -r requirements.txt` from the command line to install all the Python packages you need.
    - Run `pytest` from the command line to validate all datasets, including the one you just created.
4. Create a new [branch](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches) by running `git checkout -b my_cool_study`. Branches let you isolate changes you are making to the data, e.g., if you're simultaneously working on adding multiple studies–much appreciated! You should create a new branch from the `main` branch for each dataset you contribute; see [here](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) for more information.
5. Add your changes by running `git add data/my_cool_study/my_cool_study.yaml` and commit them by running `git commit -m "Add data from Someone et al. (20xx)."`. Feel free to pick another commit message if you prefer.
6. Push the dataset to your fork by running `git push origin my_cool_study`. This will send the data to GitHub, and the output of the command will include a line `Create a pull reuqest for 'my_cool_study' on GitHub by visiting: https://github.com/[your-username]/shedding-hub/pull/new/my_cool_study`. Click on the link and follow the next steps to create a new pull request.

Congratulations, you've just created your first pull request to contribute a new dataset! We'll now [review the changes](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/about-pull-request-reviews) you've made to make sure everything looks good. Once any questions have been resolved, we'll [merge your changes](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/merging-a-pull-request) into the repository. You've just contributed your first dataset to help make wastewater-based epidemiology a more quantitative public health monitoring tool–thank you!
