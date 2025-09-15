# Shedding Hub API

A production-ready FastAPI-based REST API for accessing biomarker shedding datasets from the Shedding Hub project.

## Features

- **Type-Safe API**: Full Pydantic model validation with automatic OpenAPI documentation
- **Rate Limiting**: IP-based rate limiting with configurable limits per endpoint type
- **Smart Caching**: In-memory caching with configurable TTL for optimal performance
- **Environment Configuration**: Flexible settings management with environment variables
- **Comprehensive Testing**: Full pytest test suite with async support and mocking
- **List datasets**: Get all available datasets with optional filtering
- **Get specific datasets**: Retrieve complete dataset information  
- **Search**: Text-based search across datasets with relevance scoring
- **Metadata only**: Get dataset metadata without measurements for faster responses
- **Statistics**: Overall project statistics and analytics
- **Filtering**: Filter by biomarker, specimen type, etc. with enum validation
- **Health checks**: Monitor API status with detailed health information

## Quick Start

### 1. Install Dependencies

```bash
cd api/
pip install -r requirements.txt
```

### 2. Configure Environment (Optional)

```bash
# Copy example environment file
cp .env.example .env

# Edit configuration as needed
nano .env
```

### 3. Start the API Server

```bash
# Option 1: Using the convenience script (recommended for development)
python run.py

# Option 2: Using uvicorn directly (recommended for production)
uvicorn main:app --reload --host 0.0.0.0 --port 8000

# Option 3: Direct Python execution
python main.py
```

### 4. Access the API

- **API Server**: http://localhost:8000
- **Interactive Documentation**: http://localhost:8000/docs
- **ReDoc Documentation**: http://localhost:8000/redoc
- **OpenAPI Schema**: http://localhost:8000/openapi.json

### 5. Test the API

```bash
# Run the comprehensive test suite
pytest test_api.py -v

# Run with coverage
pytest test_api.py --cov=main --cov-report=html

# Test rate limiting specifically
python test_rate_limiting.py
```

## API Endpoints

### Core Endpoints

- `GET /` - API information and available endpoints
- `GET /health` - Health check endpoint with system status
- `GET /datasets` - List all datasets with optional filtering and validation
- `GET /datasets/{dataset_id}` - Get specific dataset with optional measurement inclusion
- `GET /datasets/{dataset_id}/metadata` - Get dataset metadata only (faster)
- `GET /datasets/{dataset_id}/export` - Export dataset in CSV, JSON, or Excel format
- `GET /search` - Search datasets by text query with relevance scoring
- `GET /stats` - Overall project statistics and cache information
- `GET /filters` - Get available filter options (biomarkers, specimens)
- `POST /refresh-cache` - Refresh the dataset cache (rate limited)
- `GET /metrics` - Prometheus metrics endpoint for monitoring
- `GET /admin/metrics` - API metrics summary for dashboards
- `GET /admin/metrics/detailed` - Detailed metrics for admin monitoring

### Rate Limits

- **Default endpoints**: 100 requests/minute
- **Search endpoints**: 30 requests/minute
- **Expensive operations** (full dataset retrieval, exports, admin metrics): 10 requests/minute
- **Cache operations**: 5 requests/minute

### Query Parameters

#### List Datasets (`/datasets`)
- `biomarker` (optional): Filter by biomarker enum (`SARS-CoV-2`, `PMMoV`, `crAssphage`, etc.)
- `specimen` (optional): Filter by specimen type enum (`stool`, `sputum`, `urine`, etc.)
- `limit` (default: 100, max: 1000): Maximum number of datasets to return

#### Get Dataset (`/datasets/{dataset_id}`)
- `include_measurements` (default: true): Include participant measurements

#### Export Dataset (`/datasets/{dataset_id}/export`)
- `format` (required): Export format - `csv`, `json`, or `excel`
- `include_measurements` (default: true): Include participant measurements

#### Search (`/search`)
- `q` (required, min_length: 1): Search query string
- `fields` (default: ["title", "description"]): Fields to search in
- `limit` (default: 20, max: 1000): Maximum number of results

## Type Safety & Validation

The API uses comprehensive Pydantic models for request/response validation:

### Validated Enums
- **Biomarkers**: `SARS-CoV-2`, `mtDNA`, `PMMoV`, `crAssphage`, `influenza`, `sapovirus`, `SARS`
- **Specimens**: `stool`, `sputum`, `urine`, `plasma`, `nasopharyngeal_swab`, `anterior_nares_swab`, `saliva`, etc.
- **Units**: `gc/mL`, `gc/dry gram`, `gc/swab`, `gc/wet gram`, `cycle threshold`, `pfu/mL`
- **Reference Events**: `symptom onset`, `confirmation date`, `enrollment`, `hospital admission`
- **Export Formats**: `csv`, `json`, `excel`

### Validation Examples
```bash
# ✅ Valid request
curl "http://localhost:8000/datasets?biomarker=SARS-CoV-2&limit=10"

# ❌ Invalid biomarker (422 error)
curl "http://localhost:8000/datasets?biomarker=invalid"

# ❌ Invalid limit (422 error)  
curl "http://localhost:8000/datasets?limit=99999"
```

## Example Usage

### Python with requests

```python
import requests

# List all datasets with validation
response = requests.get("http://localhost:8000/datasets")
datasets = response.json()

# Get a specific dataset
response = requests.get("http://localhost:8000/datasets/woelfel2020virological")
dataset = response.json()

# Search for COVID-related datasets
response = requests.get("http://localhost:8000/search", params={"q": "COVID", "limit": 5})
results = response.json()

# Filter by validated biomarker enum
response = requests.get("http://localhost:8000/datasets", params={"biomarker": "SARS-CoV-2"})
covid_datasets = response.json()

# Get project statistics
response = requests.get("http://localhost:8000/stats")
stats = response.json()

# Export dataset to CSV
response = requests.get("http://localhost:8000/datasets/woelfel2020virological/export",
                       params={"format": "csv", "include_measurements": True})
if response.status_code == 200:
    with open("dataset.csv", "wb") as f:
        f.write(response.content)

# Get API metrics summary (admin endpoint)
response = requests.get("http://localhost:8000/admin/metrics")
metrics = response.json()

# Handle rate limiting
try:
    response = requests.get("http://localhost:8000/datasets/some_dataset")
    if response.status_code == 429:
        retry_after = response.json().get("retry_after", 60)
        print(f"Rate limited. Retry after {retry_after} seconds")
except requests.RequestException as e:
    print(f"Request failed: {e}")
```

### JavaScript/TypeScript

```javascript
// List datasets with error handling
try {
  const response = await fetch('http://localhost:8000/datasets?limit=10');
  if (!response.ok) {
    if (response.status === 422) {
      const errors = await response.json();
      console.log('Validation errors:', errors.detail);
    } else if (response.status === 429) {
      const rateLimitInfo = await response.json();
      console.log('Rate limited:', rateLimitInfo);
    }
    throw new Error(`HTTP ${response.status}`);
  }
  const datasets = await response.json();
} catch (error) {
  console.error('API request failed:', error);
}

// Search with validation
const searchParams = new URLSearchParams({
  q: 'stool',
  limit: '5'
});
const searchResponse = await fetch(`http://localhost:8000/search?${searchParams}`);
const searchResults = await searchResponse.json();

// Export dataset and handle file download
const exportResponse = await fetch('http://localhost:8000/datasets/woelfel2020virological/export?format=excel');
if (exportResponse.ok) {
  const blob = await exportResponse.blob();
  const url = window.URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'dataset.xlsx';
  a.click();
}
```

### cURL

```bash
# List all datasets
curl "http://localhost:8000/datasets"

# Get specific dataset without measurements (faster)
curl "http://localhost:8000/datasets/woelfel2020virological?include_measurements=false"

# Search with parameters
curl "http://localhost:8000/search?q=COVID&limit=5"

# Get available filter options
curl "http://localhost:8000/filters"

# Health check with system info
curl "http://localhost:8000/health"

# Export dataset as CSV
curl "http://localhost:8000/datasets/woelfel2020virological/export?format=csv" -o dataset.csv

# Export dataset as Excel without measurements (faster)
curl "http://localhost:8000/datasets/woelfel2020virological/export?format=excel&include_measurements=false" -o dataset.xlsx

# Get Prometheus metrics for monitoring
curl "http://localhost:8000/metrics"

# Get API metrics summary
curl "http://localhost:8000/admin/metrics"

# Test rate limiting (will show 429 after limit)
for i in {1..15}; do curl -s "http://localhost:8000/datasets/woelfel2020virological" -w " [%{http_code}]\n"; done
```

## Configuration

The API supports extensive configuration via environment variables:

```bash
# Environment
ENVIRONMENT=development  # development, production

# Server Configuration  
SHEDDING_HUB_HOST=0.0.0.0
SHEDDING_HUB_PORT=8000
SHEDDING_HUB_DEBUG=true

# Rate Limiting
SHEDDING_HUB_RATE_LIMIT_ENABLED=true
SHEDDING_HUB_RATE_LIMIT_DEFAULT="100/minute"
SHEDDING_HUB_RATE_LIMIT_SEARCH="30/minute"

# Caching
SHEDDING_HUB_CACHE_TTL_SECONDS=3600

# Data Source Configuration
SHEDDING_HUB_GITHUB_ENABLED=false
SHEDDING_HUB_GITHUB_REPO_OWNER=shedding-hub
SHEDDING_HUB_GITHUB_REPO_NAME=shedding-hub
SHEDDING_HUB_GITHUB_BRANCH=main
SHEDDING_HUB_GITHUB_DATA_PATH=data
SHEDDING_HUB_GITHUB_TOKEN=your_github_token  # Optional, for higher rate limits
SHEDDING_HUB_GITHUB_CACHE_TTL=600

# CORS (for production)
SHEDDING_HUB_CORS_ORIGINS="https://yourdomain.com,https://api.yourdomain.com"
```

## Project Structure

```
api/
├── __init__.py          # Makes it a Python package
├── main.py             # Main FastAPI application with all endpoints
├── models.py           # Pydantic models for validation and documentation
├── middleware.py       # Rate limiting, compression, and logging middleware
├── config.py           # Environment configuration management
├── async_utils.py      # Async utilities for data loading and export
├── metrics.py          # Prometheus metrics and monitoring
├── data_source.py      # Data source management (local/GitHub)
├── run.py              # Convenience script to start the API
├── requirements.txt    # API dependencies
├── .env.example        # Example environment configuration
├── test_api.py         # Comprehensive pytest test suite
├── test_rate_limiting.py # Rate limiting specific tests
└── README.md          # This file
```

## Response Format

All endpoints return JSON responses with consistent formatting and full type safety:

### Success Response
```json
{
  "total_found": 34,
  "returned": 10,
  "filters_applied": {
    "biomarker": "SARS-CoV-2",
    "specimen": null,
    "limit": 10
  },
  "datasets": [...] 
}
```

### Validation Error Response (422)
```json
{
  "detail": [
    {
      "type": "enum",
      "loc": ["query", "biomarker"],
      "msg": "Input should be 'SARS-CoV-2', 'mtDNA', 'PMMoV', 'crAssphage', 'influenza', 'sapovirus' or 'SARS'",
      "input": "invalid_biomarker"
    }
  ]
}
```

### Rate Limit Error Response (429)
```json
{
  "error": "Rate limit exceeded",
  "message": "Too many requests. Limit: 10 per 1 minute",
  "retry_after": 60
}
```

## Performance & Production Features

- **Smart Caching**: Dataset metadata cached in memory with configurable TTL
- **Rate Limiting**: IP-based limits prevent API abuse with per-endpoint customization
- **Request Validation**: All inputs validated with detailed error messages
- **Efficient Queries**: Use `include_measurements=false` for faster responses
- **Health Monitoring**: Detailed health checks with system status
- **Structured Logging**: Comprehensive JSON logging with request correlation
- **Environment-based Config**: Different settings for development/production
- **Type Safety**: Full Pydantic validation prevents runtime errors
- **Data Export**: Multiple export formats (CSV, JSON, Excel) with streaming
- **Prometheus Metrics**: Built-in metrics collection for monitoring and alerting
- **Compression**: Automatic response compression for large payloads
- **GitHub Integration**: Optional remote data loading from GitHub repositories
- **Admin Dashboards**: Detailed metrics endpoints for monitoring dashboards

## Security Features

- **Input Validation**: All parameters validated against strict schemas with sanitization
- **Rate Limiting**: Configurable per-endpoint limits with IP tracking
- **CORS Configuration**: Configurable origins for production security
- **Error Handling**: Detailed validation errors without information leakage
- **Request Size Limits**: Configurable maximum request sizes
- **Input Sanitization**: Automatic removal of control characters and length limits
- **Dataset ID Validation**: Strict alphanumeric validation for dataset identifiers
- **Search Query Sanitization**: Safe text processing for search queries
- **Structured Error Responses**: Consistent error format with timestamps and correlation

## Development

### Running in Development Mode

```bash
# Install dependencies
pip install -r requirements.txt

# Copy environment template
cp .env.example .env

# Run with auto-reload (recommended)
python run.py

# Or with uvicorn directly  
uvicorn main:app --reload --port 8000

# Run tests with coverage
pytest test_api.py --cov=main --cov-report=html -v
```

### Adding New Endpoints

1. Add endpoint function to `main.py` with appropriate rate limiting decorator
2. Create Pydantic models in `models.py` for request/response validation
3. Add input validation and sanitization if needed
4. Update metrics collection if applicable
5. Add comprehensive tests to `test_api.py`
6. Update this README with documentation
7. Test validation, rate limiting, and security

### Production Deployment

```bash
# Set production environment
export ENVIRONMENT=production
export SHEDDING_HUB_CORS_ORIGINS="https://yourdomain.com"
export SHEDDING_HUB_DEBUG=false
export SHEDDING_HUB_RATE_LIMIT_ENABLED=true

# Optional: Enable GitHub data source for production
export SHEDDING_HUB_GITHUB_ENABLED=true
export SHEDDING_HUB_GITHUB_TOKEN="your_production_token"

# Run with production server
gunicorn main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000

# Or with Prometheus metrics monitoring
gunicorn main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000 --access-logfile - --error-logfile -
```

## Available Datasets

The API automatically discovers all `.yaml` files in the `../data/` directory structure. The shedding-hub contains datasets from various studies related to biomarker shedding, primarily focused on SARS-CoV-2 but also including other pathogens and biomarkers.

Current dataset count: **34 datasets** with **1,800+ participants** and **comprehensive validation**.

## License

This API serves data from the Shedding Hub project. Please refer to the main project's license and citation requirements when using the data.