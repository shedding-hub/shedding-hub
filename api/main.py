from fastapi import FastAPI, HTTPException, Query, Depends, Request
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.exceptions import RequestValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException
from typing import List, Dict, Any, Optional
import os
import yaml
import pathlib
from datetime import datetime
import sys
import logging
from slowapi.errors import RateLimitExceeded
try:
    from .models import (
        DatasetListResponse, FullDataset, DatasetMetadata, SearchResponse,
        ProjectStats, HealthStatus, AvailableFilters, CacheRefreshResponse,
        ApiInfo, DatasetSummary, SearchResult, DatasetStats, ErrorDetail,
        BiomarkerType, SpecimenType, ExportFormat, ExportResponse,
        MetricsSummary, DetailedMetrics
    )
    from .middleware import (
        limiter, create_rate_limiting_middleware, custom_rate_limit_handler,
        rate_limit_default, rate_limit_search, rate_limit_expensive, rate_limit_cache_operations,
        CompressionMiddleware, LoggingMiddleware, setup_structured_logging
    )
    from .config import settings
    from .async_utils import (
        load_all_datasets_async, export_dataset_to_csv_async,
        export_dataset_to_json_async, export_dataset_to_excel_async
    )
    from .metrics import (
        MetricsMiddleware, get_prometheus_metrics, get_metrics_summary, 
        get_detailed_metrics, set_api_info, update_datasets_count, 
        record_cache_hit, record_cache_miss, CONTENT_TYPE_LATEST
    )
except ImportError:
    # Handle direct execution
    from models import (
        DatasetListResponse, FullDataset, DatasetMetadata, SearchResponse,
        ProjectStats, HealthStatus, AvailableFilters, CacheRefreshResponse,
        ApiInfo, DatasetSummary, SearchResult, DatasetStats, ErrorDetail,
        BiomarkerType, SpecimenType, ExportFormat, ExportResponse,
        MetricsSummary, DetailedMetrics
    )
    from middleware import (
        limiter, create_rate_limiting_middleware, custom_rate_limit_handler,
        rate_limit_default, rate_limit_search, rate_limit_expensive, rate_limit_cache_operations,
        CompressionMiddleware, LoggingMiddleware, setup_structured_logging
    )
    from config import settings
    from async_utils import (
        load_all_datasets_async, export_dataset_to_csv_async,
        export_dataset_to_json_async, export_dataset_to_excel_async
    )
    from metrics import (
        MetricsMiddleware, get_prometheus_metrics, get_metrics_summary, 
        get_detailed_metrics, set_api_info, update_datasets_count, 
        record_cache_hit, record_cache_miss, CONTENT_TYPE_LATEST
    )
    from data_source import data_source_manager

# Configure structured JSON logging
setup_structured_logging()
logger = logging.getLogger(__name__)

# Add the parent directory to the path so we can import shedding_hub
sys.path.append(str(pathlib.Path(__file__).parent.parent))
import shedding_hub as sh
import re

app = FastAPI(
    title=settings.api_title,
    description=settings.api_description,
    version=settings.api_version,
    docs_url="/docs",
    redoc_url="/redoc",
    debug=settings.debug,
)

# Add rate limiting state
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, custom_rate_limit_handler)

# Global exception handlers
@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Handle Pydantic validation errors."""
    logger.error(f"Validation error on {request.url.path}: {exc.errors()}")
    return JSONResponse(
        status_code=422,
        content={
            "error": "Validation Error",
            "message": "The provided data is invalid",
            "details": exc.errors(),
            "path": request.url.path,
            "timestamp": datetime.now().isoformat()
        }
    )

@app.exception_handler(StarletteHTTPException)
async def http_exception_handler(request: Request, exc: StarletteHTTPException):
    """Handle HTTP exceptions with structured error format."""
    logger.error(f"HTTP {exc.status_code} error on {request.url.path}: {exc.detail}")
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "error": f"HTTP {exc.status_code}",
            "message": exc.detail,
            "path": request.url.path,
            "timestamp": datetime.now().isoformat()
        }
    )

@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    """Handle unexpected exceptions."""
    logger.error(f"Unexpected error on {request.url.path}: {str(exc)}", exc_info=True)
    return JSONResponse(
        status_code=500,
        content={
            "error": "Internal Server Error",
            "message": "An unexpected error occurred. Please try again later.",
            "path": request.url.path,
            "timestamp": datetime.now().isoformat()
        }
    )

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    **settings.get_cors_config()
)

# Add logging middleware (first, to capture all requests)
app.add_middleware(LoggingMiddleware)

# Add metrics middleware (early, to track all requests)
app.add_middleware(MetricsMiddleware)

# Add compression middleware (before rate limiting)
app.add_middleware(CompressionMiddleware, minimum_size=1000)

# Add rate limiting middleware
if settings.rate_limit_enabled:
    app.add_middleware(create_rate_limiting_middleware())

# Global configuration is now handled by settings
# Access via settings.cache_ttl_seconds, settings.max_datasets_limit, etc.

# Global variable to store dataset metadata cache
_dataset_cache = {}
_cache_timestamp = None

def get_data_dir() -> pathlib.Path:
    """Get the data directory path."""
    return settings.get_data_directory()

def is_cache_expired() -> bool:
    """Check if the cache has expired."""
    if _cache_timestamp is None:
        return True
    
    return (datetime.now() - _cache_timestamp).total_seconds() > settings.cache_ttl_seconds


# Input sanitization and validation utilities
def sanitize_string(value: str, max_length: int = 1000) -> str:
    """Sanitize string input to prevent injection attacks."""
    if not value:
        return ""
    
    # Remove null bytes and control characters
    sanitized = re.sub(r'[\x00-\x1f\x7f]', '', value)
    
    # Limit length
    if len(sanitized) > max_length:
        sanitized = sanitized[:max_length]
    
    # Strip whitespace
    return sanitized.strip()


def validate_dataset_id(dataset_id: str) -> str:
    """Validate and sanitize dataset ID."""
    if not dataset_id:
        raise HTTPException(status_code=400, detail="Dataset ID cannot be empty")
    
    # Dataset IDs should be alphanumeric with underscores and hyphens
    if not re.match(r'^[a-zA-Z0-9_-]+$', dataset_id):
        raise HTTPException(
            status_code=400, 
            detail="Dataset ID can only contain letters, numbers, underscores, and hyphens"
        )
    
    if len(dataset_id) > 100:
        raise HTTPException(status_code=400, detail="Dataset ID too long (max 100 characters)")
    
    return dataset_id


def validate_search_query(query: str) -> str:
    """Validate and sanitize search query."""
    if not query:
        raise HTTPException(status_code=400, detail="Search query cannot be empty")
    
    sanitized = sanitize_string(query, max_length=500)
    
    if len(sanitized.strip()) < 1:
        raise HTTPException(status_code=400, detail="Search query too short")
    
    return sanitized


def validate_pagination_params(limit: int, offset: int = 0) -> tuple:
    """Validate pagination parameters."""
    if limit < 1:
        raise HTTPException(status_code=400, detail="Limit must be at least 1")
    
    if limit > settings.max_datasets_limit:
        raise HTTPException(
            status_code=400, 
            detail=f"Limit cannot exceed {settings.max_datasets_limit}"
        )
    
    if offset < 0:
        raise HTTPException(status_code=400, detail="Offset cannot be negative")
    
    return limit, offset

# Legacy function for local file loading (kept for backward compatibility)
async def load_all_datasets_local(force_reload: bool = False) -> Dict[str, Dict[str, Any]]:
    """Load all datasets from local files and cache their metadata asynchronously."""
    global _dataset_cache, _cache_timestamp
    
    if not force_reload and _dataset_cache and not is_cache_expired():
        logger.info("Using cached datasets")
        record_cache_hit("dataset")
        return _dataset_cache
    
    logger.info("Loading datasets from disk asynchronously")
    record_cache_miss("dataset")
    data_dir = get_data_dir()
    
    if not data_dir.exists():
        logger.error(f"Data directory not found: {data_dir}")
        raise HTTPException(status_code=500, detail="Data directory not found")
    
    try:
        # Load all datasets asynchronously with controlled concurrency
        datasets = await load_all_datasets_async(data_dir, max_concurrent=10)
        
        _dataset_cache = datasets
        _cache_timestamp = datetime.now()
        update_datasets_count(len(datasets))
        logger.info(f"Loaded {len(datasets)} datasets into cache asynchronously")
        return datasets
        
    except Exception as e:
        logger.error(f"Error loading datasets asynchronously: {e}")
        raise HTTPException(status_code=500, detail=f"Error loading datasets: {str(e)}")

# New unified function using data source manager
async def load_all_datasets(force_reload: bool = False) -> Dict[str, Dict[str, Any]]:
    """Load all datasets using the configured data source (local or GitHub)."""
    try:
        if settings.github_enabled:
            # For GitHub source, we need to build a dictionary from the dataset list
            dataset_ids = await data_source_manager.get_dataset_list(force_reload)
            datasets = {}
            
            # Load each dataset (this will be cached by the GitHub data source)
            for dataset_id in dataset_ids:
                try:
                    dataset = await data_source_manager.get_dataset(dataset_id, force_reload)
                    datasets[dataset_id] = dataset
                except Exception as e:
                    logger.error(f"Error loading dataset {dataset_id}: {e}")
                    continue
            
            update_datasets_count(len(datasets))
            return datasets
        else:
            # Use local file loading
            return await load_all_datasets_local(force_reload)
            
    except Exception as e:
        logger.error(f"Error in unified dataset loading: {e}")
        raise HTTPException(status_code=500, detail=f"Error loading datasets: {str(e)}")

@app.on_event("startup")
async def startup_event():
    """Initialize metrics and load initial data on startup."""
    set_api_info(settings.api_version)
    logger.info("API startup complete with metrics initialized")

@app.get("/", response_model=ApiInfo)
@rate_limit_default()
async def root(request: Request):
    """Root endpoint with basic API information."""
    return {
        "message": "Welcome to the Shedding Hub API",
        "version": "0.1.0",
        "description": "API for accessing biomarker shedding datasets",
        "endpoints": {
            "datasets": "/datasets - List all available datasets",
            "dataset": "/datasets/{dataset_id} - Get specific dataset",
            "metadata": "/datasets/{dataset_id}/metadata - Get dataset metadata only",
            "search": "/search - Search datasets by biomarker, specimen, etc.",
            "stats": "/stats - Get overall project statistics",
            "health": "/health - Health check endpoint",
            "metrics": "/metrics - Prometheus metrics for monitoring",
            "admin_metrics": "/admin/metrics - API metrics summary",
            "docs": "/docs - Interactive API documentation"
        },
        "cache_info": {
            "cache_ttl_seconds": settings.cache_ttl_seconds,
            "cache_expired": is_cache_expired(),
            "last_cache_update": _cache_timestamp.isoformat() if _cache_timestamp else None
        }
    }

@app.get("/datasets", response_model=DatasetListResponse)
@rate_limit_default()
async def list_datasets(
    request: Request,
    biomarker: Optional[BiomarkerType] = Query(None, description="Filter by biomarker type"),
    specimen: Optional[SpecimenType] = Query(None, description="Filter by specimen type"),
    limit: int = Query(settings.default_datasets_limit, description="Maximum number of datasets to return", le=settings.max_datasets_limit, ge=1)
):
    """List all available datasets with optional filtering."""
    try:
        all_datasets = await load_all_datasets()
        
        # Apply filters
        filtered_datasets = {}
        for dataset_id, dataset in all_datasets.items():
            include = True
            
            if biomarker:
                # Check if any analyte has the specified biomarker
                analytes = dataset.get('analytes', {})
                biomarker_found = any(
                    analyte.get('biomarker', '').lower() == biomarker.lower()
                    for analyte in analytes.values()
                )
                if not biomarker_found:
                    include = False
            
            if specimen and include:
                # Check if any analyte has the specified specimen
                analytes = dataset.get('analytes', {})
                specimen_found = False
                for analyte in analytes.values():
                    specimen_data = analyte.get('specimen', '')
                    if isinstance(specimen_data, list):
                        specimen_found = any(s.lower() == specimen.lower() for s in specimen_data)
                    else:
                        specimen_found = specimen_data.lower() == specimen.lower()
                    if specimen_found:
                        break
                if not specimen_found:
                    include = False
            
            if include:
                filtered_datasets[dataset_id] = dataset
        
        # Limit results
        limited_datasets = dict(list(filtered_datasets.items())[:limit])
        
        # Return summary information
        result = []
        for dataset_id, dataset in limited_datasets.items():
            result.append({
                'dataset_id': dataset_id,
                'title': dataset.get('title', ''),
                'doi': dataset.get('doi'),
                'url': dataset.get('url'),
                'description': dataset.get('description', '')[:200] + '...' if len(dataset.get('description', '')) > 200 else dataset.get('description', ''),
                'stats': dataset.get('stats', {}),
                'analytes': list(dataset.get('analytes', {}).keys()),
                'last_modified': dataset.get('last_modified')
            })
        
        return {
            'total_found': len(filtered_datasets),
            'returned': len(result),
            'filters_applied': {
                'biomarker': biomarker,
                'specimen': specimen,
                'limit': limit
            },
            'datasets': result
        }
        
    except Exception as e:
        logger.error(f"Error loading datasets: {e}")
        raise HTTPException(status_code=500, detail=f"Error loading datasets: {str(e)}")

@app.get("/datasets/{dataset_id}", response_model=FullDataset)
@rate_limit_expensive()
async def get_dataset(request: Request, dataset_id: str, include_measurements: bool = Query(True, description="Include participant measurements")):
    """Get a specific dataset by ID."""
    try:
        # Validate and sanitize dataset_id
        dataset_id = validate_dataset_id(dataset_id)
        
        all_datasets = await load_all_datasets()
        
        if dataset_id not in all_datasets:
            logger.warning(f"Dataset '{dataset_id}' not found")
            raise HTTPException(status_code=404, detail=f"Dataset '{dataset_id}' not found")
        
        dataset = all_datasets[dataset_id].copy()
        
        # Optionally exclude measurements to reduce payload size
        if not include_measurements:
            for participant in dataset.get('participants', []):
                measurement_count = len(participant.get('measurements', []))
                participant['measurements'] = f"[{measurement_count} measurements - use include_measurements=true to see them]"
        
        logger.info(f"Retrieved dataset: {dataset_id}")
        return dataset
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error loading dataset {dataset_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Error loading dataset: {str(e)}")

@app.get("/datasets/{dataset_id}/metadata", response_model=DatasetMetadata)
@rate_limit_default()
async def get_dataset_metadata(request: Request, dataset_id: str):
    """Get metadata for a specific dataset (without measurements)."""
    try:
        # Validate and sanitize dataset_id
        dataset_id = validate_dataset_id(dataset_id)
        
        all_datasets = await load_all_datasets()
        
        if dataset_id not in all_datasets:
            logger.warning(f"Dataset '{dataset_id}' not found")
            raise HTTPException(status_code=404, detail=f"Dataset '{dataset_id}' not found")
        
        dataset = all_datasets[dataset_id].copy()
        
        # Remove measurements to keep only metadata
        metadata = {
            'dataset_id': dataset.get('dataset_id'),
            'title': dataset.get('title'),
            'doi': dataset.get('doi'),
            'url': dataset.get('url'),
            'description': dataset.get('description'),
            'analytes': dataset.get('analytes', {}),
            'stats': dataset.get('stats'),
            'last_modified': dataset.get('last_modified'),
            'participants_attributes': []
        }
        
        # Include participant attributes but not measurements
        for participant in dataset.get('participants', []):
            if 'attributes' in participant:
                metadata['participants_attributes'].append(participant['attributes'])
        
        logger.info(f"Retrieved metadata for dataset: {dataset_id}")
        return metadata
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error loading dataset metadata {dataset_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Error loading dataset metadata: {str(e)}")

@app.get("/search", response_model=SearchResponse)
@rate_limit_search()
async def search_datasets(
    request: Request,
    q: str = Query(..., description="Search query", min_length=1),
    fields: List[str] = Query(['title', 'description'], description="Fields to search in"),
    limit: int = Query(settings.default_search_limit, description="Maximum number of results", le=settings.max_datasets_limit, ge=1)
):
    """Search datasets by text query."""
    try:
        # Validate and sanitize search query
        q = validate_search_query(q)
        
        # Validate pagination parameters
        limit, _ = validate_pagination_params(limit)
        
        # Validate search fields
        valid_fields = {'title', 'description', 'doi', 'url'}
        fields = [field for field in fields if field in valid_fields]
        if not fields:
            fields = ['title', 'description']  # Default fallback
        
        all_datasets = await load_all_datasets()
        results = []
        
        query_lower = q.lower()
        
        for dataset_id, dataset in all_datasets.items():
            match_score = 0
            
            # Search in specified fields
            for field in fields:
                if field in dataset:
                    field_value = str(dataset[field]).lower()
                    if query_lower in field_value:
                        match_score += field_value.count(query_lower)
            
            # Also search in analyte descriptions
            if 'analytes' in dataset:
                for analyte in dataset['analytes'].values():
                    if 'description' in analyte:
                        if query_lower in analyte['description'].lower():
                            match_score += 1
            
            if match_score > 0:
                results.append({
                    'dataset_id': dataset_id,
                    'title': dataset.get('title', ''),
                    'description': dataset.get('description', '')[:200] + '...' if len(dataset.get('description', '')) > 200 else dataset.get('description', ''),
                    'match_score': match_score,
                    'stats': dataset.get('stats', {})
                })
        
        # Sort by match score descending
        results.sort(key=lambda x: x['match_score'], reverse=True)
        
        logger.info(f"Search query '{q}' returned {len(results)} results")
        return {
            'query': q,
            'fields_searched': fields,
            'total_found': len(results),
            'returned': min(len(results), limit),
            'results': results[:limit]
        }
        
    except Exception as e:
        logger.error(f"Error searching datasets: {e}")
        raise HTTPException(status_code=500, detail=f"Error searching datasets: {str(e)}")

@app.get("/stats", response_model=ProjectStats)
@rate_limit_default()
async def get_stats(request: Request):
    """Get overall statistics about all datasets."""
    try:
        all_datasets = await load_all_datasets()
        
        total_datasets = len(all_datasets)
        total_participants = sum(d.get('stats', {}).get('participants_count', 0) for d in all_datasets.values())
        total_measurements = sum(d.get('stats', {}).get('total_measurements', 0) for d in all_datasets.values())
        
        # Count unique biomarkers and specimens
        biomarkers = set()
        specimens = set()
        
        for dataset in all_datasets.values():
            for analyte in dataset.get('analytes', {}).values():
                if 'biomarker' in analyte:
                    biomarkers.add(analyte['biomarker'])
                if 'specimen' in analyte:
                    specimen_data = analyte['specimen']
                    if isinstance(specimen_data, list):
                        specimens.update(specimen_data)
                    else:
                        specimens.add(specimen_data)
        
        logger.info("Retrieved project statistics")
        return {
            'total_datasets': total_datasets,
            'total_participants': total_participants,
            'total_measurements': total_measurements,
            'unique_biomarkers': sorted(list(biomarkers)),
            'unique_specimens': sorted(list(specimens)),
            'biomarker_count': len(biomarkers),
            'specimen_count': len(specimens),
            'cache_info': {
                'cache_expired': is_cache_expired(),
                'last_cache_update': _cache_timestamp.isoformat() if _cache_timestamp else None
            }
        }
        
    except Exception as e:
        logger.error(f"Error calculating stats: {e}")
        raise HTTPException(status_code=500, detail=f"Error calculating stats: {str(e)}")

@app.post("/refresh-cache", response_model=CacheRefreshResponse)
@rate_limit_cache_operations()
async def refresh_cache(request: Request):
    """Refresh the dataset cache."""
    try:
        global _dataset_cache, _cache_timestamp
        _dataset_cache = {}
        _cache_timestamp = None
        datasets = await load_all_datasets(force_reload=True)
        logger.info("Cache refreshed successfully")
        return {
            "message": "Cache refreshed successfully",
            "datasets_loaded": len(datasets),
            "cache_timestamp": _cache_timestamp.isoformat()
        }
    except Exception as e:
        logger.error(f"Error refreshing cache: {e}")
        raise HTTPException(status_code=500, detail=f"Error refreshing cache: {str(e)}")

# Health check endpoint
@app.get("/health", response_model=HealthStatus)
@rate_limit_default()
async def health_check(request: Request):
    """Health check endpoint."""
    try:
        # Try to load datasets to verify everything is working
        datasets = await load_all_datasets()
        
        return {
            "status": "healthy",
            "datasets_available": len(datasets),
            "cache_status": "fresh" if not is_cache_expired() else "expired",
            "data_directory_exists": get_data_dir().exists(),
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        raise HTTPException(status_code=503, detail=f"Service unhealthy: {str(e)}")

# Additional endpoint to get available biomarkers and specimens
@app.get("/filters", response_model=AvailableFilters)
@rate_limit_default()
async def get_available_filters(request: Request):
    """Get available filter options for biomarkers and specimens."""
    try:
        all_datasets = await load_all_datasets()
        
        biomarkers = set()
        specimens = set()
        
        for dataset in all_datasets.values():
            for analyte in dataset.get('analytes', {}).values():
                if 'biomarker' in analyte:
                    biomarkers.add(analyte['biomarker'])
                if 'specimen' in analyte:
                    specimen_data = analyte['specimen']
                    if isinstance(specimen_data, list):
                        specimens.update(specimen_data)
                    else:
                        specimens.add(specimen_data)
        
        return {
            'available_biomarkers': sorted(list(biomarkers)),
            'available_specimens': sorted(list(specimens))
        }
        
    except Exception as e:
        logger.error(f"Error getting available filters: {e}")
        raise HTTPException(status_code=500, detail=f"Error getting available filters: {str(e)}")


@app.get("/datasets/{dataset_id}/export", response_model=ExportResponse)
@rate_limit_expensive()
async def export_dataset(
    request: Request,
    dataset_id: str,
    format: ExportFormat = Query(..., description="Export format"),
    include_measurements: bool = Query(True, description="Include participant measurements")
):
    """Export a dataset in the specified format."""
    try:
        # Validate and sanitize dataset_id
        dataset_id = validate_dataset_id(dataset_id)
        
        all_datasets = await load_all_datasets()
        
        if dataset_id not in all_datasets:
            logger.warning(f"Dataset '{dataset_id}' not found for export")
            raise HTTPException(status_code=404, detail=f"Dataset '{dataset_id}' not found")
        
        dataset = all_datasets[dataset_id].copy()
        
        # Optionally exclude measurements
        if not include_measurements:
            for participant in dataset.get('participants', []):
                participant['measurements'] = []
        
        # Generate filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{dataset_id}_{timestamp}.{format.value}"
        
        # Export based on format
        if format == ExportFormat.csv:
            content = await export_dataset_to_csv_async(dataset)
            media_type = "text/csv"
            content_bytes = content.encode('utf-8')
            
        elif format == ExportFormat.json:
            content = await export_dataset_to_json_async(dataset, pretty=True)
            media_type = "application/json"
            content_bytes = content.encode('utf-8')
            
        elif format == ExportFormat.excel:
            content_bytes = await export_dataset_to_excel_async(dataset)
            media_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        
        else:
            raise HTTPException(status_code=400, detail="Unsupported export format")
        
        logger.info(f"Dataset '{dataset_id}' exported in {format.value} format")
        
        # Return file as download
        from fastapi.responses import Response
        return Response(
            content=content_bytes,
            media_type=media_type,
            headers={
                "Content-Disposition": f"attachment; filename={filename}",
                "Content-Length": str(len(content_bytes))
            }
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error exporting dataset {dataset_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Error exporting dataset: {str(e)}")


@app.get("/metrics", include_in_schema=False)
async def prometheus_metrics():
    """Prometheus metrics endpoint for monitoring and alerting."""
    from fastapi.responses import PlainTextResponse
    metrics_data = get_prometheus_metrics()
    return PlainTextResponse(content=metrics_data, media_type=CONTENT_TYPE_LATEST)


@app.get("/admin/metrics", response_model=MetricsSummary)
@rate_limit_expensive()
async def get_api_metrics_summary(request: Request):
    """Get API metrics summary for monitoring dashboards."""
    try:
        return get_metrics_summary()
    except Exception as e:
        logger.error(f"Error retrieving metrics summary: {e}")
        raise HTTPException(status_code=500, detail="Error retrieving metrics")


@app.get("/admin/metrics/detailed", response_model=DetailedMetrics)
@rate_limit_expensive()
async def get_detailed_metrics_endpoint(request: Request):
    """Get detailed metrics for admin monitoring and debugging."""
    try:
        return get_detailed_metrics()
    except Exception as e:
        logger.error(f"Error retrieving detailed metrics: {e}")
        raise HTTPException(status_code=500, detail="Error retrieving detailed metrics")


if __name__ == "__main__":
    import uvicorn
    # For direct execution, disable reload to avoid import string issues
    uvicorn.run(
        app,
        host=settings.host,
        port=settings.port,
        reload=False,  # Disabled for direct execution
        log_level=settings.log_level.lower(),
        access_log=True
    )
