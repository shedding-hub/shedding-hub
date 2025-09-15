"""
Metrics collection and monitoring for the Shedding Hub API.

Provides Prometheus-style metrics and system monitoring.
"""

import time
import psutil
from datetime import datetime
from typing import Dict, Any
from prometheus_client import Counter, Histogram, Gauge, Info, generate_latest, CONTENT_TYPE_LATEST
from fastapi import Request
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import Response


# Prometheus metrics
api_requests_total = Counter(
    'api_requests_total', 
    'Total API requests', 
    ['method', 'endpoint', 'status_code']
)

api_request_duration_seconds = Histogram(
    'api_request_duration_seconds', 
    'API request duration in seconds',
    ['method', 'endpoint']
)

api_cache_hits_total = Counter(
    'api_cache_hits_total', 
    'Total cache hits', 
    ['cache_type']
)

api_cache_misses_total = Counter(
    'api_cache_misses_total', 
    'Total cache misses', 
    ['cache_type']
)

api_errors_total = Counter(
    'api_errors_total', 
    'Total API errors', 
    ['error_type', 'endpoint']
)

api_active_connections = Gauge(
    'api_active_connections', 
    'Number of active connections'
)

api_memory_usage_bytes = Gauge(
    'api_memory_usage_bytes', 
    'Memory usage in bytes'
)

api_cpu_usage_percent = Gauge(
    'api_cpu_usage_percent', 
    'CPU usage percentage'
)

api_datasets_loaded = Gauge(
    'api_datasets_loaded', 
    'Number of datasets currently loaded'
)

api_rate_limit_hits = Counter(
    'api_rate_limit_hits_total', 
    'Total rate limit violations', 
    ['endpoint', 'client_ip']
)

# API info metric
api_info = Info(
    'api_info', 
    'API version and build information'
)


class MetricsMiddleware(BaseHTTPMiddleware):
    """Middleware to collect request metrics."""
    
    async def dispatch(self, request: Request, call_next):
        start_time = time.time()
        method = request.method
        endpoint = request.url.path
        
        # Track active connections
        api_active_connections.inc()
        
        try:
            response = await call_next(request)
            status_code = str(response.status_code)
            
            # Record request metrics
            api_requests_total.labels(
                method=method, 
                endpoint=endpoint, 
                status_code=status_code
            ).inc()
            
            # Record request duration
            duration = time.time() - start_time
            api_request_duration_seconds.labels(
                method=method,
                endpoint=endpoint
            ).observe(duration)
            
            # Track rate limit hits
            if response.status_code == 429:
                client_ip = request.client.host if request.client else "unknown"
                api_rate_limit_hits.labels(
                    endpoint=endpoint,
                    client_ip=client_ip
                ).inc()
            
            return response
            
        except Exception as e:
            # Record error metrics
            error_type = type(e).__name__
            api_errors_total.labels(
                error_type=error_type,
                endpoint=endpoint
            ).inc()
            raise
        finally:
            # Decrease active connections
            api_active_connections.dec()


def update_system_metrics():
    """Update system-level metrics."""
    try:
        # Memory usage
        process = psutil.Process()
        memory_info = process.memory_info()
        api_memory_usage_bytes.set(memory_info.rss)
        
        # CPU usage (averaged over 1 second)
        cpu_percent = process.cpu_percent()
        api_cpu_usage_percent.set(cpu_percent)
        
    except Exception:
        # Silently fail if psutil is not available or fails
        pass


def record_cache_hit(cache_type: str = "dataset"):
    """Record a cache hit."""
    api_cache_hits_total.labels(cache_type=cache_type).inc()


def record_cache_miss(cache_type: str = "dataset"):
    """Record a cache miss."""
    api_cache_misses_total.labels(cache_type=cache_type).inc()


def update_datasets_count(count: int):
    """Update the number of loaded datasets."""
    api_datasets_loaded.set(count)


def set_api_info(version: str, build_date: str = None):
    """Set API information metrics."""
    info_data = {
        'version': version,
        'build_date': build_date or datetime.now().isoformat()
    }
    api_info.info(info_data)


def get_prometheus_metrics() -> str:
    """Get Prometheus metrics in text format."""
    # Update system metrics before generating output
    update_system_metrics()
    return generate_latest().decode('utf-8')


def get_metrics_summary() -> Dict[str, Any]:
    """Get a summary of key metrics for the health endpoint."""
    try:
        process = psutil.Process()
        memory_mb = process.memory_info().rss / (1024 * 1024)
        cpu_percent = process.cpu_percent()
        
        return {
            "memory_usage_mb": round(memory_mb, 2),
            "cpu_usage_percent": round(cpu_percent, 2),
            "uptime_seconds": int(time.time() - psutil.boot_time()),
            "datasets_loaded": int(api_datasets_loaded._value._value),
            "total_requests": int(sum(api_requests_total._value.values())),
            "cache_hit_rate": _calculate_cache_hit_rate(),
            "error_rate": _calculate_error_rate()
        }
    except Exception:
        return {
            "memory_usage_mb": 0,
            "cpu_usage_percent": 0,
            "uptime_seconds": 0,
            "datasets_loaded": 0,
            "total_requests": 0,
            "cache_hit_rate": 0.0,
            "error_rate": 0.0
        }


def _calculate_cache_hit_rate() -> float:
    """Calculate cache hit rate."""
    try:
        hits = sum(api_cache_hits_total._value.values())
        misses = sum(api_cache_misses_total._value.values())
        total = hits + misses
        return round(hits / total, 3) if total > 0 else 0.0
    except Exception:
        return 0.0


def _calculate_error_rate() -> float:
    """Calculate error rate based on HTTP status codes."""
    try:
        total_requests = sum(api_requests_total._value.values())
        error_requests = sum(
            count for labels, count in api_requests_total._value.items()
            if labels[2].startswith(('4', '5'))  # 4xx and 5xx status codes
        )
        return round(error_requests / total_requests, 3) if total_requests > 0 else 0.0
    except Exception:
        return 0.0


def get_detailed_metrics() -> Dict[str, Any]:
    """Get detailed metrics for admin endpoints."""
    return {
        "requests": {
            "total": int(sum(api_requests_total._value.values())),
            "by_endpoint": _get_requests_by_endpoint(),
            "by_status": _get_requests_by_status()
        },
        "performance": {
            "avg_response_time_ms": _get_avg_response_time(),
            "cache_hit_rate": _calculate_cache_hit_rate(),
            "error_rate": _calculate_error_rate()
        },
        "system": get_metrics_summary(),
        "rate_limiting": {
            "violations_total": int(sum(api_rate_limit_hits.labels('', '')._value.values())),
            "by_endpoint": _get_rate_limits_by_endpoint()
        }
    }


def _get_requests_by_endpoint() -> Dict[str, int]:
    """Get request counts by endpoint."""
    endpoint_counts = {}
    for labels, count in api_requests_total._value.items():
        endpoint = labels[1]  # endpoint is the second label
        endpoint_counts[endpoint] = endpoint_counts.get(endpoint, 0) + count
    return endpoint_counts


def _get_requests_by_status() -> Dict[str, int]:
    """Get request counts by status code."""
    status_counts = {}
    for labels, count in api_requests_total._value.items():
        status = labels[2]  # status_code is the third label
        status_counts[status] = status_counts.get(status, 0) + count
    return status_counts


def _get_avg_response_time() -> float:
    """Get average response time in milliseconds."""
    try:
        total_time = 0
        total_count = 0
        
        for labels, histogram in api_request_duration_seconds._value.items():
            total_time += histogram._sum._value
            total_count += histogram._count._value
        
        if total_count > 0:
            avg_seconds = total_time / total_count
            return round(avg_seconds * 1000, 2)  # Convert to milliseconds
        return 0.0
    except Exception:
        return 0.0


def _get_rate_limits_by_endpoint() -> Dict[str, int]:
    """Get rate limit violations by endpoint."""
    endpoint_violations = {}
    for labels, count in api_rate_limit_hits._value.items():
        endpoint = labels[0]  # endpoint is the first label
        endpoint_violations[endpoint] = endpoint_violations.get(endpoint, 0) + count
    return endpoint_violations