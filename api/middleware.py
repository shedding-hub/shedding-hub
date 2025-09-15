"""
Middleware components for the Shedding Hub API.
"""

from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded
from slowapi.middleware import SlowAPIMiddleware
from fastapi import Request
from fastapi.responses import JSONResponse
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import Response
from starlette.types import ASGIApp
import logging
import gzip
import json
import uuid
import time
from contextvars import ContextVar

logger = logging.getLogger(__name__)

# Context variable for request correlation ID
correlation_id: ContextVar[str] = ContextVar("correlation_id", default="")

# Create limiter instance
limiter = Limiter(key_func=get_remote_address)


def create_rate_limiting_middleware():
    """Create and configure rate limiting middleware."""
    return SlowAPIMiddleware


def custom_rate_limit_handler(request: Request, exc: RateLimitExceeded):
    """Custom rate limit exceeded handler."""
    client_ip = get_remote_address(request)
    logger.warning(f"Rate limit exceeded for IP: {client_ip}, Path: {request.url.path}")
    
    # Get retry_after safely
    retry_after = getattr(exc, 'retry_after', 60)  # Default to 60 seconds
    
    return JSONResponse(
        status_code=429,
        content={
            "error": "Rate limit exceeded",
            "message": f"Too many requests. Limit: {exc.detail}",
            "retry_after": retry_after
        }
    )


# Rate limiting decorators for different endpoint types
def rate_limit_default():
    """Default rate limit: 100 requests per minute"""
    return limiter.limit("100/minute")


def rate_limit_search():
    """Search endpoint rate limit: 30 requests per minute"""
    return limiter.limit("30/minute")


def rate_limit_expensive():
    """Expensive operations: 10 requests per minute"""
    return limiter.limit("10/minute")


def rate_limit_cache_operations():
    """Cache operations: 5 requests per minute"""
    return limiter.limit("5/minute")


class CompressionMiddleware(BaseHTTPMiddleware):
    """Custom compression middleware for JSON responses."""
    
    def __init__(self, app: ASGIApp, minimum_size: int = 1000):
        super().__init__(app)
        self.minimum_size = minimum_size
    
    async def dispatch(self, request: Request, call_next):
        response = await call_next(request)
        
        # Only compress JSON responses that are large enough
        if (
            response.status_code == 200 and
            response.headers.get("content-type", "").startswith("application/json") and
            "content-encoding" not in response.headers and
            "gzip" in request.headers.get("accept-encoding", "")
        ):
            # Get response body
            response_body = b""
            async for chunk in response.body_iterator:
                response_body += chunk
            
            # Only compress if body is large enough
            if len(response_body) >= self.minimum_size:
                # Compress the response
                compressed_body = gzip.compress(response_body)
                
                # Create new response with compressed content
                return Response(
                    content=compressed_body,
                    status_code=response.status_code,
                    headers={
                        **dict(response.headers),
                        "content-encoding": "gzip",
                        "content-length": str(len(compressed_body)),
                    },
                    media_type=response.headers.get("content-type")
                )
            else:
                # Return original response for small payloads
                return Response(
                    content=response_body,
                    status_code=response.status_code,
                    headers=dict(response.headers),
                    media_type=response.headers.get("content-type")
                )
        
        return response


class LoggingMiddleware(BaseHTTPMiddleware):
    """Structured logging middleware with correlation IDs."""
    
    async def dispatch(self, request: Request, call_next):
        # Generate correlation ID
        request_id = str(uuid.uuid4())
        correlation_id.set(request_id)
        
        # Add correlation ID to request headers for downstream services
        request.headers.__dict__["_list"].append((b"x-correlation-id", request_id.encode()))
        
        start_time = time.time()
        client_ip = request.client.host if request.client else "unknown"
        
        # Log request
        logger.info(
            "Request started",
            extra={
                "correlation_id": request_id,
                "method": request.method,
                "url": str(request.url),
                "path": request.url.path,
                "query_params": dict(request.query_params),
                "client_ip": client_ip,
                "user_agent": request.headers.get("user-agent", "unknown"),
                "timestamp": time.time(),
                "event_type": "request_started"
            }
        )
        
        # Process request
        response = await call_next(request)
        
        # Calculate duration
        duration = time.time() - start_time
        
        # Add correlation ID to response headers
        response.headers["x-correlation-id"] = request_id
        
        # Log response
        logger.info(
            "Request completed",
            extra={
                "correlation_id": request_id,
                "method": request.method,
                "url": str(request.url),
                "path": request.url.path,
                "status_code": response.status_code,
                "duration_ms": round(duration * 1000, 2),
                "client_ip": client_ip,
                "timestamp": time.time(),
                "event_type": "request_completed"
            }
        )
        
        return response


class StructuredFormatter(logging.Formatter):
    """Custom formatter for structured JSON logging."""
    
    def format(self, record):
        # Base log structure
        log_entry = {
            "timestamp": self.formatTime(record, self.datefmt),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        
        # Add correlation ID if available
        try:
            request_id = correlation_id.get()
            if request_id:
                log_entry["correlation_id"] = request_id
        except LookupError:
            pass  # No correlation ID in context
        
        # Add extra fields from the record
        if hasattr(record, "correlation_id"):
            log_entry["correlation_id"] = record.correlation_id
        if hasattr(record, "method"):
            log_entry["method"] = record.method
        if hasattr(record, "url"):
            log_entry["url"] = record.url
        if hasattr(record, "path"):
            log_entry["path"] = record.path
        if hasattr(record, "status_code"):
            log_entry["status_code"] = record.status_code
        if hasattr(record, "duration_ms"):
            log_entry["duration_ms"] = record.duration_ms
        if hasattr(record, "client_ip"):
            log_entry["client_ip"] = record.client_ip
        if hasattr(record, "user_agent"):
            log_entry["user_agent"] = record.user_agent
        if hasattr(record, "query_params"):
            log_entry["query_params"] = record.query_params
        if hasattr(record, "event_type"):
            log_entry["event_type"] = record.event_type
        
        # Add exception info if present
        if record.exc_info:
            log_entry["exception"] = self.formatException(record.exc_info)
        
        return json.dumps(log_entry, ensure_ascii=False)


def setup_structured_logging():
    """Configure structured JSON logging."""
    formatter = StructuredFormatter()
    
    # Get root logger
    root_logger = logging.getLogger()
    
    # Clear existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Create console handler with JSON formatting
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    
    # Add handler to root logger
    root_logger.addHandler(console_handler)
    root_logger.setLevel(logging.INFO)
    
    return formatter