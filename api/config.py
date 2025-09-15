"""
Configuration management for the Shedding Hub API.

Handles environment variables, settings, and configuration validation.
"""

import os
from typing import List, Optional
from pydantic_settings import BaseSettings
from pydantic import Field, field_validator, ConfigDict
from pathlib import Path


class Settings(BaseSettings):
    """Application settings with environment variable support."""
    
    # API Configuration
    api_title: str = Field("Shedding Hub API", description="API title")
    api_description: str = Field(
        "API for accessing biomarker shedding datasets from the Shedding Hub project",
        description="API description"
    )
    api_version: str = Field("0.1.0", description="API version")
    debug: bool = Field(False, description="Enable debug mode")
    
    # Server Configuration
    host: str = Field("0.0.0.0", description="Server host")
    port: int = Field(8000, description="Server port")
    reload: bool = Field(False, description="Auto-reload on code changes")
    
    # CORS Configuration
    cors_origins: List[str] = Field(["*"], description="Allowed CORS origins")
    cors_credentials: bool = Field(True, description="Allow CORS credentials")
    cors_methods: List[str] = Field(["*"], description="Allowed CORS methods")
    cors_headers: List[str] = Field(["*"], description="Allowed CORS headers")
    
    # Cache Configuration
    cache_ttl_seconds: int = Field(3600, description="Cache TTL in seconds (1 hour)")
    cache_max_size: int = Field(1000, description="Maximum cache size")
    
    # Rate Limiting Configuration
    rate_limit_enabled: bool = Field(True, description="Enable rate limiting")
    rate_limit_default: str = Field("100/minute", description="Default rate limit")
    rate_limit_search: str = Field("30/minute", description="Search endpoint rate limit")
    rate_limit_expensive: str = Field("10/minute", description="Expensive operations rate limit")
    rate_limit_cache_ops: str = Field("5/minute", description="Cache operations rate limit")
    
    # Data Configuration
    data_directory: Optional[Path] = Field(None, description="Custom data directory path")
    max_datasets_limit: int = Field(1000, description="Maximum datasets per request")
    default_datasets_limit: int = Field(100, description="Default datasets per request")
    default_search_limit: int = Field(20, description="Default search results limit")
    
    # GitHub Integration
    github_enabled: bool = Field(False, description="Enable GitHub data source")
    github_repo_owner: str = Field("shedding-hub", description="GitHub repository owner")
    github_repo_name: str = Field("shedding-hub", description="GitHub repository name")  
    github_branch: str = Field("main", description="GitHub branch to fetch from")
    github_data_path: str = Field("data", description="Path to data directory in GitHub repo")
    github_token: Optional[str] = Field(None, description="GitHub token for higher rate limits")
    github_cache_ttl: int = Field(600, description="GitHub cache TTL in seconds (10 minutes)")
    
    # Logging Configuration
    log_level: str = Field("INFO", description="Logging level")
    log_format: str = Field(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        description="Log format string"
    )
    log_file: Optional[Path] = Field(None, description="Log file path")
    
    # Security Configuration
    allowed_hosts: List[str] = Field(["*"], description="Allowed host headers")
    max_request_size: int = Field(16777216, description="Max request size in bytes (16MB)")
    
    # Performance Configuration
    worker_processes: int = Field(1, description="Number of worker processes")
    max_connections: int = Field(1000, description="Maximum concurrent connections")
    keepalive_timeout: int = Field(5, description="Keep-alive timeout in seconds")
    
    @field_validator('cors_origins', mode='before')
    @classmethod
    def parse_cors_origins(cls, v):
        """Parse CORS origins from string or list."""
        if isinstance(v, str):
            return [origin.strip() for origin in v.split(',')]
        return v
    
    @field_validator('cors_methods', mode='before')
    @classmethod
    def parse_cors_methods(cls, v):
        """Parse CORS methods from string or list."""
        if isinstance(v, str):
            return [method.strip() for method in v.split(',')]
        return v
    
    @field_validator('cors_headers', mode='before')
    @classmethod
    def parse_cors_headers(cls, v):
        """Parse CORS headers from string or list."""
        if isinstance(v, str):
            return [header.strip() for header in v.split(',')]
        return v
    
    @field_validator('allowed_hosts', mode='before')
    @classmethod
    def parse_allowed_hosts(cls, v):
        """Parse allowed hosts from string or list."""
        if isinstance(v, str):
            return [host.strip() for host in v.split(',')]
        return v
    
    @field_validator('data_directory', mode='before')
    @classmethod
    def parse_data_directory(cls, v):
        """Parse data directory path."""
        if v is None:
            return None
        return Path(v)
    
    @field_validator('log_file', mode='before')
    @classmethod
    def parse_log_file(cls, v):
        """Parse log file path."""
        if v is None:
            return None
        return Path(v)
    
    @field_validator('port')
    @classmethod
    def validate_port(cls, v):
        """Validate port number."""
        if not (1 <= v <= 65535):
            raise ValueError('Port must be between 1 and 65535')
        return v
    
    @field_validator('cache_ttl_seconds')
    @classmethod
    def validate_cache_ttl(cls, v):
        """Validate cache TTL."""
        if v < 0:
            raise ValueError('Cache TTL must be non-negative')
        return v
    
    @field_validator('log_level')
    @classmethod
    def validate_log_level(cls, v):
        """Validate log level."""
        valid_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
        if v.upper() not in valid_levels:
            raise ValueError(f'Log level must be one of: {valid_levels}')
        return v.upper()
    
    def get_data_directory(self) -> Path:
        """Get the data directory path."""
        if self.data_directory:
            return self.data_directory
        
        # Default to ../data relative to the API directory
        api_dir = Path(__file__).parent
        return api_dir.parent / "data"
    
    def is_production(self) -> bool:
        """Check if running in production mode."""
        return not self.debug and not self.reload
    
    def get_cors_config(self) -> dict:
        """Get CORS configuration dictionary."""
        return {
            "allow_origins": self.cors_origins,
            "allow_credentials": self.cors_credentials,
            "allow_methods": self.cors_methods,
            "allow_headers": self.cors_headers,
        }
    
    def get_rate_limits(self) -> dict:
        """Get rate limiting configuration."""
        return {
            "default": self.rate_limit_default,
            "search": self.rate_limit_search,
            "expensive": self.rate_limit_expensive,
            "cache_ops": self.rate_limit_cache_ops,
        }

    model_config = ConfigDict(
        env_file=".env",
        env_file_encoding="utf-8", 
        case_sensitive=False,
        env_prefix="SHEDDING_HUB_"
    )


class DevelopmentSettings(Settings):
    """Development-specific settings."""
    debug: bool = True
    reload: bool = True
    log_level: str = "DEBUG"
    cors_origins: List[str] = ["http://localhost:3000", "http://127.0.0.1:3000"]


class ProductionSettings(Settings):
    """Production-specific settings."""
    debug: bool = False
    reload: bool = False
    log_level: str = "INFO"
    cors_origins: List[str] = []  # Must be explicitly set in production
    
    @field_validator('cors_origins')
    @classmethod
    def validate_production_cors(cls, v):
        """Validate CORS origins in production."""
        if not v or "*" in v:
            raise ValueError("CORS origins must be explicitly set in production")
        return v


def get_settings() -> Settings:
    """Get application settings based on environment."""
    environment = os.getenv("ENVIRONMENT", "development").lower()
    
    if environment == "production":
        return ProductionSettings()
    elif environment == "development":
        return DevelopmentSettings()
    else:
        return Settings()


# Create global settings instance
settings = get_settings()