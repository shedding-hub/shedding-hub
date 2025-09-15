"""
Data source abstraction layer for the Shedding Hub API.

This module provides a unified interface for fetching data from different sources
(local files or GitHub repository), allowing the API to seamlessly switch between
data sources based on configuration.
"""

import asyncio
from typing import Dict, List, Any, Optional
from datetime import datetime
import logging

from config import settings

logger = logging.getLogger(__name__)

class DataSourceManager:
    """Manages different data sources for the API."""
    
    def __init__(self):
        self._github_datasource = None
        self._local_datasource = None
        
    async def _get_github_datasource(self):
        """Get GitHub data source instance."""
        if self._github_datasource is None:
            from github_datasource import GitHubDataSource
            self._github_datasource = GitHubDataSource(
                repo_owner=settings.github_repo_owner,
                repo_name=settings.github_repo_name,
                branch=settings.github_branch,
                data_path=settings.github_data_path,
                github_token=settings.github_token
            )
            # Set cache TTL from settings
            from datetime import timedelta
            self._github_datasource._cache_ttl = timedelta(seconds=settings.github_cache_ttl)
        return self._github_datasource
    
    async def _get_local_datasource(self):
        """Get local file data source instance."""
        if self._local_datasource is None:
            from async_utils import load_all_datasets_async
            self._local_datasource = load_all_datasets_async
        return self._local_datasource
    
    async def get_dataset_list(self, force_refresh: bool = False) -> List[str]:
        """Get list of available datasets from configured source."""
        if settings.github_enabled:
            logger.info("Using GitHub data source")
            github_source = await self._get_github_datasource()
            datasets_dict = await github_source.get_dataset_list(force_refresh)
            return datasets_dict
        else:
            logger.info("Using local file data source")
            # For local files, we'll use the existing load_all_datasets logic
            # but just return the keys
            try:
                from main import load_all_datasets  # Import from main module
                datasets_dict = await load_all_datasets(force_refresh)
                return list(datasets_dict.keys())
            except Exception as e:
                logger.error(f"Error loading local datasets: {e}")
                raise
    
    async def get_dataset(self, dataset_id: str, force_refresh: bool = False) -> Dict[str, Any]:
        """Get full dataset from configured source."""
        if settings.github_enabled:
            logger.info(f"Fetching dataset {dataset_id} from GitHub")
            github_source = await self._get_github_datasource()
            return await github_source.get_dataset(dataset_id, force_refresh)
        else:
            logger.info(f"Fetching dataset {dataset_id} from local files")
            try:
                from main import load_all_datasets
                datasets_dict = await load_all_datasets(force_refresh)
                if dataset_id not in datasets_dict:
                    raise KeyError(f"Dataset {dataset_id} not found")
                return datasets_dict[dataset_id]
            except Exception as e:
                logger.error(f"Error loading local dataset {dataset_id}: {e}")
                raise
    
    async def get_dataset_metadata(self, dataset_id: str, force_refresh: bool = False) -> Dict[str, Any]:
        """Get dataset metadata from configured source."""
        if settings.github_enabled:
            logger.info(f"Fetching metadata for {dataset_id} from GitHub")
            github_source = await self._get_github_datasource()
            return await github_source.get_dataset_metadata(dataset_id, force_refresh)
        else:
            logger.info(f"Fetching metadata for {dataset_id} from local files")
            # Get full dataset and convert to metadata format
            dataset = await self.get_dataset(dataset_id, force_refresh)
            
            # Create metadata version (without full measurements)
            metadata = dataset.copy()
            
            if "participants" in metadata:
                # Extract participant attributes only
                metadata["participants_attributes"] = [
                    p.get("attributes", {}) for p in metadata["participants"]
                ]
                # Remove full participants data
                del metadata["participants"]
            
            return metadata
    
    async def get_data_source_info(self) -> Dict[str, Any]:
        """Get information about the current data source."""
        if settings.github_enabled:
            github_source = await self._get_github_datasource()
            github_info = await github_source.get_repository_info()
            return {
                "source_type": "GitHub Repository",
                "github_info": github_info,
                "cache_ttl_seconds": settings.github_cache_ttl
            }
        else:
            return {
                "source_type": "Local Files",
                "data_directory": str(settings.get_data_directory()),
                "cache_ttl_seconds": settings.cache_ttl_seconds
            }
    
    async def refresh_cache(self) -> Dict[str, Any]:
        """Refresh data source cache."""
        if settings.github_enabled:
            logger.info("Refreshing GitHub cache")
            github_source = await self._get_github_datasource()
            return await github_source.refresh_cache()
        else:
            logger.info("Refreshing local file cache")
            try:
                from main import load_all_datasets
                datasets = await load_all_datasets(force_reload=True)
                return {
                    "message": "Local cache refreshed successfully", 
                    "datasets_loaded": len(datasets),
                    "cache_timestamp": datetime.now().isoformat(),
                    "source": "Local Files"
                }
            except Exception as e:
                logger.error(f"Error refreshing local cache: {e}")
                raise

# Global data source manager instance
data_source_manager = DataSourceManager()