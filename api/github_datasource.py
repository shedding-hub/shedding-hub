"""
GitHub data source for fetching datasets from the shedding-hub repository.

This module provides functionality to fetch dataset files directly from GitHub,
enabling the API to always serve the most up-to-date data without requiring
local file storage.
"""

import asyncio
import base64
import yaml
from typing import Dict, List, Any, Optional
from datetime import datetime, timedelta
import httpx
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

class GitHubDataSource:
    """GitHub data source for fetching shedding-hub datasets."""
    
    def __init__(self, 
                 repo_owner: str = "shedding-hub",
                 repo_name: str = "shedding-hub", 
                 branch: str = "main",
                 data_path: str = "data",
                 github_token: Optional[str] = None):
        """Initialize GitHub data source.
        
        Args:
            repo_owner: GitHub repository owner
            repo_name: GitHub repository name  
            branch: Git branch to fetch from
            data_path: Path to data directory in repo
            github_token: Optional GitHub token for higher rate limits
        """
        self.repo_owner = repo_owner
        self.repo_name = repo_name
        self.branch = branch
        self.data_path = data_path
        self.github_token = github_token
        
        self.base_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}"
        self.headers = {
            "Accept": "application/vnd.github.v3+json",
            "User-Agent": "ShedingHub-API/1.0"
        }
        
        if github_token:
            self.headers["Authorization"] = f"token {github_token}"
            
        # Cache for dataset listings and content
        self._dataset_list_cache: Optional[List[str]] = None
        self._dataset_content_cache: Dict[str, Dict] = {}
        self._cache_timestamp: Optional[datetime] = None
        self._cache_ttl = timedelta(minutes=10)  # Cache for 10 minutes
    
    async def _make_request(self, url: str, client: httpx.AsyncClient) -> Dict[str, Any]:
        """Make async HTTP request to GitHub API."""
        try:
            response = await client.get(url, headers=self.headers)
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            logger.error(f"GitHub API error: {e.response.status_code} for {url}")
            raise
        except Exception as e:
            logger.error(f"Request error for {url}: {e}")
            raise
    
    async def _get_directory_contents(self, path: str = "") -> List[Dict[str, Any]]:
        """Get contents of a directory from GitHub."""
        full_path = f"{self.data_path}/{path}".strip("/") if path else self.data_path
        url = f"{self.base_url}/contents/{full_path}?ref={self.branch}"
        
        async with httpx.AsyncClient(timeout=30.0) as client:
            return await self._make_request(url, client)
    
    async def _get_file_content(self, file_path: str) -> str:
        """Get file content from GitHub (base64 decoded)."""
        url = f"{self.base_url}/contents/{file_path}?ref={self.branch}"
        
        async with httpx.AsyncClient(timeout=30.0) as client:
            file_data = await self._make_request(url, client)
            
            if file_data.get("encoding") == "base64":
                content = base64.b64decode(file_data["content"]).decode("utf-8")
                return content
            else:
                return file_data.get("content", "")
    
    def _is_cache_expired(self) -> bool:
        """Check if cache has expired."""
        if self._cache_timestamp is None:
            return True
        return datetime.now() - self._cache_timestamp > self._cache_ttl
    
    async def get_dataset_list(self, force_refresh: bool = False) -> List[str]:
        """Get list of available datasets."""
        if not force_refresh and not self._is_cache_expired() and self._dataset_list_cache:
            logger.info("Using cached dataset list")
            return self._dataset_list_cache
        
        logger.info("Fetching dataset list from GitHub")
        
        try:
            contents = await self._get_directory_contents()
            
            # Filter for directories (datasets)
            datasets = [
                item["name"] for item in contents 
                if item["type"] == "dir" and not item["name"].startswith(".")
            ]
            
            self._dataset_list_cache = sorted(datasets)
            self._cache_timestamp = datetime.now()
            
            logger.info(f"Found {len(datasets)} datasets in GitHub repository")
            return self._dataset_list_cache
            
        except Exception as e:
            logger.error(f"Error fetching dataset list: {e}")
            # Return cached data if available
            if self._dataset_list_cache:
                logger.warning("Returning stale cached dataset list")
                return self._dataset_list_cache
            raise
    
    async def get_dataset(self, dataset_id: str, force_refresh: bool = False) -> Dict[str, Any]:
        """Get full dataset content."""
        cache_key = f"{dataset_id}_full"
        
        if (not force_refresh and 
            cache_key in self._dataset_content_cache and 
            not self._is_cache_expired()):
            logger.info(f"Using cached dataset: {dataset_id}")
            return self._dataset_content_cache[cache_key]
        
        logger.info(f"Fetching dataset from GitHub: {dataset_id}")
        
        try:
            # Get the YAML file content
            yaml_path = f"{self.data_path}/{dataset_id}/{dataset_id}.yaml"
            yaml_content = await self._get_file_content(yaml_path)
            
            # Parse YAML
            dataset = yaml.safe_load(yaml_content)
            
            # Add metadata
            dataset["dataset_id"] = dataset_id
            dataset["last_modified"] = datetime.now().isoformat()
            dataset["file_path"] = yaml_path
            
            # Calculate and add stats
            dataset["stats"] = self._calculate_dataset_stats(dataset)
            
            # Cache the result
            self._dataset_content_cache[cache_key] = dataset
            
            logger.info(f"Successfully loaded dataset: {dataset_id}")
            return dataset
            
        except Exception as e:
            logger.error(f"Error fetching dataset {dataset_id}: {e}")
            # Return cached data if available
            if cache_key in self._dataset_content_cache:
                logger.warning(f"Returning stale cached dataset: {dataset_id}")
                return self._dataset_content_cache[cache_key]
            raise
    
    async def get_dataset_metadata(self, dataset_id: str, force_refresh: bool = False) -> Dict[str, Any]:
        """Get dataset metadata without full participant data."""
        # For now, we'll get the full dataset and strip measurements
        # In future, could optimize to only fetch metadata
        dataset = await self.get_dataset(dataset_id, force_refresh)
        
        # Remove measurements to create metadata-only version
        metadata = dataset.copy()
        
        if "participants" in metadata:
            participants_metadata = []
            for participant in metadata["participants"]:
                participant_meta = participant.copy()
                participant_meta["measurements"] = []  # Remove measurements
                participants_metadata.append(participant_meta)
            
            metadata["participants_attributes"] = [
                p.get("attributes", {}) for p in metadata["participants"]
            ]
            del metadata["participants"]  # Remove full participants
        
        return metadata
    
    async def get_repository_info(self) -> Dict[str, Any]:
        """Get repository information."""
        url = self.base_url
        
        async with httpx.AsyncClient(timeout=30.0) as client:
            repo_info = await self._make_request(url, client)
            
            return {
                "repository": f"{self.repo_owner}/{self.repo_name}",
                "description": repo_info.get("description", ""),
                "url": repo_info.get("html_url", ""),
                "last_updated": repo_info.get("updated_at", ""),
                "stars": repo_info.get("stargazers_count", 0),
                "forks": repo_info.get("forks_count", 0),
                "branch": self.branch,
                "data_path": self.data_path
            }
    
    def _calculate_dataset_stats(self, dataset: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate dataset statistics."""
        participants = dataset.get("participants", [])
        analytes = dataset.get("analytes", {})
        
        participants_count = len(participants)
        analytes_count = len(analytes)
        
        # Count total measurements
        total_measurements = 0
        for participant in participants:
            measurements = participant.get("measurements", [])
            total_measurements += len(measurements)
        
        return {
            "participants_count": participants_count,
            "analytes_count": analytes_count,
            "total_measurements": total_measurements
        }
    
    async def refresh_cache(self) -> Dict[str, Any]:
        """Force refresh all cached data."""
        logger.info("Refreshing GitHub data cache")
        
        # Clear existing cache
        self._dataset_list_cache = None
        self._dataset_content_cache.clear()
        self._cache_timestamp = None
        
        # Reload dataset list
        datasets = await self.get_dataset_list(force_refresh=True)
        
        return {
            "message": "Cache refreshed successfully",
            "datasets_found": len(datasets),
            "cache_timestamp": datetime.now().isoformat(),
            "source": "GitHub Repository"
        }

# Global instance
github_datasource = GitHubDataSource()

async def get_github_datasets() -> List[str]:
    """Get list of datasets from GitHub."""
    return await github_datasource.get_dataset_list()

async def get_github_dataset(dataset_id: str) -> Dict[str, Any]:
    """Get full dataset from GitHub."""
    return await github_datasource.get_dataset(dataset_id)

async def get_github_dataset_metadata(dataset_id: str) -> Dict[str, Any]:
    """Get dataset metadata from GitHub."""
    return await github_datasource.get_dataset_metadata(dataset_id)