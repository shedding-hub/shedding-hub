"""
Comprehensive pytest test suite for the Shedding Hub API.
"""

import pytest
import asyncio
from httpx import AsyncClient, ASGITransport
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock
import yaml
import tempfile
import pathlib
import json
from datetime import datetime

# Import the main app
try:
    from main import app, load_all_datasets, get_data_dir
except ImportError:
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from main import app, load_all_datasets, get_data_dir


class TestApiEndpoints:
    """Test all API endpoints."""

    @pytest.fixture
    def client(self):
        """Create test client."""
        return TestClient(app)

    @pytest.fixture
    async def async_client(self):
        """Create async test client."""
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as ac:
            yield ac

    @pytest.fixture
    def mock_dataset(self):
        """Create a mock dataset for testing."""
        return {
            "title": "Test Dataset",
            "description": "A test dataset for unit testing",
            "doi": "10.1000/test",
            "analytes": {
                "test_analyte": {
                    "description": "Test analyte description",
                    "biomarker": "SARS-CoV-2",
                    "specimen": "stool",
                    "unit": "gc/mL",
                    "reference_event": "symptom onset",
                    "limit_of_detection": 100.0,
                    "limit_of_quantification": 1000.0
                }
            },
            "participants": [
                {
                    "attributes": {
                        "sex": "female",
                        "age": 30
                    },
                    "measurements": [
                        {
                            "time": 1.0,
                            "value": 5000.0,
                            "analyte": "test_analyte"
                        }
                    ]
                }
            ]
        }

    @pytest.fixture
    def mock_data_dir(self, mock_dataset):
        """Create temporary data directory with test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            data_dir = pathlib.Path(tmpdir)
            
            # Create test dataset file
            dataset_dir = data_dir / "test_dataset"
            dataset_dir.mkdir()
            
            dataset_file = dataset_dir / "test_dataset.yaml"
            with dataset_file.open('w') as f:
                yaml.dump(mock_dataset, f)
            
            # Create schema file
            schema_file = data_dir / ".schema.yaml"
            with schema_file.open('w') as f:
                yaml.dump({"type": "object"}, f)
            
            with patch('main.get_data_dir', return_value=data_dir):
                yield data_dir

    def test_root_endpoint(self, client):
        """Test root endpoint."""
        response = client.get("/")
        assert response.status_code == 200
        data = response.json()
        assert "message" in data
        assert "version" in data
        assert "endpoints" in data
        assert data["version"] == "0.1.0"

    def test_health_endpoint(self, client, mock_data_dir):
        """Test health check endpoint."""
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
        assert "datasets_available" in data
        assert "cache_status" in data
        assert data["data_directory_exists"] is True

    def test_stats_endpoint(self, client, mock_data_dir):
        """Test stats endpoint."""
        response = client.get("/stats")
        assert response.status_code == 200
        data = response.json()
        assert "total_datasets" in data
        assert "total_participants" in data
        assert "total_measurements" in data
        assert "unique_biomarkers" in data
        assert "unique_specimens" in data
        assert data["total_datasets"] >= 0

    def test_filters_endpoint(self, client, mock_data_dir):
        """Test filters endpoint."""
        response = client.get("/filters")
        assert response.status_code == 200
        data = response.json()
        assert "available_biomarkers" in data
        assert "available_specimens" in data
        assert isinstance(data["available_biomarkers"], list)
        assert isinstance(data["available_specimens"], list)

    def test_list_datasets(self, client, mock_data_dir):
        """Test dataset listing endpoint."""
        response = client.get("/datasets")
        assert response.status_code == 200
        data = response.json()
        assert "datasets" in data
        assert "total_found" in data
        assert "returned" in data
        assert isinstance(data["datasets"], list)

    def test_list_datasets_with_filters(self, client, mock_data_dir):
        """Test dataset listing with filters."""
        # Test biomarker filter
        response = client.get("/datasets?biomarker=SARS-CoV-2")
        assert response.status_code == 200
        
        # Test specimen filter
        response = client.get("/datasets?specimen=stool")
        assert response.status_code == 200
        
        # Test limit
        response = client.get("/datasets?limit=5")
        assert response.status_code == 200

    def test_get_dataset(self, client, mock_data_dir):
        """Test getting specific dataset."""
        # First get list of available datasets
        response = client.get("/datasets")
        datasets = response.json()["datasets"]
        
        if datasets:
            dataset_id = datasets[0]["dataset_id"]
            
            # Test getting full dataset
            response = client.get(f"/datasets/{dataset_id}")
            assert response.status_code == 200
            data = response.json()
            assert "dataset_id" in data
            assert "title" in data
            assert "participants" in data

    def test_get_dataset_metadata(self, client, mock_data_dir):
        """Test getting dataset metadata only."""
        response = client.get("/datasets")
        datasets = response.json()["datasets"]
        
        if datasets:
            dataset_id = datasets[0]["dataset_id"]
            
            response = client.get(f"/datasets/{dataset_id}/metadata")
            assert response.status_code == 200
            data = response.json()
            assert "dataset_id" in data
            assert "title" in data
            assert "participants_attributes" in data

    def test_get_nonexistent_dataset(self, client, mock_data_dir):
        """Test getting non-existent dataset."""
        response = client.get("/datasets/nonexistent_dataset")
        assert response.status_code == 404

    def test_search_datasets(self, client, mock_data_dir):
        """Test dataset search."""
        response = client.get("/search?q=test")
        assert response.status_code == 200
        data = response.json()
        assert "results" in data
        assert "query" in data
        assert "total_found" in data
        assert data["query"] == "test"

    def test_search_validation(self, client):
        """Test search parameter validation."""
        # Empty query should fail
        response = client.get("/search?q=")
        assert response.status_code == 422

    def test_refresh_cache(self, client, mock_data_dir):
        """Test cache refresh endpoint."""
        response = client.post("/refresh-cache")
        assert response.status_code == 200
        data = response.json()
        assert "message" in data
        assert "datasets_loaded" in data

    def test_parameter_validation(self, client):
        """Test parameter validation."""
        # Invalid limit (too high)
        response = client.get("/datasets?limit=10000")
        assert response.status_code == 422
        
        # Invalid limit (negative)
        response = client.get("/datasets?limit=-1")
        assert response.status_code == 422

    @pytest.mark.asyncio
    async def test_rate_limiting(self, mock_data_dir):
        """Test rate limiting functionality."""
        # Make multiple requests quickly to trigger rate limiting
        # Note: This test may be flaky depending on rate limit settings
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as client:
            tasks = []
            for i in range(20):
                tasks.append(client.get("/"))
            
            responses = await asyncio.gather(*tasks, return_exceptions=True)
            
            # Check that some requests succeeded
            success_count = sum(1 for r in responses if hasattr(r, 'status_code') and r.status_code == 200)
            assert success_count > 0


class TestDataLoading:
    """Test data loading functionality."""

    @pytest.mark.asyncio
    async def test_load_all_datasets_empty_dir(self):
        """Test loading datasets from empty directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch('main.get_data_dir', return_value=pathlib.Path(tmpdir)):
                datasets = await load_all_datasets(force_reload=True)
                assert isinstance(datasets, dict)
                assert len(datasets) == 0

    @pytest.mark.asyncio
    async def test_load_all_datasets_with_invalid_yaml(self):
        """Test loading datasets with invalid YAML."""
        with tempfile.TemporaryDirectory() as tmpdir:
            data_dir = pathlib.Path(tmpdir)
            
            # Create invalid YAML file
            dataset_dir = data_dir / "invalid_dataset"
            dataset_dir.mkdir()
            
            dataset_file = dataset_dir / "invalid_dataset.yaml"
            with dataset_file.open('w') as f:
                f.write("invalid: yaml: content: !!!")
            
            with patch('main.get_data_dir', return_value=data_dir):
                datasets = await load_all_datasets(force_reload=True)
                # Should skip invalid files and return empty dict
                assert isinstance(datasets, dict)

    def test_get_data_dir(self):
        """Test get_data_dir function."""
        data_dir = get_data_dir()
        assert isinstance(data_dir, pathlib.Path)
        # Should point to the data directory relative to the API directory
        assert "data" in str(data_dir)


class TestErrorHandling:
    """Test error handling scenarios."""

    @pytest.fixture
    def client(self):
        return TestClient(app)

    def test_missing_data_directory(self, client):
        """Test behavior when data directory is missing."""
        with patch('main.get_data_dir', return_value=pathlib.Path("/nonexistent/path")):
            response = client.get("/datasets")
            assert response.status_code == 500

    def test_file_permission_error(self, client):
        """Test handling of file permission errors."""
        # This test would require creating files with restricted permissions
        # Implementation depends on the testing environment
        pass

    def test_yaml_parsing_error(self, client):
        """Test handling of YAML parsing errors."""
        with tempfile.TemporaryDirectory() as tmpdir:
            data_dir = pathlib.Path(tmpdir)
            dataset_dir = data_dir / "bad_dataset"
            dataset_dir.mkdir()
            
            # Create malformed YAML
            dataset_file = dataset_dir / "bad_dataset.yaml"
            with dataset_file.open('w') as f:
                f.write("[invalid yaml content")
            
            with patch('main.get_data_dir', return_value=data_dir):
                response = client.get("/datasets")
                # Should handle the error gracefully
                assert response.status_code in [200, 500]


class TestModelValidation:
    """Test Pydantic model validation."""

    def test_dataset_summary_validation(self):
        """Test DatasetSummary model validation."""
        try:
            from models import DatasetSummary, DatasetStats
        except ImportError:
            from .models import DatasetSummary, DatasetStats
        
        valid_data = {
            "dataset_id": "test_dataset",
            "title": "Test Dataset",
            "description": "Test description",
            "stats": {
                "participants_count": 10,
                "analytes_count": 2,
                "total_measurements": 50
            },
            "analytes": ["analyte1", "analyte2"],
            "last_modified": datetime.now().isoformat()
        }
        
        # Should not raise validation error
        summary = DatasetSummary(**valid_data)
        assert summary.dataset_id == "test_dataset"
        assert summary.stats.participants_count == 10

    def test_search_result_validation(self):
        """Test SearchResult model validation."""
        try:
            from models import SearchResult, DatasetStats
        except ImportError:
            from .models import SearchResult, DatasetStats
        
        valid_data = {
            "dataset_id": "test",
            "title": "Test",
            "description": "Test description",
            "match_score": 5,
            "stats": {
                "participants_count": 1,
                "analytes_count": 1,
                "total_measurements": 1
            }
        }
        
        result = SearchResult(**valid_data)
        assert result.match_score == 5


# Test configuration
@pytest.fixture(scope="session", autouse=True)
def setup_test_environment():
    """Set up test environment."""
    # Any global test setup can go here
    pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])