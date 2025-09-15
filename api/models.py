"""
Pydantic models for the Shedding Hub API.

These models provide request/response validation, serialization,
and automatic API documentation generation.
"""

from pydantic import BaseModel, Field, ConfigDict
from typing import List, Dict, Any, Optional, Union
from datetime import datetime
from enum import Enum


class SpecimenType(str, Enum):
    """Valid specimen types from the schema."""
    stool = "stool"
    sputum = "sputum"
    urine = "urine"
    plasma = "plasma"
    oropharyngeal_swab = "oropharyngeal_swab"
    nasopharyngeal_swab = "nasopharyngeal_swab"
    anterior_nares_swab = "anterior_nares_swab"
    saliva = "saliva"
    rectal_swab = "rectal_swab"
    unknown = "unknown"
    bronchoalveolar_lavage_fluid = "bronchoalveolar_lavage_fluid"
    serum = "serum"


class BiomarkerType(str, Enum):
    """Valid biomarker types from the schema."""
    sars_cov_2 = "SARS-CoV-2"
    mtdna = "mtDNA"
    pmmov = "PMMoV"
    crassphage = "crAssphage"
    influenza = "influenza"
    sapovirus = "sapovirus"
    sars = "SARS"


class UnitType(str, Enum):
    """Valid measurement units from the schema."""
    gc_dry_gram = "gc/dry gram"
    gc_ml = "gc/mL"
    gc_swab = "gc/swab"
    gc_wet_gram = "gc/wet gram"
    cycle_threshold = "cycle threshold"
    pfu_ml = "pfu/mL"


class ReferenceEventType(str, Enum):
    """Valid reference events from the schema."""
    symptom_onset = "symptom onset"
    confirmation_date = "confirmation date"
    enrollment = "enrollment"
    hospital_admission = "hospital admission"


# Base response model
class BaseResponse(BaseModel):
    """Base response model with common metadata."""
    model_config = ConfigDict(extra="allow")


# Request models
class DatasetListRequest(BaseModel):
    """Request parameters for listing datasets."""
    biomarker: Optional[BiomarkerType] = Field(None, description="Filter by biomarker type")
    specimen: Optional[SpecimenType] = Field(None, description="Filter by specimen type")
    limit: int = Field(100, ge=1, le=1000, description="Maximum number of datasets to return")


class SearchRequest(BaseModel):
    """Request parameters for searching datasets."""
    q: str = Field(..., min_length=1, description="Search query string")
    fields: List[str] = Field(["title", "description"], description="Fields to search in")
    limit: int = Field(20, ge=1, le=1000, description="Maximum number of results to return")


# Data models
class DatasetStats(BaseModel):
    """Statistics for a dataset."""
    participants_count: int = Field(..., description="Number of participants")
    analytes_count: int = Field(..., description="Number of analytes")
    total_measurements: int = Field(..., description="Total number of measurements")


class AnalyteSpecification(BaseModel):
    """Analyte specification model."""
    description: str = Field(..., description="Description of the analyte")
    biomarker: BiomarkerType = Field(..., description="Biomarker being measured")
    specimen: Union[SpecimenType, List[SpecimenType]] = Field(..., description="Specimen type(s)")
    unit: UnitType = Field(..., description="Measurement unit")
    reference_event: ReferenceEventType = Field(..., description="Reference event for time measurements")
    gene_target: Optional[str] = Field(None, description="Gene target for genomic biomarkers")
    limit_of_detection: Union[float, str] = Field(..., description="Limit of detection")
    limit_of_quantification: Union[float, str] = Field(..., description="Limit of quantification")
    limit_of_blank: Optional[Union[float, str]] = Field(None, description="Limit of blank")


class ParticipantAttributes(BaseModel):
    """Participant attributes model."""
    sex: Optional[str] = Field(None, description="Sex at birth")
    age: Optional[Union[int, str]] = Field(None, description="Age in years")
    race: Optional[str] = Field(None, description="Race")
    ethnicity: Optional[str] = Field(None, description="Ethnicity")
    vaccinated: Optional[Union[bool, str]] = Field(None, description="Vaccination status")
    lineage: Optional[str] = Field(None, description="Pathogen lineage")
    variant: Optional[str] = Field(None, description="Pathogen variant")


class Measurement(BaseModel):
    """Individual measurement model."""
    time: Union[float, str] = Field(..., description="Temporal offset")
    value: Union[float, str] = Field(..., description="Measurement value")
    analyte: str = Field(..., description="Analyte identifier")
    sample_id: Optional[Union[str, int]] = Field(None, description="Sample identifier")
    limit_of_quantification: Optional[Union[float, str]] = Field(None, description="Sample-specific LOQ")
    limit_of_blank: Optional[Union[float, str]] = Field(None, description="Sample-specific LOB")


class Participant(BaseModel):
    """Participant model."""
    attributes: Optional[ParticipantAttributes] = Field(None, description="Participant attributes")
    measurements: List[Measurement] = Field(..., description="List of measurements")


# Response models
class DatasetSummary(BaseModel):
    """Summary information for a dataset in list responses."""
    dataset_id: str = Field(..., description="Unique dataset identifier")
    title: str = Field(..., description="Dataset title")
    doi: Optional[str] = Field(None, description="Digital object identifier")
    url: Optional[str] = Field(None, description="Dataset URL")
    description: str = Field(..., description="Truncated description")
    stats: DatasetStats = Field(..., description="Dataset statistics")
    analytes: List[str] = Field(..., description="List of analyte identifiers")
    last_modified: str = Field(..., description="Last modification timestamp")


class DatasetListResponse(BaseResponse):
    """Response for dataset listing."""
    total_found: int = Field(..., description="Total datasets matching filters")
    returned: int = Field(..., description="Number of datasets returned")
    filters_applied: Dict[str, Any] = Field(..., description="Applied filters")
    datasets: List[DatasetSummary] = Field(..., description="List of dataset summaries")


class FullDataset(BaseModel):
    """Complete dataset model."""
    dataset_id: str = Field(..., description="Unique dataset identifier")
    title: str = Field(..., description="Dataset title")
    doi: Optional[str] = Field(None, description="Digital object identifier")
    url: Optional[str] = Field(None, description="Dataset URL")
    description: str = Field(..., description="Full description")
    analytes: Dict[str, AnalyteSpecification] = Field(..., description="Analyte specifications")
    participants: List[Participant] = Field(..., description="Participant data")
    stats: DatasetStats = Field(..., description="Dataset statistics")
    last_modified: str = Field(..., description="Last modification timestamp")
    file_path: str = Field(..., description="Source file path")


class DatasetMetadata(BaseModel):
    """Dataset metadata without measurements."""
    dataset_id: str = Field(..., description="Unique dataset identifier")
    title: str = Field(..., description="Dataset title")
    doi: Optional[str] = Field(None, description="Digital object identifier")
    url: Optional[str] = Field(None, description="Dataset URL")
    description: str = Field(..., description="Full description")
    analytes: Dict[str, AnalyteSpecification] = Field(..., description="Analyte specifications")
    stats: DatasetStats = Field(..., description="Dataset statistics")
    last_modified: str = Field(..., description="Last modification timestamp")
    participants_attributes: List[ParticipantAttributes] = Field(..., description="Participant attributes only")


class SearchResult(BaseModel):
    """Single search result."""
    dataset_id: str = Field(..., description="Dataset identifier")
    title: str = Field(..., description="Dataset title")
    description: str = Field(..., description="Truncated description")
    match_score: int = Field(..., description="Search relevance score")
    stats: DatasetStats = Field(..., description="Dataset statistics")


class SearchResponse(BaseResponse):
    """Response for search queries."""
    query: str = Field(..., description="Original search query")
    fields_searched: List[str] = Field(..., description="Fields that were searched")
    total_found: int = Field(..., description="Total results found")
    returned: int = Field(..., description="Number of results returned")
    results: List[SearchResult] = Field(..., description="Search results")


class ProjectStats(BaseResponse):
    """Overall project statistics."""
    total_datasets: int = Field(..., description="Total number of datasets")
    total_participants: int = Field(..., description="Total number of participants")
    total_measurements: int = Field(..., description="Total number of measurements")
    unique_biomarkers: List[str] = Field(..., description="List of unique biomarkers")
    unique_specimens: List[str] = Field(..., description="List of unique specimen types")
    biomarker_count: int = Field(..., description="Number of unique biomarkers")
    specimen_count: int = Field(..., description="Number of unique specimen types")
    cache_info: Dict[str, Any] = Field(..., description="Cache status information")


class HealthStatus(BaseResponse):
    """Health check response."""
    status: str = Field(..., description="Overall health status")
    datasets_available: int = Field(..., description="Number of available datasets")
    cache_status: str = Field(..., description="Cache status")
    data_directory_exists: bool = Field(..., description="Whether data directory exists")
    timestamp: str = Field(..., description="Health check timestamp")


class AvailableFilters(BaseResponse):
    """Available filter options."""
    available_biomarkers: List[str] = Field(..., description="Available biomarker types")
    available_specimens: List[str] = Field(..., description="Available specimen types")


class CacheRefreshResponse(BaseResponse):
    """Cache refresh response."""
    message: str = Field(..., description="Status message")
    datasets_loaded: int = Field(..., description="Number of datasets loaded")
    cache_timestamp: str = Field(..., description="Cache update timestamp")


class ApiInfo(BaseResponse):
    """Root endpoint API information."""
    message: str = Field(..., description="Welcome message")
    version: str = Field(..., description="API version")
    description: str = Field(..., description="API description")
    endpoints: Dict[str, str] = Field(..., description="Available endpoints")
    cache_info: Dict[str, Any] = Field(..., description="Cache information")


class ErrorDetail(BaseModel):
    """Error response model."""
    detail: str = Field(..., description="Error message")
    error_code: Optional[str] = Field(None, description="Error code")
    timestamp: Optional[str] = Field(None, description="Error timestamp")


class ExportFormat(str, Enum):
    """Supported export formats."""
    csv = "csv"
    json = "json"
    excel = "excel"


class ExportResponse(BaseResponse):
    """Dataset export response."""
    dataset_id: str = Field(..., description="Dataset ID that was exported")
    format: ExportFormat = Field(..., description="Export format")
    filename: str = Field(..., description="Generated filename")
    size_bytes: int = Field(..., description="File size in bytes")
    export_timestamp: str = Field(..., description="Export timestamp")


class MetricsSummary(BaseResponse):
    """API metrics summary."""
    memory_usage_mb: float = Field(..., description="Memory usage in MB")
    cpu_usage_percent: float = Field(..., description="CPU usage percentage")
    uptime_seconds: int = Field(..., description="API uptime in seconds")
    datasets_loaded: int = Field(..., description="Number of datasets loaded")
    total_requests: int = Field(..., description="Total API requests")
    cache_hit_rate: float = Field(..., description="Cache hit rate (0-1)")
    error_rate: float = Field(..., description="Error rate (0-1)")


class DetailedMetrics(BaseResponse):
    """Detailed API metrics for monitoring."""
    requests: Dict[str, Any] = Field(..., description="Request statistics")
    performance: Dict[str, Any] = Field(..., description="Performance metrics")
    system: Dict[str, Any] = Field(..., description="System resource metrics")
    rate_limiting: Dict[str, Any] = Field(..., description="Rate limiting statistics")