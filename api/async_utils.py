"""
Async utilities for the Shedding Hub API.

Provides async file operations and other performance optimizations.
"""

import aiofiles
import asyncio
import yaml
import pathlib
from typing import Dict, Any, List
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


async def read_yaml_file_async(file_path: pathlib.Path) -> Dict[str, Any]:
    """
    Asynchronously read and parse a YAML file.
    
    Args:
        file_path: Path to the YAML file
        
    Returns:
        Parsed YAML content as dictionary
        
    Raises:
        FileNotFoundError: If file doesn't exist
        yaml.YAMLError: If YAML parsing fails
    """
    try:
        async with aiofiles.open(file_path, 'r', encoding='utf-8') as f:
            content = await f.read()
            return yaml.safe_load(content)
    except Exception as e:
        logger.error(f"Error reading YAML file {file_path}: {e}")
        raise


async def load_dataset_async(file_path: pathlib.Path, dataset_id: str) -> Dict[str, Any]:
    """
    Asynchronously load a single dataset with metadata.
    
    Args:
        file_path: Path to the dataset YAML file
        dataset_id: Dataset identifier
        
    Returns:
        Dataset dictionary with added metadata
    """
    try:
        # Read the YAML file asynchronously
        dataset = await read_yaml_file_async(file_path)
        
        # Add metadata
        dataset['dataset_id'] = dataset_id
        dataset['file_path'] = str(file_path)
        
        # Get file stats asynchronously
        stat_info = file_path.stat()
        dataset['last_modified'] = datetime.fromtimestamp(stat_info.st_mtime).isoformat()
        
        # Extract basic stats
        participants_count = len(dataset.get('participants', []))
        analytes_count = len(dataset.get('analytes', {}))
        
        # Count total measurements
        total_measurements = 0
        for participant in dataset.get('participants', []):
            total_measurements += len(participant.get('measurements', []))
        
        dataset['stats'] = {
            'participants_count': participants_count,
            'analytes_count': analytes_count,
            'total_measurements': total_measurements
        }
        
        return dataset
        
    except Exception as e:
        logger.error(f"Error loading dataset {dataset_id} from {file_path}: {e}")
        return None


async def load_all_datasets_async(data_dir: pathlib.Path, max_concurrent: int = 10) -> Dict[str, Dict[str, Any]]:
    """
    Asynchronously load all datasets with controlled concurrency.
    
    Args:
        data_dir: Directory containing dataset files
        max_concurrent: Maximum number of concurrent file operations
        
    Returns:
        Dictionary of loaded datasets keyed by dataset_id
    """
    if not data_dir.exists():
        logger.error(f"Data directory not found: {data_dir}")
        return {}
    
    # Find all YAML files in subdirectories
    yaml_files = list(data_dir.glob("*/*.yaml"))
    yaml_files = [f for f in yaml_files if f.name != ".schema.yaml"]
    
    logger.info(f"Found {len(yaml_files)} YAML files to load")
    
    # Create semaphore to limit concurrent operations
    semaphore = asyncio.Semaphore(max_concurrent)
    
    async def load_with_semaphore(file_path: pathlib.Path) -> tuple:
        """Load a single dataset with semaphore control."""
        async with semaphore:
            dataset_id = file_path.stem
            dataset = await load_dataset_async(file_path, dataset_id)
            return dataset_id, dataset
    
    # Load all datasets concurrently with controlled concurrency
    tasks = [load_with_semaphore(yaml_file) for yaml_file in yaml_files]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    # Process results
    datasets = {}
    successful_loads = 0
    failed_loads = 0
    
    for result in results:
        if isinstance(result, Exception):
            logger.error(f"Dataset load failed with exception: {result}")
            failed_loads += 1
        else:
            dataset_id, dataset = result
            if dataset is not None:
                datasets[dataset_id] = dataset
                successful_loads += 1
            else:
                failed_loads += 1
    
    logger.info(f"Loaded {successful_loads} datasets successfully, {failed_loads} failed")
    return datasets


async def write_file_async(file_path: pathlib.Path, content: str, encoding: str = 'utf-8') -> None:
    """
    Asynchronously write content to a file.
    
    Args:
        file_path: Path to write to
        content: Content to write
        encoding: File encoding (default: utf-8)
    """
    async with aiofiles.open(file_path, 'w', encoding=encoding) as f:
        await f.write(content)


async def export_dataset_to_csv_async(dataset: Dict[str, Any]) -> str:
    """
    Convert a dataset to CSV format asynchronously.
    
    Args:
        dataset: Dataset dictionary
        
    Returns:
        CSV content as string
    """
    import io
    import csv
    
    output = io.StringIO()
    
    # Get all unique measurement fields across participants
    all_fields = set(['participant_id'])
    
    # Add participant attribute fields
    for i, participant in enumerate(dataset.get('participants', [])):
        if 'attributes' in participant:
            for key in participant['attributes'].keys():
                all_fields.add(f'participant_{key}')
    
    # Add measurement fields
    for participant in dataset.get('participants', []):
        for measurement in participant.get('measurements', []):
            for key in measurement.keys():
                all_fields.add(f'measurement_{key}')
    
    # Convert to sorted list for consistent output
    fieldnames = sorted(list(all_fields))
    
    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    
    # Write data rows
    for i, participant in enumerate(dataset.get('participants', [])):
        participant_id = f"participant_{i}"
        
        # Get participant attributes
        base_row = {'participant_id': participant_id}
        if 'attributes' in participant:
            for key, value in participant['attributes'].items():
                base_row[f'participant_{key}'] = value
        
        # Write a row for each measurement
        for measurement in participant.get('measurements', []):
            row = base_row.copy()
            for key, value in measurement.items():
                row[f'measurement_{key}'] = value
            writer.writerow(row)
    
    return output.getvalue()


async def export_dataset_to_json_async(dataset: Dict[str, Any], pretty: bool = True) -> str:
    """
    Convert a dataset to JSON format asynchronously.
    
    Args:
        dataset: Dataset dictionary
        pretty: Whether to format JSON with indentation
        
    Returns:
        JSON content as string
    """
    import json
    
    if pretty:
        return json.dumps(dataset, indent=2, ensure_ascii=False, default=str)
    else:
        return json.dumps(dataset, ensure_ascii=False, default=str)


async def export_dataset_to_excel_async(dataset: Dict[str, Any]) -> bytes:
    """
    Convert a dataset to Excel format asynchronously.
    
    Args:
        dataset: Dataset dictionary
        
    Returns:
        Excel content as bytes
    """
    import pandas as pd
    import io
    
    # Prepare data for Excel export
    rows = []
    
    for i, participant in enumerate(dataset.get('participants', [])):
        participant_id = f"participant_{i}"
        
        # Get participant attributes
        base_data = {'participant_id': participant_id}
        if 'attributes' in participant:
            for key, value in participant['attributes'].items():
                base_data[f'participant_{key}'] = value
        
        # Add measurements
        for measurement in participant.get('measurements', []):
            row_data = base_data.copy()
            for key, value in measurement.items():
                row_data[f'measurement_{key}'] = value
            rows.append(row_data)
    
    # Create DataFrame
    df = pd.DataFrame(rows)
    
    # Create Excel file in memory
    output = io.BytesIO()
    
    # Create Excel writer with multiple sheets
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # Main data sheet
        df.to_excel(writer, sheet_name='Data', index=False)
        
        # Metadata sheet
        metadata = {
            'Field': ['Title', 'Description', 'DOI', 'URL', 'Dataset ID', 'Participants Count', 'Analytes Count', 'Total Measurements'],
            'Value': [
                dataset.get('title', ''),
                dataset.get('description', ''),
                dataset.get('doi', ''),
                dataset.get('url', ''),
                dataset.get('dataset_id', ''),
                dataset.get('stats', {}).get('participants_count', 0),
                dataset.get('stats', {}).get('analytes_count', 0),
                dataset.get('stats', {}).get('total_measurements', 0)
            ]
        }
        metadata_df = pd.DataFrame(metadata)
        metadata_df.to_excel(writer, sheet_name='Metadata', index=False)
        
        # Analytes sheet
        if 'analytes' in dataset:
            analytes_data = []
            for analyte_id, analyte in dataset['analytes'].items():
                analyte_row = {'analyte_id': analyte_id}
                analyte_row.update(analyte)
                analytes_data.append(analyte_row)
            
            if analytes_data:
                analytes_df = pd.DataFrame(analytes_data)
                analytes_df.to_excel(writer, sheet_name='Analytes', index=False)
    
    return output.getvalue()