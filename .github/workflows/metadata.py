import pandas as pd
import requests
import yaml
import json
import shedding_hub as sh
import os
import glob
from pathlib import Path


def create_summary(schema: dict) -> pd.DataFrame:
    """
    Create data structure of data summary using schema file.

    Args:
        schema: Schema file loaded from GitHub repository.

    Returns:
        Empty data summary dataframe with column names.
    """
    # Initialize column names as a list with ID column
    column_names = ["ID"]
    # Add required fields of analyte into it;
    column_names = column_names + schema["$defs"]["analyte_specification"]["required"]
    # Add summary statistics into it;
    column_names = column_names + [
        "n_samples",
        "n_unique_participants",
        "n_negative",
        "n_positive",
        "n_quantified",
    ]
    return pd.DataFrame(columns=column_names)


def safe_compare(a, b, operator):
    def compare(x, y):
        # Allow int and float to compare with each other
        if not (isinstance(x, (int, float)) and isinstance(y, (int, float))):
            return False

        if operator == "==":
            return x == y
        elif operator == "<":
            return x < y
        elif operator == ">":
            return x > y
        elif operator == "<=":
            return x <= y
        elif operator == ">=":
            return x >= y
        else:
            raise ValueError(f"Unsupported operator: {operator}")

    if isinstance(a, list):
        return [compare(item, b) for item in a]
    else:
        return compare(a, b)


def append_summary(summary: pd.DataFrame, dataset: str) -> pd.DataFrame:
    """
    Extract data from yaml file and append it to the data summary. This function loads the dataset using the shedding-hub module, processes the measurements, and calculates summary statistics for each analyte. The results are appended to the data summary dataframe.

    Args:
        summary: data summary dataframe with column names.
        dataset: Dataset identifier, e.g., :code:`woelfel2020virological`.

    Returns:
        data summary dataframe with extracted data appended.
    """
    data = sh.load_dataset(dataset)

    # create measurements dataframe for the new dataset
    measurements = pd.DataFrame(columns=["participant_ID", "analyte", "time", "value"])
    counter = 1
    for participant in data["participants"]:
        new_data = pd.DataFrame.from_dict(participant["measurements"])
        new_data["participant_ID"] = counter
        if len(measurements) == 0:
            measurements = new_data[["participant_ID", "analyte", "time", "value"]]
        else:
            measurements = pd.concat(
                [measurements, new_data], axis=0, join="inner", ignore_index=True
            )
        counter += 1

    for analyte in data["analytes"]:
        measurements_filtered = measurements.loc[measurements["analyte"] == analyte]
        # calculate summary statistics
        n_samples = len(measurements_filtered)
        n_unique_participants = len(set(measurements_filtered["participant_ID"]))
        n_negative = len(
            measurements_filtered.loc[
                [
                    (value == "negative")
                    | (
                        (data["analytes"][analyte]["unit"] != "cycle threshold")
                        & (
                            safe_compare(
                                value,
                                data["analytes"][analyte]["limit_of_detection"],
                                "<",
                            )
                        )
                    )
                    | (
                        (data["analytes"][analyte]["unit"] == "cycle threshold")
                        & (
                            safe_compare(
                                value,
                                data["analytes"][analyte]["limit_of_detection"],
                                ">",
                            )
                        )
                    )
                    for value in measurements_filtered["value"]
                ]
            ]
        )
        n_positive = len(
            measurements_filtered.loc[
                [
                    (value == "positive")
                    | (
                        (data["analytes"][analyte]["unit"] != "cycle threshold")
                        & (
                            safe_compare(
                                value,
                                data["analytes"][analyte]["limit_of_quantification"],
                                "<=",
                            )
                        )
                    )
                    | (
                        (data["analytes"][analyte]["unit"] == "cycle threshold")
                        & (
                            safe_compare(
                                value,
                                data["analytes"][analyte]["limit_of_detection"],
                                "<=",
                            )
                            | (
                                (
                                    data["analytes"][analyte]["limit_of_detection"]
                                    == "unknown"
                                )
                                & (safe_compare(value, 0, ">"))
                            )
                        )
                    )
                    for value in measurements_filtered["value"]
                ]
            ]
        )
        n_quantified = len(
            measurements_filtered.loc[
                [
                    (
                        (data["analytes"][analyte]["unit"] != "cycle threshold")
                        & (
                            data["analytes"][analyte]["limit_of_quantification"]
                            == "unknown"
                        )
                        & (data["analytes"][analyte]["limit_of_detection"] == "unknown")
                        & (safe_compare(value, 0, ">"))
                    )
                    | (
                        (data["analytes"][analyte]["unit"] != "cycle threshold")
                        & (
                            data["analytes"][analyte]["limit_of_quantification"]
                            == "unknown"
                        )
                        & (data["analytes"][analyte]["limit_of_detection"] != "unknown")
                        & (
                            safe_compare(
                                value,
                                data["analytes"][analyte]["limit_of_detection"],
                                ">=",
                            )
                        )
                    )
                    | (
                        (data["analytes"][analyte]["unit"] != "cycle threshold")
                        & (
                            data["analytes"][analyte]["limit_of_quantification"]
                            != "unknown"
                        )
                        & (
                            safe_compare(
                                value,
                                data["analytes"][analyte]["limit_of_quantification"],
                                ">=",
                            )
                        )
                    )
                    for value in measurements_filtered["value"]
                ]
            ]
        )
        new_line = (
            [dataset]
            + [
                data["analytes"][analyte][key]
                for key in summary.columns
                if key in list(data["analytes"][analyte].keys())
            ]
            + [n_samples, n_unique_participants, n_negative, n_positive, n_quantified]
        )
        summary.loc[len(summary)] = new_line
    return summary


def generate_jsonld(
    dataset: str,
    *,
    repo: str = "shedding-hub/shedding-hub",
    ref: str = "main",
    version: str = "0.1.0",
) -> dict:
    """
    Generate JSON-LD metadata dictionary for a dataset. This function uses the shedding-hub module to load dataset and create a structured JSON-LD representation.

    Parameters:
    - dataset (str): Dataset identifier (used in URLs and @id). e.g., :code:`woelfel2020virological`.
    - repo (str): GitHub repo in format "owner/repo".
    - ref (str): Git reference (branch, tag, or commit) to use for the dataset.
    - version (str): Dataset version string.
    - sh (module): shedding-hub module with get_publication_date(doi) function.

    Returns:
    - dict: JSON-LD metadata.
    """
    data = sh.load_dataset(dataset)
    doi = data.get("doi")

    # Validate required fields
    for field in ["title", "description"]:
        if field not in data or not data[field]:
            raise ValueError(f"Missing required metadata field: {field}")

    # Get publication date with fallback
    pub_date = sh.get_publication_date(doi) if doi else "unknown"

    # Validate version fallback
    if not version:
        version = "unknown"

    json_ld = {
        "@context": "https://schema.org",
        "@type": "Dataset",
        "@id": dataset,
        "name": data["title"],
        "description": data["description"],
        "url": f"https://github.com/{repo}/tree/{ref}/data/{dataset}",
        "identifier": doi,
        "keywords": [
            "pathogen",
            "shedding",
            "biomarker",
            "infectious disease",
            "wastewater surveillance",
        ],
        "creator": [
            {
                "@type": "Person",
                "name": "Yuke Wang",
                "affiliation": {"@type": "Organization", "name": "Example University"},
                "@id": "https://orcid.org/0000-0002-9615-7859",
            },
            {
                "@type": "Person",
                "name": "Till Hoffmann",
                "affiliation": {"@type": "Organization", "name": "Harvard University"},
                "@id": "https://orcid.org/0000-0003-4403-0722",
            },
        ],
        "publisher": {"@type": "Organization", "name": "Shedding Hub Organization"},
        "datePublished": pub_date,
        "license": "https://opensource.org/licenses/MIT",
        "version": version,
        "distribution": {
            "@type": "DataDownload",
            "encodingFormat": "application/yaml",
            "contentUrl": f"https://raw.githubusercontent.com/{repo}/{ref}/data/{dataset}/{dataset}.yaml",
        },
    }

    return json_ld


def save_jsonld(
    json_ld: dict, directory: str, filename: str = "metadata.jsonld"
) -> None:
    """
    Save JSON-LD dictionary to a JSON file in the specified directory.

    Parameters:
    - json_ld (dict): JSON-LD metadata dictionary.
    - directory (str): Path to the directory to save the file.
    - filename (str): Name of the JSON file (default: "metadata.jsonld").

    Returns:
    - str: Full path to the saved JSON file.
    """
    # Create directory if not exists
    os.makedirs(directory, exist_ok=True)

    # Full file path
    file_path = os.path.join(directory, filename)

    # Write JSON-LD to file with indentation for readability
    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(json_ld, f, ensure_ascii=False, indent=2)

    return file_path


def __main__() -> None:
    """
    Main function to generate metadata (JSON-LD) for all datasets in the repository.
    It also fetches the schema from the GitHub repository, processes each dataset, and saves the summary of datasets to a YAML file.
    """
    repo = "shedding-hub/shedding-hub"
    ref = "main"
    version = "0.1.0"
    # create summary of datasets
    print("Create summary for all datasets in the repository...")
    url = f"https://raw.githubusercontent.com/{repo}/{ref}/data/.schema.yaml"
    response = requests.get(url)
    schema = yaml.safe_load(response.text)
    summary = create_summary(schema)

    # append all the data into the summary and create JSON-LD metadata for each dataset
    print("Append data to the summary and create JSON-LD metadata for each dataset...")
    for filename in [Path(file).stem for file in Path("data").glob("*/*.yaml")]:
        print(f"Load the data: {filename}")
        summary = append_summary(summary, filename)
        # create JSON-LD metadata
        metadata = generate_jsonld(
            dataset=filename, repo=repo, ref=ref, version=version
        )
        # save JSON-LD metadata
        print(f"Save JSON-LD metadata for {filename}...")
        save_jsonld(
            json_ld=metadata, directory=f"data/{filename}", filename="metadata.jsonld"
        )

    # remove "\n" from description
    summary["description"] = [item.replace("\n", "") for item in summary["description"]]

    # save summary to YAML file
    with open("data/summary.yaml", "w") as f:
        yaml.safe_dump(summary.to_dict(orient="records"), f)


if __name__ == "__main__":
    __main__()
