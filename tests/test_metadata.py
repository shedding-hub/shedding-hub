import pytest
import pandas as pd
import yaml
import requests
import shedding_hub as sh


def create_meta(schema: dict) -> pd.DataFrame:
    """
    Create data structure of metadata using schema file.

    Args:
        schema: Schema file loaded from GitHub repository.

    Returns:
        Empty Metadata dataframe with column names.
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


@pytest.mark.parametrize(
    "schema, expected",
    [
        (
            yaml.safe_load(
                requests.get(
                    "https://raw.githubusercontent.com/shedding-hub/shedding-hub/refs/heads/main/data/.schema.yaml"
                ).text
            ),
            {
                "length": 0,
                "names": [
                    "ID",
                    "biomarker",
                    "description",
                    "limit_of_detection",
                    "limit_of_quantification",
                    "reference_event",
                    "specimen",
                    "unit",
                    "n_samples",
                    "n_unique_participants",
                    "n_negative",
                    "n_positive",
                    "n_quantified",
                ],
            },
        )
    ],
)
def test_create_meta(schema: dict, expected: dict) -> None:
    temp = create_meta(schema)
    assert len(temp) == expected["length"]
    assert list(temp.columns) == expected["names"]


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


@pytest.mark.parametrize(
    "a, b, operator, expected",
    [
        (10, 5, ">", True),
        (5.5, 10, "<=", True),
        ("5.5", 10, "<", False),
        (35, "unknown", "<", False),
        ("nagtive", 100, "<=", False),
        ("positive", "unknown", "<", False),
    ],
)
def test_safe_compare(
    a: int | float | str, b: int | float | str, operator: str, expected: bool
) -> FileNotFoundError:
    assert safe_compare(a, b, operator) == expected


def append_meta(metadata: pd.DataFrame, dataset: str) -> pd.DataFrame:
    """
    Extract data from yaml file and append it to the metadata.

    Args:
        metadata: Metadata dataframe with column names.
        dataset: Dataset identifier, e.g., :code:`woelfel2020virological`.

    Returns:
        Metadata dataframe with extracted data appended.
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
                for key in metadata.columns
                if key in list(data["analytes"][analyte].keys())
            ]
            + [n_samples, n_unique_participants, n_negative, n_positive, n_quantified]
        )
        metadata.loc[len(metadata)] = new_line
    return metadata


@pytest.mark.parametrize(
    "schema, dataset, expected",
    [
        (
            yaml.safe_load(
                requests.get(
                    "https://raw.githubusercontent.com/shedding-hub/shedding-hub/refs/heads/main/data/.schema.yaml"
                ).text
            ),
            "woelfel2020virological",
            pd.DataFrame(
                {
                    "ID": ["woelfel2020virological"] * 3,
                    "biomarker": ["SARS-CoV-2"] * 3,
                    "description": [
                        'RNA gene copy concentration in stool samples. The authors report that "stool samples were taken and shipped in native conditions," suggesting that results reported as gene copies per gram refer to wet weight.',
                        "RNA gene copy concentration in sputum samples. Results are reported as gene copies per mL.",
                        "Number of gene copies per throat swab.",
                    ],
                    "limit_of_detection": ["unknown"] * 3,
                    "limit_of_quantification": [100] * 3,
                    "reference_event": ["symptom onset"] * 3,
                    "specimen": ["stool", "sputum", "oropharyngeal_swab"],
                    "unit": ["gc/mL", "gc/mL", "gc/swab"],
                    "n_samples": [82, 147, 153],
                    "n_unique_participants": [9, 9, 9],
                    "n_negative": [13, 24, 57],
                    "n_positive": [2, 1, 40],
                    "n_quantified": [67, 122, 56],
                }
            ),
        )
    ],
)
def test_append_meta(schema: dict, dataset: str, expected: pd.DataFrame) -> None:
    temp = create_meta(schema)
    temp = append_meta(temp, dataset)
    for i in range(len(temp)):
        assert all(temp.loc[i] == expected.loc[i])
