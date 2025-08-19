import matplotlib

matplotlib.use("Agg")

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
import shedding_hub as sh


# Sample minimal datasets for testing
@pytest.fixture
def minimal_dataset():
    return {
        "dataset_id": "test_dataset",
        "analytes": {
            "A": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "A", "value": 1.0, "time": 0},
                    {"analyte": "A", "value": 2.0, "time": 1},
                    {"analyte": "A", "value": "negative", "time": 2},
                ]
            },
            {
                "measurements": [
                    {"analyte": "A", "value": "negative", "time": 0},
                    {"analyte": "A", "value": 3.0, "time": 1},
                ]
            },
        ],
    }


@pytest.fixture
def minimal_dataset_2():
    return {
        "dataset_id": "test_dataset_2",
        "analytes": {
            "B": {
                "specimen": "nasal",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
            },
            "C": {
                "specimen": "serum",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
            },
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "B", "value": 1.0, "time": 0},
                    {"analyte": "B", "value": "negative", "time": 5},
                    {"analyte": "C", "value": 2.0, "time": 1},
                    {"analyte": "C", "value": 3.0, "time": 3},
                ]
            },
            {
                "measurements": [
                    {"analyte": "B", "value": 2.0, "time": 1},
                    {"analyte": "B", "value": 3.0, "time": 3},
                    {"analyte": "C", "value": "negative", "time": 0},
                    {"analyte": "C", "value": 1.0, "time": 2},
                ]
            },
        ],
    }


def test_calc_shedding_duration_valid(minimal_dataset):
    df = sh.calc_shedding_duration(minimal_dataset)
    assert not df.empty
    assert "shedding_duration_min" in df.columns
    assert df["dataset_id"].iloc[0] == "test_dataset"


def test_calc_shedding_duration_invalid():
    with pytest.raises(ValueError):
        sh.calc_shedding_duration({})
    with pytest.raises(ValueError):
        sh.calc_shedding_duration({"foo": "bar"})


@patch("shedding_hub.shedding_duration.sh.load_dataset")
def test_calc_shedding_durations_valid(mock_load_dataset, minimal_dataset):
    mock_load_dataset.return_value = minimal_dataset
    df = sh.calc_shedding_durations(["test_dataset"])
    assert not df.empty
    assert "shedding_duration_mean" in df.columns


@patch("shedding_hub.shedding_duration.sh.load_dataset")
def test_calc_shedding_durations_all_fail(mock_load_dataset):
    mock_load_dataset.side_effect = Exception("fail")
    with pytest.raises(Exception):
        df = sh.calc_shedding_durations(["bad_dataset"])


def test_plot_shedding_duration(minimal_dataset):
    fig = sh.calc_shedding_duration(minimal_dataset, plotting=True)
    assert fig is not None


def test_plot_shedding_duration_empty():
    with pytest.raises(ValueError):
        sh.plot_shedding_duration(pd.DataFrame())


@patch("shedding_hub.shedding_duration.sh.load_dataset")
def test_plot_shedding_durations_with_ids(
    mock_load_dataset, minimal_dataset, minimal_dataset_2
):
    mock_load_dataset.side_effect = [minimal_dataset, minimal_dataset_2]
    df1 = sh.calc_shedding_duration(minimal_dataset)
    df2 = sh.calc_shedding_duration(minimal_dataset_2)
    df = pd.concat([df1, df2], ignore_index=True)
    fig = sh.plot_shedding_durations(df, biomarker="SARS-CoV-2")
    import matplotlib.figure

    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_durations_empty():
    with pytest.raises(ValueError):
        sh.plot_shedding_durations(pd.DataFrame(), biomarker="SARS-CoV-2")
