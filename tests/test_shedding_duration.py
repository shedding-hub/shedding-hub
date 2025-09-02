import matplotlib

matplotlib.use("Agg")

import pytest
import pandas as pd
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
    assert "first_detect" in df.columns
    assert "last_sample" in df.columns
    assert df["dataset_id"].iloc[0] == "test_dataset"


def test_calc_shedding_duration_invalid():
    with pytest.raises(ValueError):
        sh.calc_shedding_duration({})
    with pytest.raises(ValueError):
        sh.calc_shedding_duration({"foo": "bar"})


def test_calc_shedding_durations_valid():
    df = sh.calc_shedding_durations(["woelfel2020virological"])
    assert not df.empty

    # Check first row values
    expected_first_row = pd.Series(
        {
            "dataset_id": "woelfel2020virological",
            "biomarker": "SARS-CoV-2",
            "specimen": "oropharyngeal_swab",
            "reference_event": "symptom onset",
            "shedding_duration_min": 3,
            "shedding_duration_max": 23,
            "shedding_duration_mean": 13.111111,
            "n_sample": 153,
            "n_participant": 9,
        }
    )

    # Get first row as Series and check each value
    first_row = df.iloc[0].copy()
    first_row.name = None  # Reset the index name
    pd.testing.assert_series_equal(
        first_row[expected_first_row.index], expected_first_row, check_names=True
    )


def test_calc_shedding_durations_all_fail():
    with pytest.raises(Exception):
        sh.calc_shedding_durations(["invalid_dataset_id"])


def test_plot_shedding_duration(minimal_dataset):
    # First get the individual level data
    df_individual = sh.calc_shedding_duration(minimal_dataset, output="individual")
    # Then plot it
    fig = sh.plot_shedding_duration(df_individual)
    assert fig is not None


def test_plot_shedding_duration_empty():
    with pytest.raises(ValueError):
        sh.plot_shedding_duration(pd.DataFrame())


def test_plot_shedding_durations_with_ids(minimal_dataset, minimal_dataset_2):
    # Get summary data for both datasets
    df1 = sh.calc_shedding_duration(minimal_dataset, output="summary")
    df2 = sh.calc_shedding_duration(minimal_dataset_2, output="summary")
    df = pd.concat([df1, df2], ignore_index=True)
    fig = sh.plot_shedding_durations(df, biomarker="SARS-CoV-2")
    import matplotlib.figure

    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_durations_empty():
    with pytest.raises(ValueError):
        sh.plot_shedding_durations(pd.DataFrame(), biomarker="SARS-CoV-2")
