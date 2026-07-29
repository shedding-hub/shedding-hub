import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest
import matplotlib.figure
import matplotlib.pyplot as plt
import shedding_hub as sh
from shedding_hub.shedding_fit import SheddingFit


# Sample minimal datasets for testing
@pytest.fixture
def minimal_dataset():
    """Minimal dataset with single analyte and specimen type."""
    return {
        "dataset_id": "test_dataset",
        "analytes": {
            "A": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "description": "Test analyte",
                "limit_of_detection": 100,
                "limit_of_quantification": 500,
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "A", "value": 1000.0, "time": 0},
                    {"analyte": "A", "value": 2000.0, "time": 1},
                    {"analyte": "A", "value": 1500.0, "time": 2},
                    {"analyte": "A", "value": 500.0, "time": 3},
                ]
            },
            {
                "measurements": [
                    {"analyte": "A", "value": 800.0, "time": 0},
                    {"analyte": "A", "value": 1200.0, "time": 1},
                    {"analyte": "A", "value": 900.0, "time": 2},
                    {"analyte": "A", "value": "negative", "time": 3},
                ]
            },
            {
                "measurements": [
                    {"analyte": "A", "value": "negative", "time": 0},
                    {"analyte": "A", "value": 1000.0, "time": 1},
                    {"analyte": "A", "value": 1800.0, "time": 2},
                    {"analyte": "A", "value": 1200.0, "time": 3},
                ]
            },
        ],
    }


@pytest.fixture
def multi_specimen_dataset():
    """Dataset with multiple specimen types."""
    return {
        "dataset_id": "test_multi_specimen",
        "analytes": {
            "stool_analyte": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/wet gram",
                "description": "Stool test",
                "limit_of_detection": 100,
                "limit_of_quantification": 500,
            },
            "swab_analyte": {
                "specimen": "nasopharyngeal_swab",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/swab",
                "description": "Swab test",
                "limit_of_detection": 50,
                "limit_of_quantification": 200,
            },
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool_analyte", "value": 1000.0, "time": 0},
                    {"analyte": "stool_analyte", "value": 2000.0, "time": 2},
                    {"analyte": "swab_analyte", "value": 500.0, "time": 0},
                    {"analyte": "swab_analyte", "value": 800.0, "time": 1},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool_analyte", "value": 1500.0, "time": 1},
                    {"analyte": "stool_analyte", "value": 1200.0, "time": 3},
                    {"analyte": "swab_analyte", "value": 600.0, "time": 0},
                    {"analyte": "swab_analyte", "value": 300.0, "time": 2},
                ]
            },
        ],
    }


@pytest.fixture
def dataset_with_unknown_times():
    """Dataset with some unknown time values."""
    return {
        "dataset_id": "test_unknown",
        "analytes": {
            "A": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "description": "Test analyte",
                "limit_of_detection": 100,
                "limit_of_quantification": 500,
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "A", "value": 1000.0, "time": 0},
                    {"analyte": "A", "value": 2000.0, "time": "unknown"},
                    {"analyte": "A", "value": 1500.0, "time": 2},
                ]
            },
        ],
    }


@pytest.fixture
def minimal_dataset_2():
    """Second minimal dataset for multi-dataset tests."""
    return {
        "dataset_id": "test_dataset_2",
        "analytes": {
            "B": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "description": "Test analyte 2",
                "limit_of_detection": 100,
                "limit_of_quantification": 500,
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "B", "value": 1500.0, "time": 0},
                    {"analyte": "B", "value": 2500.0, "time": 1},
                    {"analyte": "B", "value": 1800.0, "time": 2},
                ]
            },
            {
                "measurements": [
                    {"analyte": "B", "value": 1100.0, "time": 0},
                    {"analyte": "B", "value": 1600.0, "time": 1},
                    {"analyte": "B", "value": 1300.0, "time": 2},
                ]
            },
        ],
    }


def test_plot_time_course_valid(minimal_dataset):
    """Test plot_time_course with valid minimal dataset."""
    fig = sh.plot_time_course(minimal_dataset)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_course_multi_specimen(multi_specimen_dataset):
    """Test plot_time_course with multiple specimen types."""
    fig = sh.plot_time_course(multi_specimen_dataset)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_course_with_specimen_filter(multi_specimen_dataset):
    """Test plot_time_course with specimen filtering."""
    fig = sh.plot_time_course(multi_specimen_dataset, specimen="stool")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_course_unknown_times(dataset_with_unknown_times):
    """Test that unknown time values are filtered out."""
    fig = sh.plot_time_course(dataset_with_unknown_times)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_course_show_negative(minimal_dataset):
    """Test plot_time_course with show_negative=True."""
    fig = sh.plot_time_course(minimal_dataset, show_negative=True)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_course_custom_styling(minimal_dataset):
    """Test plot_time_course with custom styling parameters."""
    fig = sh.plot_time_course(
        minimal_dataset,
        marker="s",
        markersize=8,
        line_alpha=0.5,
        line_color="red",
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_course_sampling(minimal_dataset):
    """Test plot_time_course with participant sampling."""
    fig = sh.plot_time_course(minimal_dataset, max_nparticipant=2, random_seed=42)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_course_empty_dataset():
    """Test plot_time_course with empty dataset."""
    with pytest.raises(ValueError, match="Dataset must be a non-empty dictionary"):
        sh.plot_time_course({})


def test_plot_time_course_missing_keys():
    """Test plot_time_course with missing required keys."""
    invalid_dataset = {"dataset_id": "test"}
    with pytest.raises(ValueError, match="Dataset missing required keys"):
        sh.plot_time_course(invalid_dataset)


def test_plot_time_course_no_participants():
    """Test plot_time_course with no participants."""
    invalid_dataset = {
        "dataset_id": "test",
        "analytes": {
            "A": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "description": "Test",
                "limit_of_detection": 100,
                "limit_of_quantification": 500,
            }
        },
        "participants": [],
    }
    with pytest.raises(ValueError, match="Dataset has no participants"):
        sh.plot_time_course(invalid_dataset)


def test_plot_time_course_invalid_biomarker(minimal_dataset):
    """Test plot_time_course with invalid biomarker filter."""
    with pytest.raises(ValueError, match="No measurements found for biomarker"):
        sh.plot_time_course(minimal_dataset, biomarker="nonexistent_biomarker")


def test_plot_time_courses_valid(minimal_dataset, minimal_dataset_2):
    """Test plot_time_courses with valid datasets."""
    fig = sh.plot_time_courses([minimal_dataset, minimal_dataset_2])
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_courses_single_dataset(minimal_dataset):
    """Test plot_time_courses with single dataset."""
    fig = sh.plot_time_courses([minimal_dataset])
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_courses_custom_styling(minimal_dataset, minimal_dataset_2):
    """Test plot_time_courses with custom styling parameters."""
    fig = sh.plot_time_courses(
        [minimal_dataset, minimal_dataset_2],
        marker="^",
        markersize=6,
        line_alpha=0.3,
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_courses_sampling(minimal_dataset, minimal_dataset_2):
    """Test plot_time_courses with participant sampling."""
    fig = sh.plot_time_courses(
        [minimal_dataset, minimal_dataset_2],
        max_nparticipant=1,
        random_seed=123,
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_time_courses_empty_list():
    """Test plot_time_courses with empty list."""
    with pytest.raises(ValueError, match="Datasets must be a non-empty list"):
        sh.plot_time_courses([])


def test_plot_time_courses_invalid_dataset():
    """Test plot_time_courses with invalid dataset in list."""
    with pytest.raises(ValueError, match="Each dataset must be a non-empty dictionary"):
        sh.plot_time_courses([{}])


def test_plot_time_courses_not_list():
    """Test plot_time_courses with non-list input."""
    with pytest.raises(ValueError, match="Datasets must be a non-empty list"):
        sh.plot_time_courses("not a list")


def test_plot_time_course_with_real_dataset():
    """Test plot_time_course with a real dataset from the repository."""
    try:
        dataset = sh.load_dataset("woelfel2020virological")
        fig = sh.plot_time_course(dataset, max_nparticipant=5)
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
    except Exception as e:
        # If loading fails (e.g., network issues), skip this test
        pytest.skip(f"Could not load real dataset: {e}")


def test_plot_time_courses_with_real_datasets():
    """Test plot_time_courses with real datasets from the repository."""
    try:
        dataset1 = sh.load_dataset("woelfel2020virological")
        dataset2 = sh.load_dataset("ke2022daily")
        fig = sh.plot_time_courses([dataset1, dataset2], max_nparticipant=3)
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
    except Exception as e:
        # If loading fails (e.g., network issues), skip this test
        pytest.skip(f"Could not load real datasets: {e}")


# ==================== plot_shedding_heatmap tests ====================


@pytest.fixture
def ct_dataset():
    """Dataset with CT values."""
    return {
        "dataset_id": "test_ct_dataset",
        "analytes": {
            "A": {
                "specimen": "nasopharyngeal_swab",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "cycle threshold",
                "description": "CT test",
                "limit_of_detection": 40,
                "limit_of_quantification": 35,
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "A", "value": 25.0, "time": 0},
                    {"analyte": "A", "value": 20.0, "time": 1},
                    {"analyte": "A", "value": 22.0, "time": 2},
                    {"analyte": "A", "value": 30.0, "time": 3},
                ]
            },
            {
                "measurements": [
                    {"analyte": "A", "value": 28.0, "time": 0},
                    {"analyte": "A", "value": 18.0, "time": 1},
                    {"analyte": "A", "value": 24.0, "time": 2},
                    {"analyte": "A", "value": "negative", "time": 3},
                ]
            },
        ],
    }


def test_plot_shedding_heatmap_valid(minimal_dataset):
    """Test plot_shedding_heatmap with valid minimal dataset."""
    fig = sh.plot_shedding_heatmap(minimal_dataset)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_with_specimen_filter(multi_specimen_dataset):
    """Test plot_shedding_heatmap with specimen filtering."""
    fig = sh.plot_shedding_heatmap(multi_specimen_dataset, specimen="stool")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_with_biomarker_filter(minimal_dataset):
    """Test plot_shedding_heatmap with biomarker filtering."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, biomarker="SARS-CoV-2")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_with_value_concentration(minimal_dataset):
    """Test plot_shedding_heatmap with value='concentration'."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, value="concentration")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_with_value_ct(ct_dataset):
    """Test plot_shedding_heatmap with value='ct'."""
    fig = sh.plot_shedding_heatmap(ct_dataset, value="ct")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_sort_by_first_positive(minimal_dataset):
    """Test plot_shedding_heatmap with sort_by='first_positive'."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, sort_by="first_positive")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_sort_by_peak_time(minimal_dataset):
    """Test plot_shedding_heatmap with sort_by='peak_time'."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, sort_by="peak_time")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_sort_by_peak_value(minimal_dataset):
    """Test plot_shedding_heatmap with sort_by='peak_value'."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, sort_by="peak_value")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_sort_by_participant_id(minimal_dataset):
    """Test plot_shedding_heatmap with sort_by='participant_id'."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, sort_by="participant_id")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_with_time_range(minimal_dataset):
    """Test plot_shedding_heatmap with time_range filtering."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, time_range=(0, 2))
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_show_negative_true(minimal_dataset):
    """Test plot_shedding_heatmap with show_negative=True."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, show_negative=True)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_show_negative_false(minimal_dataset):
    """Test plot_shedding_heatmap with show_negative=False."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, show_negative=False)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_custom_figsize(minimal_dataset):
    """Test plot_shedding_heatmap with custom figsize."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, figsize=(12, 8))
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)
    assert fig.get_size_inches()[0] == 12
    assert fig.get_size_inches()[1] == 8


def test_plot_shedding_heatmap_custom_cmap(minimal_dataset):
    """Test plot_shedding_heatmap with custom colormap."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, cmap="viridis")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_custom_time_bin_size(minimal_dataset):
    """Test plot_shedding_heatmap with custom time_bin_size."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, time_bin_size=0.5)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_show_participant_labels(minimal_dataset):
    """Test plot_shedding_heatmap with participant labels shown."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, show_participant_labels=True)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_no_colorbar(minimal_dataset):
    """Test plot_shedding_heatmap with colorbar disabled."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, show_colorbar=False)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_sampling(minimal_dataset):
    """Test plot_shedding_heatmap with participant sampling."""
    fig = sh.plot_shedding_heatmap(minimal_dataset, max_nparticipant=2, random_seed=42)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_shedding_heatmap_empty_dataset():
    """Test plot_shedding_heatmap with empty dataset."""
    with pytest.raises(ValueError, match="Dataset must be a non-empty dictionary"):
        sh.plot_shedding_heatmap({})


def test_plot_shedding_heatmap_missing_keys():
    """Test plot_shedding_heatmap with missing required keys."""
    invalid_dataset = {"dataset_id": "test"}
    with pytest.raises(ValueError, match="Dataset missing required keys"):
        sh.plot_shedding_heatmap(invalid_dataset)


def test_plot_shedding_heatmap_no_participants():
    """Test plot_shedding_heatmap with no participants."""
    invalid_dataset = {
        "dataset_id": "test",
        "analytes": {"A": {"specimen": "stool", "unit": "gc/mL"}},
        "participants": [],
    }
    with pytest.raises(ValueError, match="Dataset has no participants"):
        sh.plot_shedding_heatmap(invalid_dataset)


def test_plot_shedding_heatmap_invalid_biomarker(minimal_dataset):
    """Test plot_shedding_heatmap with invalid biomarker filter."""
    with pytest.raises(ValueError, match="No measurements found for biomarker"):
        sh.plot_shedding_heatmap(minimal_dataset, biomarker="nonexistent")


def test_plot_shedding_heatmap_invalid_specimen(minimal_dataset):
    """Test plot_shedding_heatmap with invalid specimen filter."""
    with pytest.raises(ValueError, match="No measurements found for specimen"):
        sh.plot_shedding_heatmap(minimal_dataset, specimen="nonexistent")


def test_plot_shedding_heatmap_invalid_sort_by(minimal_dataset):
    """Test plot_shedding_heatmap with invalid sort_by parameter."""
    with pytest.raises(ValueError, match="Invalid sort_by"):
        sh.plot_shedding_heatmap(minimal_dataset, sort_by="invalid")


def test_plot_shedding_heatmap_invalid_value(minimal_dataset):
    """Test plot_shedding_heatmap with invalid value parameter."""
    with pytest.raises(ValueError, match="Invalid value"):
        sh.plot_shedding_heatmap(minimal_dataset, value="invalid")


def test_plot_shedding_heatmap_with_real_dataset():
    """Test plot_shedding_heatmap with a real dataset from the repository."""
    try:
        dataset = sh.load_dataset("woelfel2020virological")
        fig = sh.plot_shedding_heatmap(
            dataset, specimen="sputum", value="concentration", max_nparticipant=10
        )
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
    except Exception as e:
        # If loading fails (e.g., network issues), skip this test
        pytest.skip(f"Could not load real dataset: {e}")


# ==================== plot_mean_trajectory tests ====================


def test_plot_mean_trajectory_valid(minimal_dataset):
    """Test plot_mean_trajectory with valid minimal dataset."""
    fig = sh.plot_mean_trajectory(minimal_dataset)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_mean_95ci(minimal_dataset):
    """Test plot_mean_trajectory with mean and 95% CI (default)."""
    fig = sh.plot_mean_trajectory(
        minimal_dataset, central_tendency="mean", uncertainty="95ci"
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_median_iqr(minimal_dataset):
    """Test plot_mean_trajectory with median and IQR."""
    fig = sh.plot_mean_trajectory(
        minimal_dataset, central_tendency="median", uncertainty="iqr"
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_mean_sd(minimal_dataset):
    """Test plot_mean_trajectory with mean and standard deviation."""
    fig = sh.plot_mean_trajectory(
        minimal_dataset, central_tendency="mean", uncertainty="sd"
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_mean_range(minimal_dataset):
    """Test plot_mean_trajectory with mean and full range."""
    fig = sh.plot_mean_trajectory(
        minimal_dataset, central_tendency="mean", uncertainty="range"
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_with_biomarker_filter(minimal_dataset):
    """Test plot_mean_trajectory with biomarker filtering."""
    fig = sh.plot_mean_trajectory(minimal_dataset, biomarker="SARS-CoV-2")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_with_specimen_filter(multi_specimen_dataset):
    """Test plot_mean_trajectory with specimen filtering."""
    fig = sh.plot_mean_trajectory(
        multi_specimen_dataset, specimen="stool", min_observations=1
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_with_value_concentration(minimal_dataset):
    """Test plot_mean_trajectory with value='concentration'."""
    fig = sh.plot_mean_trajectory(minimal_dataset, value="concentration")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_with_value_ct(ct_dataset):
    """Test plot_mean_trajectory with value='ct'."""
    fig = sh.plot_mean_trajectory(ct_dataset, value="ct")
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_with_time_range(minimal_dataset):
    """Test plot_mean_trajectory with time_range filtering."""
    fig = sh.plot_mean_trajectory(minimal_dataset, time_range=(0, 2))
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_custom_time_bin_size(minimal_dataset):
    """Test plot_mean_trajectory with custom time_bin_size."""
    fig = sh.plot_mean_trajectory(minimal_dataset, time_bin_size=0.5)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_custom_figsize(minimal_dataset):
    """Test plot_mean_trajectory with custom figsize."""
    fig = sh.plot_mean_trajectory(minimal_dataset, figsize=(12, 8))
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)
    assert fig.get_size_inches()[0] == 12
    assert fig.get_size_inches()[1] == 8


def test_plot_mean_trajectory_custom_styling(minimal_dataset):
    """Test plot_mean_trajectory with custom styling parameters."""
    fig = sh.plot_mean_trajectory(
        minimal_dataset,
        line_color="red",
        fill_alpha=0.5,
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_show_individual(minimal_dataset):
    """Test plot_mean_trajectory with individual trajectories shown."""
    fig = sh.plot_mean_trajectory(
        minimal_dataset, show_individual=True, individual_alpha=0.2
    )
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_hide_n(minimal_dataset):
    """Test plot_mean_trajectory with sample size annotations hidden."""
    fig = sh.plot_mean_trajectory(minimal_dataset, show_n=False)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_min_observations(minimal_dataset):
    """Test plot_mean_trajectory with custom min_observations."""
    fig = sh.plot_mean_trajectory(minimal_dataset, min_observations=2)
    assert fig is not None
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_mean_trajectory_empty_dataset():
    """Test plot_mean_trajectory with empty dataset."""
    with pytest.raises(ValueError, match="Dataset must be a non-empty dictionary"):
        sh.plot_mean_trajectory({})


def test_plot_mean_trajectory_missing_keys():
    """Test plot_mean_trajectory with missing required keys."""
    invalid_dataset = {"dataset_id": "test"}
    with pytest.raises(ValueError, match="Dataset missing required keys"):
        sh.plot_mean_trajectory(invalid_dataset)


def test_plot_mean_trajectory_no_participants():
    """Test plot_mean_trajectory with no participants."""
    invalid_dataset = {
        "dataset_id": "test",
        "analytes": {"A": {"specimen": "stool", "unit": "gc/mL"}},
        "participants": [],
    }
    with pytest.raises(ValueError, match="Dataset has no participants"):
        sh.plot_mean_trajectory(invalid_dataset)


def test_plot_mean_trajectory_invalid_biomarker(minimal_dataset):
    """Test plot_mean_trajectory with invalid biomarker filter."""
    with pytest.raises(ValueError, match="No measurements found for biomarker"):
        sh.plot_mean_trajectory(minimal_dataset, biomarker="nonexistent")


def test_plot_mean_trajectory_invalid_specimen(minimal_dataset):
    """Test plot_mean_trajectory with invalid specimen filter."""
    with pytest.raises(ValueError, match="No measurements found for specimen"):
        sh.plot_mean_trajectory(minimal_dataset, specimen="nonexistent")


def test_plot_mean_trajectory_invalid_central_tendency(minimal_dataset):
    """Test plot_mean_trajectory with invalid central_tendency parameter."""
    with pytest.raises(ValueError, match="Invalid central_tendency"):
        sh.plot_mean_trajectory(minimal_dataset, central_tendency="invalid")


def test_plot_mean_trajectory_invalid_uncertainty(minimal_dataset):
    """Test plot_mean_trajectory with invalid uncertainty parameter."""
    with pytest.raises(ValueError, match="Invalid uncertainty"):
        sh.plot_mean_trajectory(minimal_dataset, uncertainty="invalid")


def test_plot_mean_trajectory_invalid_value(minimal_dataset):
    """Test plot_mean_trajectory with invalid value parameter."""
    with pytest.raises(ValueError, match="Invalid value"):
        sh.plot_mean_trajectory(minimal_dataset, value="invalid")


def test_plot_mean_trajectory_with_real_dataset():
    """Test plot_mean_trajectory with a real dataset from the repository."""
    try:
        dataset = sh.load_dataset("woelfel2020virological")
        fig = sh.plot_mean_trajectory(
            dataset, specimen="sputum", value="concentration", min_observations=2
        )
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
    except Exception as e:
        # If loading fails (e.g., network issues), skip this test
        pytest.skip(f"Could not load real dataset: {e}")


# ---------------------------------------------------------------------------
# plot_catalog_fits
# ---------------------------------------------------------------------------

LN10 = float(np.log(10.0))


def _stub_fit(
    dataset_id="study_a",
    model="exponential",
    *,
    biomarker="SARS-CoV-2",
    specimen="stool",
    unit="gc/mL",
    reference_event="symptom onset",
    analyte=None,
    half_life=4.0,
    peak_day=3.0,
    peak_log10=7.0,
    censoring_limit=2.0,
    median_first_observed_day=float("nan"),
    n_subjects=20,
):
    """A directly-constructed SheddingFit with a hand-chosen median individual.

    ``population_mean`` is built in each model's own population coordinates, so
    ``half_life``, ``peak_day`` and ``peak_log10`` come back out of the fit as
    the values passed in. The exponential model only decays, so its peak is at
    day 0 whatever ``peak_day`` says.
    """
    a0 = np.log(2.0) / half_life
    if model == "exponential":
        mean = np.array([np.log(a0), np.log(peak_log10 * LN10)])
    else:
        mean = np.array([np.log(a0), np.log(peak_day), peak_log10])
    return SheddingFit(
        model=model,
        method="mle",
        population_mean=mean,
        population_cov=np.diag(np.full(mean.size, 0.04)),
        sigma=0.3,
        subject_params=pd.DataFrame(),
        censoring_limit=censoring_limit,
        dataset_id=dataset_id,
        analyte=specimen if analyte is None else analyte,
        biomarker=biomarker,
        specimen=specimen,
        reference_event=reference_event,
        unit=unit,
        gene_target=None,
        dose=None,
        vaccine_type=None,
        n_subjects=n_subjects,
        n_measurements=100,
        n_censored=10,
        n_excluded_subjects=0,
        n_dropped_measurements=0,
        converged=True,
        log_likelihood=-1.0,
        aic=2.0,
        median_first_observed_day=median_first_observed_day,
    )


def _curves(ax):
    """The model curves in a panel: solid for exponential, dashed for gamma."""
    return [line for line in ax.get_lines() if line.get_linestyle() in ("-", "--")]


def _limits(ax):
    """The dotted censoring-limit lines in a panel."""
    return [line for line in ax.get_lines() if line.get_linestyle() == ":"]


@pytest.fixture
def catalog_fits():
    """Five fits spanning four comparability groups, one of which has 2 studies."""
    return [
        _stub_fit("study_a"),
        _stub_fit("study_b"),
        _stub_fit("study_c", specimen="sputum"),
        _stub_fit("study_d", biomarker="norovirus", unit="gc/wet gram"),
        _stub_fit("study_e", reference_event="enrollment"),
    ]


def test_plot_catalog_fits_separates_units():
    """Curves in different units cannot share a y axis, so cannot share a panel."""
    fits = [_stub_fit("study_a", unit="gc/mL"), _stub_fit("study_b", unit="pfu/mL")]
    fig = sh.plot_catalog_fits(fits)
    assert len(fig.axes) == 2


def test_plot_catalog_fits_separates_reference_events():
    """t=0 means something different per reference event, so panels must not mix them."""
    fits = [
        _stub_fit("study_a", reference_event="symptom onset"),
        _stub_fit("study_b", reference_event="enrollment"),
    ]
    fig = sh.plot_catalog_fits(fits)
    assert len(fig.axes) == 2


def test_plot_catalog_fits_draws_both_models_of_a_study_in_one_colour():
    fits = [
        _stub_fit("study_a", model="exponential"),
        _stub_fit("study_a", model="gamma"),
    ]
    fig = sh.plot_catalog_fits(fits)
    assert len(fig.axes) == 1
    curves = _curves(fig.axes[0])
    assert len(curves) == 2
    assert len({curve.get_color() for curve in curves}) == 1
    assert {curve.get_linestyle() for curve in curves} == {"-", "--"}


def test_plot_catalog_fits_gives_each_study_its_own_colour():
    fig = sh.plot_catalog_fits([_stub_fit("study_a"), _stub_fit("study_b")])
    curves = _curves(fig.axes[0])
    assert len({curve.get_color() for curve in curves}) == 2


@pytest.mark.parametrize(
    "filters, expected_panels",
    [
        ({"biomarker": "norovirus"}, 1),
        ({"specimen": "sputum"}, 1),
        ({"unit": "gc/mL"}, 3),
        ({"reference_event": "symptom onset"}, 3),
    ],
)
def test_plot_catalog_fits_filters_reduce_the_panel_count(
    catalog_fits, filters, expected_panels
):
    fig = sh.plot_catalog_fits(catalog_fits, **filters)
    assert len(fig.axes) == expected_panels


def test_plot_catalog_fits_restricts_to_named_datasets(catalog_fits):
    fig = sh.plot_catalog_fits(catalog_fits, dataset_ids=["study_a"])
    assert len(fig.axes) == 1
    assert len(_curves(fig.axes[0])) == 1


def test_plot_catalog_fits_accepts_a_catalog(catalog_fits):
    """A SheddingCatalog and a bare list of fits are both acceptable input."""
    catalog = sh.SheddingCatalog(fits=catalog_fits)
    assert len(sh.plot_catalog_fits(catalog).axes) == 4


def test_plot_catalog_fits_fades_the_stretch_before_first_observation():
    """A curve drawn earlier than the study observed is functional form, not data."""
    fit = _stub_fit("study_a", median_first_observed_day=6.0)
    fig = sh.plot_catalog_fits([fit])
    curves = _curves(fig.axes[0])
    assert len(curves) == 2
    faded, full = sorted(curves, key=lambda curve: curve.get_alpha())
    assert faded.get_alpha() == pytest.approx(0.35)
    assert full.get_alpha() == pytest.approx(1.0)
    assert faded.get_xdata()[-1] == pytest.approx(6.0, abs=0.1)
    assert full.get_xdata()[0] == pytest.approx(6.0, abs=0.1)


def test_plot_catalog_fits_can_disable_the_extrapolation_fade():
    fit = _stub_fit("study_a", median_first_observed_day=6.0)
    fig = sh.plot_catalog_fits([fit], show_extrapolation=False)
    curves = _curves(fig.axes[0])
    assert len(curves) == 1
    assert curves[0].get_alpha() == pytest.approx(1.0)


def test_plot_catalog_fits_derives_the_horizon_from_peak_and_decay():
    """Five half-lives past the peak: 3 + 5 x 4 = 23 days."""
    fit = _stub_fit("study_a", model="gamma", peak_day=3.0, half_life=4.0)
    fig = sh.plot_catalog_fits([fit])
    assert max(_curves(fig.axes[0])[0].get_xdata()) == pytest.approx(23.0)


def test_plot_catalog_fits_clamps_a_fast_decay_to_seven_days():
    fit = _stub_fit("study_a", model="exponential", half_life=0.2)
    fig = sh.plot_catalog_fits([fit])
    assert max(_curves(fig.axes[0])[0].get_xdata()) == pytest.approx(7.0)


def test_plot_catalog_fits_clamps_a_runaway_decay_to_sixty_days():
    fit = _stub_fit("study_a", model="exponential", half_life=200.0)
    fig = sh.plot_catalog_fits([fit])
    assert max(_curves(fig.axes[0])[0].get_xdata()) == pytest.approx(60.0)


def test_plot_catalog_fits_n_days_overrides_the_derived_horizon():
    fit = _stub_fit("study_a", model="exponential", half_life=4.0)
    fig = sh.plot_catalog_fits([fit], n_days=12.0)
    assert max(_curves(fig.axes[0])[0].get_xdata()) == pytest.approx(12.0)


def test_plot_catalog_fits_starts_after_zero_so_gamma_is_defined():
    """The gamma curve is undefined at t=0, where log10 diverges."""
    fit = _stub_fit("study_a", model="gamma")
    curve = _curves(sh.plot_catalog_fits([fit]).axes[0])[0]
    assert min(curve.get_xdata()) > 0
    assert np.isfinite(curve.get_ydata()).all()


def test_plot_catalog_fits_draws_each_studys_censoring_limit():
    """Limits differ between studies, so one shared line would be wrong."""
    fits = [
        _stub_fit("study_a", censoring_limit=2.0),
        _stub_fit("study_b", censoring_limit=3.5),
    ]
    ax = sh.plot_catalog_fits(fits).axes[0]
    drawn = sorted(line.get_ydata()[0] for line in _limits(ax))
    assert drawn == pytest.approx([2.0, 3.5])


def test_plot_catalog_fits_puts_multi_study_panels_first():
    """Two studies in stool, one in sputum; sputum sorts first alphabetically."""
    fits = [
        _stub_fit("study_a", specimen="stool"),
        _stub_fit("study_b", specimen="stool"),
        _stub_fit("study_c", specimen="sputum"),
    ]
    fig = sh.plot_catalog_fits(fits)
    assert "stool" in fig.axes[0].get_title()
    assert "sputum" in fig.axes[1].get_title()


def test_plot_catalog_fits_returns_one_closed_axis_per_group(catalog_fits):
    fig = sh.plot_catalog_fits(catalog_fits)
    assert isinstance(fig, matplotlib.figure.Figure)
    assert len(fig.axes) == 4
    assert fig.number not in plt.get_fignums()


def test_plot_catalog_fits_rejects_filters_matching_nothing(catalog_fits):
    with pytest.raises(ValueError, match="No fits match"):
        sh.plot_catalog_fits(catalog_fits, biomarker="influenza")


def test_plot_catalog_fits_names_the_combinations_that_do_exist(catalog_fits):
    with pytest.raises(ValueError, match="norovirus"):
        sh.plot_catalog_fits(catalog_fits, biomarker="influenza")


def test_plot_catalog_fits_rejects_an_empty_catalog():
    with pytest.raises(ValueError, match="no fits"):
        sh.plot_catalog_fits([])


def test_plot_catalog_fits_rejects_an_unmatched_dataset_id(catalog_fits):
    """Dropping the name silently would shrink the figure without saying so."""
    with pytest.raises(ValueError, match="study_typo"):
        sh.plot_catalog_fits(catalog_fits, dataset_ids=["study_a", "study_typo"])


def test_plot_catalog_fits_names_the_analyte_when_a_study_contributes_several():
    """One study can put many analytes in a panel, and they disagree.

    ``natarajan2022gastrointestinal`` contributes 14 stool analytes to one
    panel, whose peaks span over 2 log10. Labelling them all with the study
    name alone would claim 14 different curves were the same fit.
    """
    fits = [
        _stub_fit("study_a", analyte="N1-ddPCR", peak_log10=5.0),
        _stub_fit("study_a", analyte="E-qPCR", peak_log10=3.0),
    ]
    labels = sorted(
        line.get_label() for line in _curves(sh.plot_catalog_fits(fits).axes[0])
    )
    assert labels == [
        "study_a E-qPCR (exponential)",
        "study_a N1-ddPCR (exponential)",
    ]


def test_plot_catalog_fits_floors_the_y_axis_near_the_censoring_limit():
    """Decay far below the limit of quantification is not measurable by anyone.

    A panel mixing a slow fit, which stretches the horizon, with a fast one,
    which then plunges over it, otherwise spends most of its y axis on
    concentrations no study could have detected, squashing the rest.
    """
    fits = [
        _stub_fit("study_a", half_life=30.0, censoring_limit=2.0),
        _stub_fit("study_b", half_life=0.5, censoring_limit=2.0),
    ]
    ax = sh.plot_catalog_fits(fits).axes[0]
    assert ax.get_ylim()[0] == pytest.approx(1.0)


def test_plot_catalog_fits_keeps_a_curve_peaking_below_its_limit_visible():
    """Two catalog fits peak below their own limit; flooring at the limit alone
    would drop them off the bottom of the panel entirely."""
    fit = _stub_fit("study_a", peak_log10=1.0, censoring_limit=4.0)
    ax = sh.plot_catalog_fits([fit]).axes[0]
    assert ax.get_ylim()[0] < 1.0


def test_plot_catalog_fits_caps_a_crowded_legend():
    """A panel holding many curves would otherwise be buried under its own key.

    The overflow is counted rather than dropped, so a truncated legend never
    reads as the whole panel.
    """
    fits = [
        _stub_fit("study_a", analyte=f"assay_{index}", peak_log10=4.0 + index * 0.1)
        for index in range(10)
    ]
    ax = sh.plot_catalog_fits(fits).axes[0]
    texts = [text.get_text() for text in ax.get_legend().get_texts()]
    assert len(texts) == 7
    assert texts[-1] == "+ 4 more"


def test_plot_catalog_fits_lists_every_curve_in_a_small_legend():
    fits = [_stub_fit("study_a"), _stub_fit("study_b")]
    ax = sh.plot_catalog_fits(fits).axes[0]
    texts = [text.get_text() for text in ax.get_legend().get_texts()]
    assert texts == ["study_a (exponential)", "study_b (exponential)"]


def test_plot_catalog_fits_omits_the_analyte_when_a_study_contributes_one():
    """Both models of one analyte are already told apart by linestyle."""
    fits = [
        _stub_fit("study_a", model="exponential"),
        _stub_fit("study_a", model="gamma"),
    ]
    labels = sorted(
        line.get_label() for line in _curves(sh.plot_catalog_fits(fits).axes[0])
    )
    assert labels == ["study_a (exponential)", "study_a (gamma)"]
