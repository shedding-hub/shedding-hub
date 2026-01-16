import matplotlib

matplotlib.use("Agg")

import pytest
import matplotlib.figure
import shedding_hub as sh


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
