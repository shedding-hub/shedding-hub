import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from shedding_hub.shedding_fit import (
    SheddingDataError,
    prepare_observations,
)


@pytest.fixture
def simple_dataset():
    """Two subjects, one analyte, with a censored point and a qualitative one."""
    return {
        "dataset_id": "test_study",
        "analytes": {
            "stool": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": 1e6},
                    {"analyte": "stool", "time": 2, "value": 1e5},
                    {"analyte": "stool", "time": 3, "value": "negative"},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": 1e7},
                    {"analyte": "stool", "time": 2, "value": 1e6},
                    {"analyte": "stool", "time": 3, "value": "positive"},
                    {"analyte": "stool", "time": "unknown", "value": 1e4},
                ]
            },
        ],
    }


def test_extracts_values_on_log10_scale(simple_dataset):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.n_subjects == 2
    np.testing.assert_allclose(obs.values[0], 6.0)
    np.testing.assert_allclose(obs.values[1], 5.0)


def test_negative_becomes_censored(simple_dataset):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censored.sum() == 1
    assert np.isnan(obs.values[obs.censored][0])


def test_censoring_limit_uses_declared_loq_when_below_smallest_positive(
    simple_dataset,
):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    # LOQ 100 -> log10 = 2, which is below the smallest observed positive (5.0)
    assert obs.censoring_limit == pytest.approx(2.0)


def test_censoring_limit_falls_back_below_smallest_positive(simple_dataset):
    # Declare a limit above every observed value; the fallback must kick in.
    simple_dataset["analytes"]["stool"]["limit_of_quantification"] = 1e8
    with pytest.warns(UserWarning, match="censoring limit"):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censoring_limit == pytest.approx(5.0 - 0.01)


def test_censoring_limit_uses_lod_when_loq_is_unusable(simple_dataset):
    # LOQ is "unknown" (unusable); LOD is a valid number below the smallest
    # observed positive (5.0), so the declared LOD should resolve cleanly,
    # with no fallback-below-smallest-positive warning.
    simple_dataset["analytes"]["stool"]["limit_of_quantification"] = "unknown"
    simple_dataset["analytes"]["stool"]["limit_of_detection"] = 10
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censoring_limit == pytest.approx(1.0)


def test_qualitative_and_unknown_time_are_dropped_with_warning(simple_dataset):
    with pytest.warns(UserWarning):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    # 7 measurements total, minus 1 qualitative "positive", minus 1 unknown time
    assert obs.times.size == 5
    assert obs.n_dropped_measurements == 2


def test_ct_analyte_is_rejected():
    dataset = {
        "dataset_id": "ct_study",
        "analytes": {
            "swab": {
                "specimen": "saliva",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "cycle threshold",
                "limit_of_quantification": "unknown",
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "swab", "time": 1, "value": 20.0},
                    {"analyte": "swab", "time": 2, "value": 25.0},
                ]
            }
        ],
    }
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(dataset, "swab", "exponential")
    assert excinfo.value.reason == "ct_units"


def test_non_pathogen_biomarker_is_rejected():
    # Plenty of clean, numeric, positive-time data for both subjects, so a
    # rejection here can only be the biomarker check, not a data shortage.
    dataset = {
        "dataset_id": "crassphage_study",
        "analytes": {
            "stool_crAssphage": {
                "specimen": "stool",
                "biomarker": "crAssphage",
                "reference_event": "symptom onset",
                "unit": "gc/dry gram",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool_crAssphage", "time": 1, "value": 1e6},
                    {"analyte": "stool_crAssphage", "time": 2, "value": 1e5},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool_crAssphage", "time": 1, "value": 1e7},
                    {"analyte": "stool_crAssphage", "time": 2, "value": 1e6},
                ]
            },
        ],
    }
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(dataset, "stool_crAssphage", "exponential")
    assert excinfo.value.reason == "non_pathogen_biomarker"


def test_vaccine_strain_biomarker_is_not_rejected():
    # Live-attenuated vaccine shedding has a real trajectory (reference event
    # "vaccination"), unlike the fecal-strength/normalization indicators in
    # NON_PATHOGEN_BIOMARKERS. Guards against that set later being "tidied"
    # into a broader substring/category match that would catch this too.
    dataset = {
        "dataset_id": "vaccine_study",
        "analytes": {
            "stool_rotavirus": {
                "specimen": "stool",
                "biomarker": "rotavirus vaccine",
                "reference_event": "vaccination",
                "unit": "gc/mL",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool_rotavirus", "time": 1, "value": 1e6},
                    {"analyte": "stool_rotavirus", "time": 2, "value": 1e5},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool_rotavirus", "time": 1, "value": 1e7},
                    {"analyte": "stool_rotavirus", "time": 2, "value": 1e6},
                ]
            },
        ],
    }
    obs = prepare_observations(dataset, "stool_rotavirus", "exponential")
    assert obs.n_subjects == 2


def test_no_positive_measurements_raises():
    # Both subjects have enough censored measurements to clear
    # min_observations on their own, so they are retained rather than
    # excluded — the failure must come from having no positives at all.
    dataset = {
        "dataset_id": "all_negative_study",
        "analytes": {
            "stool": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": "negative"},
                    {"analyte": "stool", "time": 2, "value": "negative"},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": "negative"},
                    {"analyte": "stool", "time": 2, "value": "negative"},
                ]
            },
        ],
    }
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(dataset, "stool", "exponential")
    assert excinfo.value.reason == "no_positive_measurements"


def test_gamma_drops_non_positive_times(simple_dataset):
    simple_dataset["participants"][0]["measurements"].append(
        {"analyte": "stool", "time": 0, "value": 1e6}
    )
    simple_dataset["participants"][0]["measurements"].append(
        {"analyte": "stool", "time": -2, "value": 1e6}
    )
    with pytest.warns(UserWarning):
        gamma_obs = prepare_observations(simple_dataset, "stool", "gamma")
    exponential_obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert (gamma_obs.times > 0).all()
    # Two effects separate exponential from gamma here: the two appended
    # non-positive times are dropped outright under gamma, and subject 2 (only
    # 2 usable measurements) falls below gamma's 3-observation minimum and is
    # excluded entirely, while it clears exponential's 2-observation minimum.
    assert exponential_obs.times.size == gamma_obs.times.size + 4


def test_subjects_below_min_observations_are_excluded(simple_dataset):
    simple_dataset["participants"].append(
        {"measurements": [{"analyte": "stool", "time": 1, "value": 1e5}]}
    )
    with pytest.warns(UserWarning, match="excluded"):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.n_subjects == 2
    assert obs.n_excluded_subjects == 1


def test_no_usable_subject_raises(simple_dataset):
    for participant in simple_dataset["participants"]:
        participant["measurements"] = participant["measurements"][:1]
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(simple_dataset, "stool", "exponential")
    assert excinfo.value.reason == "too_few_subjects"


def test_unknown_analyte_raises(simple_dataset):
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(simple_dataset, "sputum", "exponential")
    assert excinfo.value.reason == "unknown_analyte"


def test_subject_index_is_contiguous_after_exclusions(simple_dataset):
    simple_dataset["participants"].insert(
        0, {"measurements": [{"analyte": "stool", "time": 1, "value": 1e5}]}
    )
    with pytest.warns(UserWarning):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert set(obs.subject_index.tolist()) == {0, 1}
