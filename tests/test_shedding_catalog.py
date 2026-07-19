import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from shedding_hub.shedding_catalog import SheddingCatalog, fit_shedding_models


@pytest.fixture
def two_study_catalog(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    cov = np.diag([0.04, 0.04])
    a = make_synthetic_dataset(
        "exponential", mu, cov, n_subjects=20, seed=1, dataset_id="study_a"
    )
    b = make_synthetic_dataset(
        "exponential", mu, cov, n_subjects=20, seed=2, dataset_id="study_b"
    )
    return fit_shedding_models([a, b], models=("exponential",))


def test_table_has_one_row_per_fit(two_study_catalog):
    table = two_study_catalog.table
    assert len(table) == 2
    assert set(table["dataset_id"]) == {"study_a", "study_b"}


def test_table_reports_medians_and_derived_quantities(two_study_catalog):
    table = two_study_catalog.table
    for column in [
        "dataset_id",
        "analyte",
        "biomarker",
        "specimen",
        "reference_event",
        "unit",
        "model",
        "n_subjects",
        "n_measurements",
        "pct_censored",
        "a_median",
        "sigma",
        "peak_day",
        "peak_log10",
        "half_life_days",
        "aic",
        "converged",
    ]:
        assert column in table.columns
    fit = two_study_catalog.fits[0]
    row = table.iloc[0]
    assert row["a_median"] == pytest.approx(np.exp(fit.population_mean[0]))
    assert row["half_life_days"] == pytest.approx(np.log(2.0) / row["a_median"])


def test_gamma_table_has_b_median_and_peak_day(make_synthetic_dataset):
    mu = np.array([np.log(0.5), np.log(2.0), np.log(12.0)])
    dataset = make_synthetic_dataset("gamma", mu, np.diag([0.04] * 3), n_subjects=20)
    catalog = fit_shedding_models([dataset], models=("gamma",))
    row = catalog.table.iloc[0]
    assert row["peak_day"] == pytest.approx(row["b_median"] / row["a_median"])


def test_select_returns_one_fit(two_study_catalog):
    fit = two_study_catalog.select(dataset_id="study_a")
    assert fit.dataset_id == "study_a"


def test_select_raises_on_ambiguous_match(two_study_catalog):
    with pytest.raises(ValueError, match="matched 2"):
        two_study_catalog.select(analyte="stool")


def test_select_raises_on_no_match(two_study_catalog):
    with pytest.raises(ValueError, match="matched no"):
        two_study_catalog.select(dataset_id="study_z")


def test_ct_analyte_is_recorded_in_skipped():
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
    catalog = fit_shedding_models([dataset], models=("exponential",))
    assert catalog.table.empty
    assert (catalog.skipped["reason"] == "ct_units").all()
    assert set(catalog.skipped["dataset_id"]) == {"ct_study"}


def test_cross_sectional_study_is_skipped():
    dataset = {
        "dataset_id": "cross_sectional",
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
            {"measurements": [{"analyte": "stool", "time": 0, "value": 1e5}]}
            for _ in range(30)
        ],
    }
    catalog = fit_shedding_models([dataset], models=("exponential",))
    assert catalog.table.empty
    assert (catalog.skipped["reason"] == "too_few_subjects").all()


def test_empty_table_has_same_columns_as_populated_table(two_study_catalog):
    empty = SheddingCatalog().table
    populated = two_study_catalog.table
    assert list(empty.columns) == list(populated.columns)


def test_exponential_only_catalog_still_has_b_median_column(two_study_catalog):
    table = two_study_catalog.table
    assert "b_median" in table.columns
    assert table["b_median"].isna().all()


def test_round_trip_serialization(two_study_catalog, tmp_path):
    payload = two_study_catalog.to_dict()
    restored = SheddingCatalog.from_dict(payload)
    assert len(restored.fits) == len(two_study_catalog.fits)
    original = two_study_catalog.fits[0]
    copy = restored.select(dataset_id=original.dataset_id, model=original.model)
    np.testing.assert_allclose(copy.population_mean, original.population_mean)
    np.testing.assert_allclose(copy.population_cov, original.population_cov)
    assert copy.sigma == pytest.approx(original.sigma)
    assert copy.censoring_limit == pytest.approx(original.censoring_limit)
    assert copy.model == original.model
    assert copy.biomarker == original.biomarker
    assert copy.specimen == original.specimen
    assert copy.reference_event == original.reference_event
    assert copy.unit == original.unit
    assert copy.aic == pytest.approx(original.aic)
    assert copy.converged == original.converged
    assert copy.n_subjects == original.n_subjects
    assert copy.subject_params is None


def test_restored_fit_can_still_simulate(two_study_catalog):
    from shedding_hub.shedding_simulate import simulate_shedding

    restored = SheddingCatalog.from_dict(two_study_catalog.to_dict())
    traj = simulate_shedding(
        restored.fits[0], n_individuals=5, times=[1.0, 2.0], seed=0
    )
    assert len(traj) == 10
