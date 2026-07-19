import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from shedding_hub.shedding_catalog import SheddingCatalog, fit_shedding_models
from shedding_hub.shedding_fit import SheddingFit


def _fit_with_degenerates(n_degenerate):
    """A directly-constructed fit carrying a chosen ``n_degenerate_subjects``.

    Built by hand rather than fitted so the count is exactly what the test
    intends, independent of what any optimizer happens to do.
    """
    return SheddingFit(
        model="exponential",
        method="mle",
        population_mean=np.array([np.log(0.6), np.log(18.0)]),
        population_cov=np.diag([0.04, 0.04]),
        sigma=0.3,
        subject_params=None,
        censoring_limit=2.0,
        dataset_id="synthetic",
        analyte="stool",
        biomarker="SARS-CoV-2",
        specimen="stool",
        reference_event="symptom onset",
        unit="gc/mL",
        gene_target=None,
        dose=None,
        vaccine_type=None,
        n_subjects=10,
        n_measurements=100,
        n_censored=5,
        n_excluded_subjects=0,
        n_dropped_measurements=0,
        converged=True,
        log_likelihood=-10.0,
        aic=42.0,
        n_degenerate_subjects=n_degenerate,
        pct_subjects_with_rise=62.5,
        median_first_observed_day=3.0,
    )


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
    assert copy.n_excluded_subjects == original.n_excluded_subjects
    assert copy.n_degenerate_subjects == original.n_degenerate_subjects
    assert copy.pct_subjects_with_rise == pytest.approx(
        original.pct_subjects_with_rise, nan_ok=True
    )
    assert copy.median_first_observed_day == pytest.approx(
        original.median_first_observed_day, nan_ok=True
    )
    assert copy.n_dropped_measurements == original.n_dropped_measurements
    assert copy.n_measurements == original.n_measurements
    assert copy.n_censored == original.n_censored
    assert copy.log_likelihood == pytest.approx(original.log_likelihood)
    assert copy.method == original.method
    assert copy.subject_params is None


def test_round_trip_preserves_a_non_zero_degenerate_count():
    """``n_degenerate_subjects`` must survive serialization on its own merits.

    ``test_round_trip_serialization`` compares fits whose count happens to be
    zero, which a field dropped from the payload would also satisfy — the
    default is zero. Round-tripping a non-zero count is what actually pins it.
    """
    fit = _fit_with_degenerates(3)
    restored = SheddingCatalog.from_dict(SheddingCatalog(fits=[fit]).to_dict())
    assert restored.fits[0].n_degenerate_subjects == 3


def test_fit_from_payload_defaults_missing_degenerate_count():
    """Catalogs written before degeneracy detection existed must still load."""
    payload = SheddingCatalog(fits=[_fit_with_degenerates(3)]).to_dict()
    del payload["fits"][0]["n_degenerate_subjects"]
    assert SheddingCatalog.from_dict(payload).fits[0].n_degenerate_subjects == 0


def test_round_trip_preserves_the_rise_percentage():
    fit = _fit_with_degenerates(0)
    restored = SheddingCatalog.from_dict(SheddingCatalog(fits=[fit]).to_dict())
    assert restored.fits[0].pct_subjects_with_rise == pytest.approx(62.5)


def test_fit_from_payload_defaults_missing_rise_percentage_to_nan():
    """A catalog predating the rise gate must read as unknown, not as zero.

    Zero would assert that no subject rose, which is a claim the old catalog
    never made.
    """
    payload = SheddingCatalog(fits=[_fit_with_degenerates(0)]).to_dict()
    del payload["fits"][0]["pct_subjects_with_rise"]
    assert np.isnan(SheddingCatalog.from_dict(payload).fits[0].pct_subjects_with_rise)


def test_round_trip_preserves_the_median_first_observed_day():
    fit = _fit_with_degenerates(0)
    restored = SheddingCatalog.from_dict(SheddingCatalog(fits=[fit]).to_dict())
    assert restored.fits[0].median_first_observed_day == pytest.approx(3.0)


def test_fit_from_payload_defaults_missing_median_first_observed_day_to_nan():
    """A catalog predating this column must read as unknown, not as day zero.

    Zero would be the strongest possible claim -- that the study sampled the
    reference event itself -- which is exactly the claim being audited.
    """
    payload = SheddingCatalog(fits=[_fit_with_degenerates(0)]).to_dict()
    del payload["fits"][0]["median_first_observed_day"]
    assert np.isnan(
        SheddingCatalog.from_dict(payload).fits[0].median_first_observed_day
    )


def test_catalog_refuses_a_two_subject_fit_as_a_population(make_synthetic_dataset):
    """The population gate is applied by the builder, so no such row is published."""
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=2
    )
    catalog = fit_shedding_models([dataset], models=("exponential",))
    assert catalog.table.empty
    assert (catalog.skipped["reason"] == "too_few_subjects_for_population").all()


def test_catalog_publishes_a_three_subject_fit(make_synthetic_dataset):
    """One subject more than parameters is published, so the gate is not stricter."""
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=3
    )
    catalog = fit_shedding_models([dataset], models=("exponential",))
    assert len(catalog.table) == 1
    assert catalog.skipped.empty


def test_table_exposes_the_median_first_observed_day(make_synthetic_dataset):
    """The extrapolation behind peak_log10 must be auditable from the table."""
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=20
    )
    catalog = fit_shedding_models([dataset], models=("exponential",))
    row = catalog.table.iloc[0]
    assert "median_first_observed_day" in catalog.table.columns
    # The fixture samples days 1..14, so nothing was observed at t = 0 -- which
    # is precisely where the exponential model reports its peak.
    assert row["median_first_observed_day"] == pytest.approx(1.0)


def test_table_exposes_the_rise_percentage_for_both_models(make_synthetic_dataset):
    """The gamma gate must be auditable from the table, not just from absences."""
    mu = np.array([np.log(0.5), np.log(2.0), np.log(12.0)])
    dataset = make_synthetic_dataset("gamma", mu, np.diag([0.04] * 3), n_subjects=20)
    catalog = fit_shedding_models([dataset], models=("gamma", "exponential"))
    table = catalog.table
    assert "pct_subjects_with_rise" in table.columns
    # This synthetic truth peaks at day 4 within a 1..14 window, so essentially
    # every subject rises -- and the exponential row reports it too.
    assert (table["pct_subjects_with_rise"] > 50).all()
    assert set(table["model"]) == {"gamma", "exponential"}


def test_restored_fit_can_still_simulate(two_study_catalog):
    from shedding_hub.shedding_simulate import simulate_shedding

    restored = SheddingCatalog.from_dict(two_study_catalog.to_dict())
    traj = simulate_shedding(
        restored.fits[0], n_individuals=5, times=[1.0, 2.0], seed=0
    )
    assert len(traj) == 10


def test_shipped_catalog_covers_every_dataset():
    """CI staleness check: adding a dataset without regenerating must fail."""
    import pathlib

    from shedding_hub.shedding_catalog import load_shedding_catalog

    data_dir = pathlib.Path(__file__).parent.parent / "data"
    on_disk = {
        path.name
        for path in data_dir.iterdir()
        if path.is_dir() and not path.name.startswith(".")
    }
    catalog = load_shedding_catalog()
    accounted = set(catalog.table["dataset_id"]) | set(catalog.skipped["dataset_id"])
    missing = on_disk - accounted
    assert not missing, (
        f"Datasets absent from the shipped catalog: {sorted(missing)}. "
        "Run `make catalog` to regenerate."
    )
