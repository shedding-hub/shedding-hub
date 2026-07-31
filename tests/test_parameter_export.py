import json

import numpy as np
import pytest

from shedding_hub.shedding_catalog import SheddingCatalog
from shedding_hub.shedding_export import catalog_to_records


@pytest.fixture
def two_model_catalog(make_synthetic_dataset):
    from shedding_hub.shedding_catalog import fit_shedding_models

    dataset = make_synthetic_dataset(
        "gamma",
        [0.0, np.log(2.0), np.log(12.0)],
        np.diag([0.04, 0.04, 0.09]),
        n_subjects=12,
        seed=3,
    )
    return fit_shedding_models([dataset], models=("exponential", "gamma"))


def test_one_record_per_fit(two_model_catalog):
    records = catalog_to_records(two_model_catalog)
    assert len(records) == len(two_model_catalog.fits)
    keys = {(r["dataset_id"], r["analyte"], r["model"]) for r in records}
    assert keys == {(f.dataset_id, f.analyte, f.model) for f in two_model_catalog.fits}


def test_parameters_are_named_by_model(two_model_catalog):
    """A reader should not have to know which model produced a row to read it."""
    records = {r["model"]: r for r in catalog_to_records(two_model_catalog)}
    assert set(records["exponential"]["parameters"]) == {"a0", "c0"}
    assert set(records["gamma"]["parameters"]) == {"a0", "b0", "c0"}


def test_record_carries_enough_to_simulate_without_refitting(two_model_catalog):
    """The point of the export: reuse, not just reading."""
    fit = two_model_catalog.fits[0]
    record = next(
        r
        for r in catalog_to_records(two_model_catalog)
        if (r["dataset_id"], r["analyte"], r["model"])
        == (fit.dataset_id, fit.analyte, fit.model)
    )
    population = record["population"]
    assert population["coordinates"] == list(fit.population_coords)
    np.testing.assert_allclose(population["mean"], fit.population_mean)
    np.testing.assert_allclose(population["covariance"], fit.population_cov)
    assert record["measurement_error_sd"] == pytest.approx(fit.sigma)
    assert record["censoring_limit_log10"] == pytest.approx(fit.censoring_limit)


def test_records_are_json_serializable(two_model_catalog):
    """numpy floats and arrays must not leak into the file."""
    text = json.dumps(catalog_to_records(two_model_catalog))
    assert json.loads(text)


def test_empty_catalog_exports_nothing():
    assert catalog_to_records(SheddingCatalog()) == []
