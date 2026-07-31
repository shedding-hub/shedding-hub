import json
import pathlib

import numpy as np
import pandas as pd
import pytest

from shedding_hub.shedding_catalog import SheddingCatalog, load_shedding_catalog
from shedding_hub.shedding_export import catalog_to_records

DOCS = pathlib.Path(__file__).parent.parent / "docs"
JSON_PATH = DOCS / "shedding_parameters.json"
CSV_PATH = DOCS / "shedding_parameters.csv"


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


def test_shipped_json_matches_the_shipped_catalog():
    """
    Staleness check: regenerating the catalog without re-exporting must fail.

    The published table is meant to be reused without this package, so a reader
    has no way to notice it has drifted from the catalog it claims to describe.
    The catalog has such a guard; the file derived from it needs one too.
    """
    catalog = load_shedding_catalog()
    # Round-tripping both sides through the same serializer makes this an
    # exact comparison rather than a float-repr one.
    fresh = json.loads(json.dumps(catalog_to_records(catalog)))
    committed = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    assert committed["fits"] == fresh, (
        "docs/shedding_parameters.json disagrees with the shipped catalog. "
        "Run `make parameters` to regenerate it."
    )
    assert committed["n_fits"] == len(fresh)
    assert committed["n_datasets"] == len({r["dataset_id"] for r in fresh})
    assert committed["models"] == sorted({r["model"] for r in fresh})


def test_shipped_csv_matches_the_shipped_catalog(tmp_path):
    """The browsing view goes stale the same way, and by the same omission."""
    catalog = load_shedding_catalog()
    fresh_path = tmp_path / "fresh.csv"
    # Written and read back exactly as the committed file was, so neither line
    # endings nor inferred dtypes can masquerade as a real difference.
    catalog.table.to_csv(fresh_path, index=False)

    pd.testing.assert_frame_equal(
        pd.read_csv(fresh_path),
        pd.read_csv(CSV_PATH),
        obj="docs/shedding_parameters.csv (run `make parameters`)",
    )
