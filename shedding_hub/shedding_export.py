"""
Flatten a fitted catalog into plain records for reference and reuse.

The catalog's own YAML is the canonical store; this is the browsing and
interchange view of it. Each record is self-contained: alongside the fitted
parameters it carries the population mean and covariance, the measurement-error
SD and the censoring limit, which is everything ``simulate_shedding`` needs. A
reader can therefore reuse an estimate without the package, and without
refitting anything.
"""

import numpy as np

from .shedding_models import PARAM_NAMES


def _plain(value):
    """Convert numpy scalars and arrays to JSON-safe Python objects."""
    if isinstance(value, np.ndarray):
        return [_plain(item) for item in value.tolist()]
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        # JSON has no NaN or Infinity. Emitting null keeps the file readable by
        # any parser rather than only by Python's permissive one.
        return None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def fit_to_record(fit) -> dict:
    """
    Represent one fit as a nested, JSON-safe record.

    Parameters are keyed by name rather than positionally, so a reader never has
    to know which model produced a record before indexing into it. The keys
    differ between models because the models differ: the exponential has no
    rise, and only ``gamma_shifted`` has an onset.
    """
    medians = fit.median_params
    return {
        "dataset_id": fit.dataset_id,
        "analyte": fit.analyte,
        "model": fit.model,
        "biomarker": fit.biomarker,
        "specimen": fit.specimen,
        "unit": fit.unit,
        "reference_event": fit.reference_event,
        "gene_target": fit.gene_target,
        # The median individual, in the model's own parameters.
        "parameters": {
            name: _plain(value) for name, value in zip(PARAM_NAMES[fit.model], medians)
        },
        # The same individual described in interpretable terms.
        "summary": {
            "peak_day": _plain(fit.peak_day),
            "peak_log10": _plain(fit.peak_log10),
            "half_life_days": _plain(fit.half_life_days),
        },
        # Everything simulate_shedding needs, so a record can be reused as-is.
        "population": {
            "coordinates": list(fit.population_coords),
            "mean": _plain(fit.population_mean),
            "covariance": _plain(fit.population_cov),
        },
        "measurement_error_sd": _plain(fit.sigma),
        "censoring_limit_log10": _plain(fit.censoring_limit),
        # What the estimate rests on. A parameter cannot be judged without it.
        "data": {
            "n_subjects": _plain(fit.n_subjects),
            "n_measurements": _plain(fit.n_measurements),
            "pct_censored": _plain(
                100.0 * fit.n_censored / fit.n_measurements
                if fit.n_measurements
                else np.nan
            ),
            "n_degenerate_subjects": _plain(fit.n_degenerate_subjects),
            "median_first_observed_day": _plain(fit.median_first_observed_day),
            "pct_subjects_with_rise": _plain(fit.pct_subjects_with_rise),
        },
        "fit": {"aic": _plain(fit.aic), "converged": _plain(fit.converged)},
    }


def catalog_to_records(catalog) -> list:
    """One record per fit, ordered by dataset, analyte, then model."""
    return [
        fit_to_record(fit)
        for fit in sorted(
            catalog.fits, key=lambda f: (f.dataset_id, f.analyte, f.model)
        )
    ]
