"""
A browsable catalog of fitted shedding models, and cross-study ensembles.

The catalog is the surface a modeller browses: one row per (dataset, analyte,
model), summarising each fit by its median individual, so a row reads as "peaks
day 4.2 at 6.8 log10 gc/mL, declines with a 1.5 day half-life" rather than as raw
log-parameters.
"""

import pathlib
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import yaml

from .shedding_fit import SheddingDataError, SheddingFit, fit_shedding_model
from .shedding_models import MODELS, PARAM_NAMES

CATALOG_PATH = pathlib.Path(__file__).parent / "data" / "shedding_catalog.yaml"

_KEY_COLUMNS = (
    "dataset_id",
    "analyte",
    "biomarker",
    "specimen",
    "reference_event",
    "unit",
    "gene_target",
    "dose",
    "vaccine_type",
    "model",
)

# Every column a fit's row can produce, in display order. Both models' medians
# (a/b/c) are always present so that an empty catalog, an exponential-only
# catalog (no b0), and a mixed catalog all present exactly the same columns —
# downstream code should never need to check which model produced a row before
# indexing into it.
_TABLE_COLUMNS = _KEY_COLUMNS + (
    "n_subjects",
    "n_measurements",
    "pct_censored",
    "pct_subjects_with_rise",
    "a_median",
    "b_median",
    "c_median",
    "sigma",
    "peak_day",
    "peak_log10",
    "half_life_days",
    "aic",
    "converged",
)


def fit_to_row(fit: SheddingFit) -> dict:
    """
    Summarize a fit as one table row describing its median individual.

    Because ``theta = log(params)`` is normal, the parameters are lognormal and
    ``exp(mu)`` is exactly their median. These are therefore labelled ``_median``,
    which is the accurate name rather than a compromise.

    All three of ``a_median``, ``b_median``, ``c_median`` are always present,
    even for models that lack one of them (the exponential model has no
    ``b0``): the missing one is ``NaN`` rather than the key being absent, so a
    row's schema does not depend on which model produced it.
    """
    row = {column: getattr(fit, column) for column in _KEY_COLUMNS}
    row.update({"a_median": np.nan, "b_median": np.nan, "c_median": np.nan})
    medians = fit.median_params
    for name, value in zip(PARAM_NAMES[fit.model], medians):
        row[f"{name[0]}_median"] = float(value)
    row.update(
        {
            "n_subjects": fit.n_subjects,
            "n_measurements": fit.n_measurements,
            "pct_censored": (
                100.0 * fit.n_censored / fit.n_measurements
                if fit.n_measurements
                else np.nan
            ),
            # Surfaced for both models so the gamma gate is auditable from the
            # table rather than implicit in which rows are missing. On an
            # exponential row it is informational: a low value there is the
            # normal case, not a warning.
            "pct_subjects_with_rise": fit.pct_subjects_with_rise,
            "sigma": fit.sigma,
            "peak_day": fit.peak_day,
            "peak_log10": fit.peak_log10,
            "half_life_days": fit.half_life_days,
            "aic": fit.aic,
            "converged": fit.converged,
        }
    )
    return row


def _fits_to_frame(fits: list[SheddingFit]) -> pd.DataFrame:
    if not fits:
        return pd.DataFrame(columns=list(_TABLE_COLUMNS))
    return pd.DataFrame([fit_to_row(fit) for fit in fits]).reindex(
        columns=list(_TABLE_COLUMNS)
    )


def _fit_to_payload(fit: SheddingFit) -> dict:
    """Serialize a fit, omitting per-subject parameters to keep the file small."""
    payload = {column: getattr(fit, column) for column in _KEY_COLUMNS}
    payload.update(
        {
            "method": fit.method,
            "population_mean": [float(v) for v in fit.population_mean],
            "population_cov": [[float(v) for v in row] for row in fit.population_cov],
            "sigma": float(fit.sigma),
            "censoring_limit": float(fit.censoring_limit),
            "n_subjects": int(fit.n_subjects),
            "n_measurements": int(fit.n_measurements),
            "n_censored": int(fit.n_censored),
            "n_excluded_subjects": int(fit.n_excluded_subjects),
            "n_degenerate_subjects": int(fit.n_degenerate_subjects),
            "pct_subjects_with_rise": float(fit.pct_subjects_with_rise),
            "n_dropped_measurements": int(fit.n_dropped_measurements),
            "converged": bool(fit.converged),
            "log_likelihood": float(fit.log_likelihood),
            "aic": float(fit.aic),
        }
    )
    return payload


def _fit_from_payload(payload: dict) -> SheddingFit:
    return SheddingFit(
        model=payload["model"],
        method=payload.get("method", "mle"),
        population_mean=np.asarray(payload["population_mean"], dtype=float),
        population_cov=np.asarray(payload["population_cov"], dtype=float),
        sigma=float(payload["sigma"]),
        subject_params=None,
        censoring_limit=float(payload["censoring_limit"]),
        dataset_id=payload["dataset_id"],
        analyte=payload["analyte"],
        biomarker=payload.get("biomarker"),
        specimen=payload.get("specimen"),
        reference_event=payload.get("reference_event"),
        unit=payload.get("unit"),
        gene_target=payload.get("gene_target"),
        dose=payload.get("dose"),
        vaccine_type=payload.get("vaccine_type"),
        n_subjects=int(payload["n_subjects"]),
        n_measurements=int(payload["n_measurements"]),
        n_censored=int(payload["n_censored"]),
        n_excluded_subjects=int(payload["n_excluded_subjects"]),
        # Defaulted, not required: catalogs written before degeneracy detection
        # existed have no such key, and every fit in them predates the concept.
        n_degenerate_subjects=int(payload.get("n_degenerate_subjects", 0)),
        # Likewise defaulted: NaN reads as "this catalog predates the rise
        # gate", which is honest, where 0.0 would assert that no subject rose.
        pct_subjects_with_rise=float(
            payload.get("pct_subjects_with_rise", float("nan"))
        ),
        n_dropped_measurements=int(payload["n_dropped_measurements"]),
        converged=bool(payload["converged"]),
        log_likelihood=float(payload["log_likelihood"]),
        aic=float(payload["aic"]),
    )


@dataclass
class SheddingCatalog:
    """A collection of fitted models with a browsable summary table."""

    fits: list[SheddingFit] = field(default_factory=list)
    skipped: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(
            columns=["dataset_id", "analyte", "model", "reason", "message"]
        )
    )

    @property
    def table(self) -> pd.DataFrame:
        """One row per fit, summarising its median individual."""
        return _fits_to_frame(self.fits)

    def select(self, **keys) -> SheddingFit:
        """
        Return the single fit matching ``keys``.

        Raises:
            ValueError: If no fit matches, or if more than one does. Never picks
                silently — the error lists the candidates and the columns that
                would tell them apart.
        """
        matches = [
            fit
            for fit in self.fits
            if all(getattr(fit, key, None) == value for key, value in keys.items())
        ]
        if not matches:
            raise ValueError(
                f"select({keys}) matched no fits. "
                f"Available combinations are listed in `catalog.table`; "
                f"{len(self.fits)} fit(s) are loaded."
            )
        if len(matches) > 1:
            candidates = _fits_to_frame(matches)
            distinguishing = [
                column
                for column in _KEY_COLUMNS
                if column not in keys and candidates[column].nunique() > 1
            ]
            raise ValueError(
                f"select({keys}) matched {len(matches)} fits. Narrow it with "
                f"{distinguishing or 'more specific keys'}. Candidates:\n"
                f"{candidates[list(_KEY_COLUMNS)].to_string(index=False)}"
            )
        return matches[0]

    def ensemble(
        self, *, dataset_ids=None, weights="n_subjects", method="mixture", **keys
    ):
        """Build a ``SheddingEnsemble`` from the matching fits. See Task 7."""
        from .shedding_ensemble import build_ensemble

        return build_ensemble(
            self, dataset_ids=dataset_ids, weights=weights, method=method, **keys
        )

    def to_dict(self) -> dict:
        return {
            "fits": [_fit_to_payload(fit) for fit in self.fits],
            "skipped": self.skipped.to_dict(orient="records"),
        }

    @classmethod
    def from_dict(cls, payload: dict) -> "SheddingCatalog":
        skipped = pd.DataFrame(
            payload.get("skipped") or [],
            columns=["dataset_id", "analyte", "model", "reason", "message"],
        )
        return cls(
            fits=[_fit_from_payload(item) for item in payload.get("fits", [])],
            skipped=skipped,
        )


def fit_shedding_models(
    datasets,
    *,
    models=MODELS,
    min_observations: int | None = None,
) -> SheddingCatalog:
    """
    Fit every analyte of every dataset, for every requested model.

    Analytes that cannot be fitted are recorded in ``catalog.skipped`` with a
    reason rather than raising, so one unsuitable analyte does not abort a
    repository-wide build — and so a missing study reads as unsuitable rather
    than as a bug.

    Args:
        datasets: Dataset dictionaries from ``load_dataset``.
        models: Model names to fit. Defaults to both.
        min_observations: Passed through to the fitter.

    Returns:
        A ``SheddingCatalog``.
    """
    fits: list[SheddingFit] = []
    skipped: list[dict] = []

    for dataset in datasets:
        dataset_id = dataset.get("dataset_id", "unknown")
        for analyte in dataset.get("analytes", {}):
            for model in models:
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", UserWarning)
                        fits.append(
                            fit_shedding_model(
                                dataset,
                                analyte=analyte,
                                model=model,
                                min_observations=min_observations,
                            )
                        )
                except SheddingDataError as error:
                    skipped.append(
                        {
                            "dataset_id": dataset_id,
                            "analyte": analyte,
                            "model": model,
                            "reason": error.reason,
                            "message": str(error),
                        }
                    )
                except (ValueError, np.linalg.LinAlgError) as error:
                    skipped.append(
                        {
                            "dataset_id": dataset_id,
                            "analyte": analyte,
                            "model": model,
                            "reason": "did_not_converge",
                            "message": str(error),
                        }
                    )

    return SheddingCatalog(
        fits=fits,
        skipped=pd.DataFrame(
            skipped, columns=["dataset_id", "analyte", "model", "reason", "message"]
        ),
    )


def load_shedding_catalog(path: str | None = None) -> SheddingCatalog:
    """
    Load the catalog of precomputed estimates shipped with the package.

    Args:
        path: Optional path to a catalog YAML. Defaults to the shipped file.

    Returns:
        A ``SheddingCatalog``. Loaded fits carry ``subject_params is None``
        because per-subject values are not serialized; everything needed to
        simulate (``mu``, ``Sigma``, ``sigma``) is present.
    """
    catalog_path = pathlib.Path(path) if path else CATALOG_PATH
    if not catalog_path.is_file():
        raise FileNotFoundError(
            f"No shedding catalog at {catalog_path}. Run `make catalog` to build it."
        )
    with catalog_path.open(encoding="utf-8") as stream:
        payload = yaml.safe_load(stream)
    return SheddingCatalog.from_dict(payload)
