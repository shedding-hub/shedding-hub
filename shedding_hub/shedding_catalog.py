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

from .shedding_fit import (
    SheddingDataError,
    SheddingFit,
    fit_shedding_model,
    require_estimable_population,
)
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
    "median_first_observed_day",
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

    Three columns qualify the rest, and a row is easy to misread without them:

    ``median_first_observed_day``
        Median across retained subjects of each subject's own first sampling
        day — deliberately the median rather than the earliest, since one
        early-enrolled subject should not make a late-starting study look
        well-observed. ``peak_log10`` is evaluated at the peak — ``t = 0`` for
        the exponential model, since that model only decays — so when this
        column is well above zero the peak is a backward extrapolation to a time
        most subjects were never observed at, not a measured concentration. A
        large ``median_first_observed_day`` beside a large ``peak_log10`` means
        precisely that, and should not be read as the study having detected
        ``10 ** peak_log10`` of anything.

    ``pct_censored``
        Share of measurements that were below the limit of detection. A very
        high value (the repository reaches the low 90s) means the analyte is
        rarely detected at all, so the row's estimates describe mostly-censored
        data — a low ``peak_log10`` there is the fitter reporting honestly on an
        analyte that is almost never found, not a defect.

    ``pct_subjects_with_rise``
        Share of adequately-sampled subjects whose highest reading was not their
        first. The gamma model is refused below 50%; on an exponential row it is
        informational, and a low value is the normal case rather than a warning.
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
            "median_first_observed_day": fit.median_first_observed_day,
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
    """Serialize a fit, omitting per-subject parameters to keep the file small.

    Delegates to ``SheddingFit.to_dict``, the single serializer for a fit, so
    the catalog's on-disk format and a fit's own persistence story never
    diverge into two implementations.
    """
    return fit.to_dict()


def _fit_from_payload(payload: dict) -> SheddingFit:
    """Deserialize a fit. Delegates to ``SheddingFit.from_dict``."""
    return SheddingFit.from_dict(payload)


@dataclass
class SheddingCatalog:
    """
    A collection of fitted models with a browsable summary table.

    Example:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> catalog.table.shape
        (126, 24)
        >>> fit = catalog.select(
        ...     dataset_id='woelfel2020virological', analyte='stool', model='gamma'
        ... )
        >>> fit.dataset_id
        'woelfel2020virological'
    """

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
        """
        Build a ``SheddingEnsemble`` from the catalog's fits matching ``keys``.

        Args:
            dataset_ids: Optional list restricting components to named
                studies. A name with no matching fit raises rather than being
                dropped, which would silently shrink the ensemble.
            weights: ``"n_subjects"`` (default), ``"equal"``, or an explicit
                array of length equal to the number of matching fits.
            method: ``"mixture"`` (default), which draws each simulated
                individual from one contributing study's own fit, preserving
                between-study heterogeneity; or ``"moment"``, which collapses
                the components into a single Gaussian.
            **keys: Attribute filters, e.g. ``biomarker="SARS-CoV-2"``.

        Returns:
            A ``SheddingEnsemble``.

        Raises:
            ValueError: If the matching fits disagree on model, unit,
                reference event, biomarker, or specimen, or if one study
                contributes more than one analyte. See ``build_ensemble`` for
                the full validation this delegates to.
        """
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
    min_time: float | None = None,
    max_peak_above_observed: float | None = None,
) -> SheddingCatalog:
    """
    Fit every analyte of every dataset, for every requested model.

    Analytes that cannot be fitted are recorded in ``catalog.skipped`` with a
    reason rather than raising, so one unsuitable analyte does not abort a
    repository-wide build — and so a missing study reads as unsuitable rather
    than as a bug.

    Beyond the reasons ``fit_shedding_model`` itself raises, this applies
    ``require_estimable_population``, so a fit with too few subjects to support
    a between-subject covariance is recorded as
    ``too_few_subjects_for_population`` instead of being published. That check
    lives here rather than in the fitter because fitting a single subject
    directly is legitimate; publishing it as a population is not.

    Args:
        datasets: Dataset dictionaries from ``load_dataset``.
        models: Model names to fit. Defaults to both.
        min_observations: Passed through to the fitter.
        min_time: Passed through to the fitter as the earliest usable reading
            time, in days from the reference event. ``None`` keeps its default.
        max_peak_above_observed: Passed through to the fitter, which uses it to
            decide when a subject's implied peak is extrapolation rather than
            estimate. ``None`` keeps the fitter's default. Lower it to rebuild
            the catalog under a stricter reading and compare.

    Returns:
        A ``SheddingCatalog``.

    Example:
        Fitting every analyte of every dataset can take minutes, so this is
        not run here; the call shape is shown for reference.

        >>> import shedding_hub as sh
        >>> data = sh.load_dataset(
        ...     'woelfel2020virological', local='./data'
        ... )  # doctest: +SKIP
        >>> catalog = sh.fit_shedding_models([data])  # doctest: +SKIP
        >>> len(catalog.fits)  # doctest: +SKIP
        4
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
                        extra = {}
                        if max_peak_above_observed is not None:
                            extra["max_peak_above_observed"] = max_peak_above_observed
                        if min_time is not None:
                            extra["min_time"] = min_time
                        fit = fit_shedding_model(
                            dataset,
                            analyte=analyte,
                            model=model,
                            min_observations=min_observations,
                            **extra,
                        )
                    # Applied here rather than inside fit_shedding_model so that
                    # fitting one subject on purpose stays possible, while a fit
                    # too thin to describe a population never reaches the
                    # catalog. Raises SheddingDataError, so the existing handler
                    # records it with a reason like any other refusal.
                    require_estimable_population(fit)
                    fits.append(fit)
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
                    # Not a convergence failure: non-convergence never raises
                    # (fit_shedding_model returns the fit with converged=False
                    # and it is published normally). Reaching here means some
                    # other, unanticipated ValueError/LinAlgError escaped
                    # every named SheddingDataError reason above -- e.g. a
                    # malformed dataset -- so it is filed as a catch-all
                    # rather than misnamed after a cause that cannot be true.
                    skipped.append(
                        {
                            "dataset_id": dataset_id,
                            "analyte": analyte,
                            "model": model,
                            "reason": "unexpected_error",
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

    Example:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> len(catalog.fits)
        126
    """
    catalog_path = pathlib.Path(path) if path else CATALOG_PATH
    if not catalog_path.is_file():
        raise FileNotFoundError(
            f"No shedding catalog at {catalog_path}. Run `make catalog` to build it."
        )
    with catalog_path.open(encoding="utf-8") as stream:
        payload = yaml.safe_load(stream)
    return SheddingCatalog.from_dict(payload)
