"""
Combine per-study shedding fits into a cross-study ensemble.

Two methods are available. ``mixture`` draws a study per simulated individual and
then draws that individual's parameters from the study's fit, preserving
between-study heterogeneity and keeping the distribution multimodal when studies
genuinely disagree. ``moment`` collapses the components into a single Gaussian
whose covariance is within-study plus between-study variance.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .shedding_catalog import SheddingCatalog, _fits_to_frame
from .shedding_fit import (
    SheddingFit,
    _scaled,
    validate_dispersion,
    _require_positive_semidefinite,
)
from .shedding_models import from_population_coords

_COMPATIBILITY_KEYS = (
    "model",
    "unit",
    "reference_event",
    "biomarker",
    "specimen",
)


@dataclass
class SheddingEnsemble:
    """
    An ensemble of per-study fits sharing a biomarker, specimen, and unit.

    Example:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> ensemble = sh.shedding_for('SARS-CoV-2', 'stool', catalog=catalog)
        >>> ensemble.model
        'gamma'
        >>> ensemble.components.shape
        (2, 25)
        >>> ensemble.median_params
        Traceback (most recent call last):
            ...
        ValueError: median_params is only defined for method='moment'. This is a
            mixture ensemble, which has no closed-form median: simulate and
            take empirical quantiles instead.
    """

    fits: list[SheddingFit]
    weights: np.ndarray
    method: str
    # Provenance, not statistics: set by shedding_for to record why this
    # combination was chosen over the alternatives. Typed loosely to keep
    # shedding_ensemble free of an import from shedding_select, which imports it.
    # Deliberately absent from to_dict: it describes a choice made against one
    # catalog, and would be misleading if restored beside fits from another.
    selection: object = None

    @property
    def components(self) -> pd.DataFrame:
        """One row per contributing fit, with the same columns as the catalog."""
        frame = _fits_to_frame(self.fits)
        frame.insert(0, "weight", self.weights)
        return frame

    @property
    def model(self) -> str:
        return self.fits[0].model

    @property
    def sigma(self) -> float:
        """
        Effective measurement-error SD across components.

        An individual is drawn from component ``s`` with probability
        ``weights[s]``, so the marginal measurement-error variance is the
        weighted mean of the component *variances*, ``sum_s w_s * sigma_s**2``
        (variances add across a mixture; standard deviations do not). This is
        therefore the root-mean-square of the component sigmas weighted by
        ``weights``, not their weighted arithmetic mean, which would
        understate the effective noise.
        """
        sigmas = np.array([fit.sigma for fit in self.fits])
        return float(np.sqrt(np.sum(self.weights * sigmas**2)))

    @property
    def censoring_limit(self) -> float:
        """The most conservative (highest) limit among the components."""
        return float(max(fit.censoring_limit for fit in self.fits))

    @property
    def reference_event(self) -> str | None:
        return self.fits[0].reference_event

    @property
    def unit(self) -> str | None:
        return self.fits[0].unit

    @property
    def biomarker(self) -> str | None:
        return self.fits[0].biomarker

    @property
    def specimen(self) -> str | None:
        return self.fits[0].specimen

    @property
    def population_mean(self) -> np.ndarray:
        """Weighted mean of the component means."""
        means = np.array([fit.population_mean for fit in self.fits])
        return np.sum(self.weights[:, None] * means, axis=0)

    @property
    def population_cov(self) -> np.ndarray:
        """
        Moment-matched covariance: within-study plus between-study.

        Defined for both methods, but only used for sampling under
        ``method="moment"``.
        """
        means = np.array([fit.population_mean for fit in self.fits])
        covs = np.array([fit.population_cov for fit in self.fits])
        within = np.sum(self.weights[:, None, None] * covs, axis=0)
        deviations = means - self.population_mean
        between = np.einsum("s,si,sj->ij", self.weights, deviations, deviations)
        return within + between

    @property
    def median_params(self) -> np.ndarray:
        """
        Parameters of the median individual.

        Only well defined for ``method="moment"``. A mixture has no closed-form
        median, so rather than report a misleading number, simulate with
        ``simulate_shedding`` and take empirical quantiles of the result — which
        is the quantity you actually want.

        Mapped back through ``from_population_coords`` rather than exponentiated:
        the components share a model, hence a coordinate space, and for the gamma
        model that space is not the log-parameters.
        """
        if self.method != "moment":
            raise ValueError(
                "median_params is only defined for method='moment'. This is a "
                "mixture ensemble, which has no closed-form median: simulate and "
                "take empirical quantiles instead."
            )
        return from_population_coords(self.model, self.population_mean[None, :])[0]

    def sample_params(
        self, rng: np.random.Generator, n: int, dispersion: float = 1.0
    ) -> tuple[np.ndarray, np.ndarray]:
        """Draw ``n`` individuals' natural-scale parameters.

        ``dispersion`` scales each covariance drawn from by ``dispersion ** 2``;
        under ``method="mixture"`` that shrinks each component around its own
        mean, leaving the between-study spread of those means intact.
        """
        if n < 1:
            raise ValueError("n_individuals must be at least 1")
        validate_dispersion(dispersion)
        if self.method == "moment":
            cov = _require_positive_semidefinite(
                self.population_cov,
                advice=(
                    "The moment-matched covariance is within-study plus "
                    "between-study variance computed from this ensemble's own "
                    "components; this usually means the components' means are "
                    "nearly collinear. Consider method='mixture' instead, which "
                    "draws each individual from one component's own (already "
                    "validated) covariance rather than a combined one."
                ),
            )
            theta = rng.multivariate_normal(
                self.population_mean, _scaled(cov, dispersion), n
            )
            return (
                from_population_coords(self.model, theta),
                np.full(n, "ensemble", dtype=object),
            )

        if len(self.fits) == 1:
            # Skip the categorical draw entirely so a one-study ensemble consumes
            # the generator exactly as the underlying fit would, making the two
            # interchangeable for a given seed.
            return self.fits[0].sample_params(rng, n, dispersion)

        choices = rng.choice(len(self.fits), size=n, p=self.weights)
        k = self.fits[0].population_mean.size
        theta = np.empty((n, k))
        sources = np.empty(n, dtype=object)
        for index, fit in enumerate(self.fits):
            mask = choices == index
            count = int(mask.sum())
            if not count:
                continue
            # Validated here, per component, rather than inherited for free the
            # way the single-component shortcut above gets it from
            # fit.sample_params: this loop calls rng.multivariate_normal
            # directly, so without this check a component with a bad
            # covariance would slip through silently. This is also what makes
            # the moment path's advice to "consider method='mixture' instead,
            # which uses ... (already validated) covariances" actually true.
            cov = _require_positive_semidefinite(
                fit.population_cov,
                advice=(
                    f"Component {fit.dataset_id!r} of this mixture ensemble "
                    "has an invalid population covariance; drop it or refit "
                    "it before building the ensemble."
                ),
            )
            theta[mask] = rng.multivariate_normal(
                fit.population_mean, _scaled(cov, dispersion), count
            )
            sources[mask] = fit.dataset_id
        # One conversion for the whole array: make_ensemble guarantees every
        # component shares a model, so every row is in the same coordinates.
        return from_population_coords(self.model, theta), sources

    def to_dict(self) -> dict:
        """
        Serialize this ensemble to a JSON/YAML-safe dict.

        Each component fit is serialized via ``SheddingFit.to_dict`` (so
        ``subject_params`` is omitted there too, for the same reason),
        alongside the resolved weights and the combination method.

        Returns:
            A dict with keys ``fits``, ``weights``, and ``method``.
        """
        return {
            "fits": [fit.to_dict() for fit in self.fits],
            "weights": [float(w) for w in self.weights],
            "method": self.method,
        }

    @classmethod
    def from_dict(cls, payload: dict) -> "SheddingEnsemble":
        """
        Reconstruct an ensemble from ``to_dict``'s output.

        Component fits come back with ``subject_params is None``, same as any
        fit loaded from a catalog. Weights are restored exactly as saved
        rather than re-resolved through ``make_ensemble``, so a round trip
        cannot perturb them.

        Args:
            payload: A dict produced by ``to_dict``.

        Returns:
            A ``SheddingEnsemble``.
        """
        return cls(
            fits=[SheddingFit.from_dict(item) for item in payload["fits"]],
            weights=np.asarray(payload["weights"], dtype=float),
            method=payload["method"],
        )


def _resolve_weights(fits: list[SheddingFit], weights) -> np.ndarray:
    # Dispatch on isinstance(weights, str) first, rather than comparing
    # `weights == "equal"` directly: once weights is a numpy array (the
    # "explicit array" the docstring documents), `weights == "equal"` is an
    # elementwise comparison that raises ValueError on any array with more
    # than one element, before the length/negativity validation below ever
    # runs.
    if isinstance(weights, str):
        if weights == "equal":
            raw = np.ones(len(fits))
        elif weights == "n_subjects":
            raw = np.array([fit.n_subjects for fit in fits], dtype=float)
        else:
            raise ValueError(
                f"weights must be 'n_subjects', 'equal', or an array of length "
                f"{len(fits)}; got {weights!r}."
            )
    else:
        raw = np.asarray(weights, dtype=float)
        if raw.shape != (len(fits),):
            raise ValueError(
                f"weights must be 'n_subjects', 'equal', or an array of length "
                f"{len(fits)}; got shape {raw.shape}."
            )
    negative = raw[raw < 0]
    if negative.size:
        raise ValueError(
            f"Ensemble weights must be non-negative; got negative value(s) "
            f"{negative.tolist()}."
        )
    if raw.sum() <= 0:
        raise ValueError("Ensemble weights must sum to a positive value.")
    return raw / raw.sum()


def make_ensemble(
    fits: list[SheddingFit],
    *,
    weights="n_subjects",
    method: str = "mixture",
) -> SheddingEnsemble:
    """
    Assemble an ensemble from explicit fits.

    Accepts fits from anywhere — the shipped catalog, a fresh fit on private
    data, or a mixture of both.

    Args:
        fits: Component fits. A single-component ensemble is legal and behaves
            identically to the underlying fit, so callers can keep one code path
            regardless of how many studies they selected.
        weights: ``"n_subjects"`` (default), ``"equal"``, or an explicit array.
        method: ``"mixture"`` (default) or ``"moment"``.

    Returns:
        A ``SheddingEnsemble``.

    Raises:
        ValueError: If the fits disagree on model, unit, reference event,
            biomarker, or specimen, or if one study contributes more than one
            analyte.

    Example:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> fits = [
        ...     fit for fit in catalog.fits
        ...     if fit.dataset_id in ('lui2020viral', 'woelfel2020virological')
        ...     and fit.model == 'gamma' and fit.specimen == 'stool'
        ... ]
        >>> ensemble = sh.make_ensemble(fits)
        >>> ensemble.method
        'mixture'
        >>> len(ensemble.fits)
        2
    """
    fits = list(fits)
    if not fits:
        raise ValueError("An ensemble needs at least one fit.")
    if method not in ("mixture", "moment"):
        raise ValueError(
            f"Unknown ensemble method {method!r}; choose 'mixture' or 'moment'."
        )

    for key in _COMPATIBILITY_KEYS:
        values = {getattr(fit, key) for fit in fits}
        if len(values) > 1:
            raise ValueError(
                f"Ensemble components disagree on {key}: {sorted(map(str, values))}. "
                "Estimates are only comparable within one of these."
            )

    dataset_ids = [fit.dataset_id for fit in fits]
    duplicated = {name for name in dataset_ids if dataset_ids.count(name) > 1}
    if duplicated:
        offenders = _fits_to_frame(
            [fit for fit in fits if fit.dataset_id in duplicated]
        )
        raise ValueError(
            f"Study/studies {sorted(duplicated)} contribute more than one analyte "
            "to this ensemble, which would enter their subjects twice. Narrow the "
            "selection (for example by gene_target or analyte). Candidates:\n"
            f"{offenders[['dataset_id', 'analyte', 'gene_target']].to_string(index=False)}"
        )

    return SheddingEnsemble(
        fits=fits, weights=_resolve_weights(fits, weights), method=method
    )


def build_ensemble(
    catalog: SheddingCatalog,
    *,
    dataset_ids=None,
    weights="n_subjects",
    method: str = "mixture",
    **keys,
) -> SheddingEnsemble:
    """
    Build an ensemble from the catalog fits matching ``keys``.

    Args:
        catalog: The catalog to draw components from.
        dataset_ids: Optional list restricting components to named studies. A
            name with no matching fit raises rather than being dropped, which
            would silently shrink the ensemble.
        weights: Passed to ``make_ensemble``.
        method: Passed to ``make_ensemble``.
        **keys: Attribute filters, e.g. ``biomarker="SARS-CoV-2"``.

    Returns:
        A ``SheddingEnsemble``.
    """
    matches = [
        fit
        for fit in catalog.fits
        if all(getattr(fit, key, None) == value for key, value in keys.items())
    ]
    if dataset_ids is not None:
        requested = list(dataset_ids)
        available = {fit.dataset_id for fit in matches}
        missing = [name for name in requested if name not in available]
        if missing:
            raise ValueError(
                f"No fit matching {keys} for dataset_id(s) {missing}. "
                f"Matching studies are {sorted(available)}."
            )
        matches = [fit for fit in matches if fit.dataset_id in set(requested)]
    if not matches:
        raise ValueError(
            f"No fits match {keys}. Browse `catalog.table` for available "
            "combinations."
        )
    return make_ensemble(matches, weights=weights, method=method)
