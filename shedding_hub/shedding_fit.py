"""
Fit shedding-curve models to Shedding Hub datasets.

Fitting is done per analyte by joint maximum likelihood over every subject's
parameters plus one shared measurement-error standard deviation, using a
left-censored normal likelihood so that ``negative`` measurements contribute the
information they carry (that the concentration was below the limit) instead of
being discarded. Roughly 37% of measurements in the repository are ``negative``;
dropping them biases decay rates slow and inflates simulated late-phase shedding.
"""

import math
import warnings
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .shedding_models import LN10, PARAM_NAMES, validate_model

CENSORING_MARGIN = 0.01
NEGATIVE_VALUE = "negative"

# Fecal-strength / normalization indicators, not pathogens shed by infected
# people. They have no time-since-infection trajectory, so fitting a shedding
# curve to them is meaningless regardless of how much data is available.
NON_PATHOGEN_BIOMARKERS = frozenset({"crAssphage", "PMMoV", "mtDNA"})


class SheddingDataError(ValueError):
    """
    Raised when an analyte cannot be fitted.

    Attributes:
        reason: Machine-readable cause, one of ``ct_units``,
            ``non_pathogen_biomarker``, ``too_few_subjects``,
            ``no_positive_measurements``, ``unknown_analyte``,
            ``degenerate_fit``. The catalog builder records this so a missing
            study reads as unsuitable, not as a bug.

            All but ``degenerate_fit`` are raised by ``prepare_observations``,
            before any fitting happens. ``degenerate_fit`` is raised by
            ``fit_shedding_model`` only after the optimizer has run, because it
            describes the fit rather than the data.
    """

    def __init__(self, message: str, reason: str):
        super().__init__(message)
        self.reason = reason


@dataclass
class Observations:
    """Model-ready observations for a single analyte."""

    subject_index: np.ndarray
    times: np.ndarray
    values: np.ndarray
    censored: np.ndarray
    censoring_limit: float
    subject_ids: list = field(default_factory=list)
    n_subjects: int = 0
    n_excluded_subjects: int = 0
    n_dropped_measurements: int = 0


def _is_ct_unit(unit: Any) -> bool:
    if unit is None:
        return False
    lowered = str(unit).strip().lower()
    return "cycle" in lowered or lowered == "ct"


def _numeric_limit(value: Any) -> float | None:
    """Return a positive numeric limit, or None if unusable (e.g. ``unknown``)."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value) if value > 0 else None


def _resolve_censoring_limit(analyte_spec: dict, observed_log10: np.ndarray) -> float:
    """
    Resolve the log10 censoring limit.

    Prefers the declared limit of quantification, then the limit of detection.
    Falls back to just below the smallest observed positive value when neither is
    usable, or when the declared limit is not strictly below every observation —
    the likelihood is only coherent if censored points really do sit below the
    limit.

    Assumes ``observed_log10`` is non-empty; ``prepare_observations`` guarantees
    this by raising ``no_positive_measurements`` itself before ever calling here.
    """
    declared = None
    for key in ("limit_of_quantification", "limit_of_detection"):
        limit = _numeric_limit(analyte_spec.get(key))
        if limit is not None:
            declared = math.log10(limit)
            break

    smallest = float(observed_log10.min())
    if declared is None or declared >= smallest:
        fallback = smallest - CENSORING_MARGIN
        warnings.warn(
            "Falling back to a censoring limit of "
            f"{fallback:.4g} (log10) because the declared limit "
            f"({'none' if declared is None else f'{declared:.4g}'}) is missing or "
            f"not below the smallest observed value ({smallest:.4g}).",
            UserWarning,
            stacklevel=2,
        )
        return fallback
    return declared


def prepare_observations(
    dataset: dict,
    analyte: str,
    model: str,
    *,
    min_observations: int | None = None,
) -> Observations:
    """
    Extract model-ready observations for one analyte of one dataset.

    Only the analyte-level ``limit_of_quantification``/``limit_of_detection``
    are used to resolve the censoring limit; no analyte in the repository that
    is eligible for fitting declares a per-measurement limit (the one analyte
    that did, crAssphage, is a non-pathogen indicator and is rejected before
    the censoring limit is ever resolved).

    Args:
        dataset: Dataset dictionary from ``load_dataset``.
        analyte: Key into ``dataset["analytes"]``.
        model: ``"exponential"`` or ``"gamma"``.
        min_observations: Minimum usable measurements a subject must have to be
            retained. Defaults to the number of per-subject parameters (3 for
            gamma, 2 for exponential). ``sigma`` is shared across subjects, so a
            subject does not need residual degrees of freedom of its own.

    Returns:
        An ``Observations`` instance with subject indices renumbered contiguously
        from zero over the retained subjects.

    Raises:
        SheddingDataError: The analyte is unknown, uses cycle-threshold units,
            is a non-pathogen indicator biomarker, has no positive
            measurements, or leaves no subject with enough data — reasons
            ``unknown_analyte``, ``ct_units``, ``non_pathogen_biomarker``,
            ``no_positive_measurements`` and ``too_few_subjects``
            respectively. The remaining reason, ``degenerate_fit``, cannot
            arise here: it describes a fit that collapsed onto the parameter
            bounds and so is raised by ``fit_shedding_model`` after optimizing.
    """
    validate_model(model)
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")
    for key in ("analytes", "participants"):
        if key not in dataset:
            raise ValueError(f"Dataset missing required key: {key}")

    analytes = dataset["analytes"]
    if analyte not in analytes:
        raise SheddingDataError(
            f"Analyte {analyte!r} not in dataset; available: {sorted(analytes)}.",
            "unknown_analyte",
        )
    analyte_spec = analytes[analyte]

    if _is_ct_unit(analyte_spec.get("unit")):
        raise SheddingDataError(
            f"Analyte {analyte!r} is reported in {analyte_spec.get('unit')!r}. "
            "Cycle-threshold values are inversely related to concentration and "
            "already on a log scale, so neither shedding model applies. Select a "
            "concentration analyte instead.",
            "ct_units",
        )

    biomarker = analyte_spec.get("biomarker")
    if biomarker in NON_PATHOGEN_BIOMARKERS:
        raise SheddingDataError(
            f"Analyte {analyte!r} measures {biomarker!r}, a fecal-strength/"
            "normalization indicator rather than a pathogen shed by infected "
            "people. It has no time-since-infection trajectory, so neither "
            "shedding model applies. Select a pathogen analyte instead.",
            "non_pathogen_biomarker",
        )

    if min_observations is None:
        min_observations = len(PARAM_NAMES[model])

    per_subject: list[dict[str, list]] = []
    subject_ids: list = []
    n_dropped = 0

    for position, participant in enumerate(dataset["participants"]):
        times: list[float] = []
        values: list[float] = []
        censored: list[bool] = []
        for measurement in participant.get("measurements") or []:
            if measurement.get("analyte") != analyte:
                continue
            time = measurement.get("time")
            if not isinstance(time, (int, float)) or isinstance(time, bool):
                n_dropped += 1
                continue
            time = float(time)
            if model == "gamma" and time <= 0:
                n_dropped += 1
                continue
            value = measurement.get("value")
            if isinstance(value, str):
                if value == NEGATIVE_VALUE:
                    times.append(time)
                    values.append(np.nan)
                    censored.append(True)
                else:
                    n_dropped += 1
                continue
            if not isinstance(value, (int, float)) or value <= 0:
                n_dropped += 1
                continue
            times.append(time)
            values.append(math.log10(float(value)))
            censored.append(False)

        if times:
            per_subject.append(
                {
                    "times": times,
                    "values": values,
                    "censored": censored,
                }
            )
            subject_ids.append(participant.get("patient_id", position + 1))

    # Filter subject_ids and per_subject together, from a single zipped list,
    # so the retention predicate appears exactly once: applying it separately
    # to each list could silently desync subject_ids from the arrays below.
    retained_pairs = [
        (subject_id, subject)
        for subject_id, subject in zip(subject_ids, per_subject)
        if len(subject["times"]) >= min_observations
    ]
    n_excluded = len(per_subject) - len(retained_pairs)
    retained_ids = [subject_id for subject_id, _ in retained_pairs]
    retained = [subject for _, subject in retained_pairs]

    if n_excluded:
        warnings.warn(
            f"{n_excluded} subject(s) excluded from the {analyte!r} fit for having "
            f"fewer than {min_observations} usable measurements.",
            UserWarning,
            stacklevel=2,
        )
    if n_dropped:
        warnings.warn(
            f"{n_dropped} measurement(s) dropped from the {analyte!r} fit "
            "(qualitative result, unknown time, or a non-positive time under the "
            "gamma model).",
            UserWarning,
            stacklevel=2,
        )
    if not retained:
        raise SheddingDataError(
            f"No subject has at least {min_observations} usable measurements for "
            f"analyte {analyte!r}.",
            "too_few_subjects",
        )

    subject_index = np.concatenate(
        [np.full(len(s["times"]), i, dtype=int) for i, s in enumerate(retained)]
    )
    times_array = np.concatenate([np.asarray(s["times"], float) for s in retained])
    values_array = np.concatenate([np.asarray(s["values"], float) for s in retained])
    censored_array = np.concatenate([np.asarray(s["censored"], bool) for s in retained])

    observed = values_array[~censored_array]
    if observed.size == 0:
        raise SheddingDataError(
            f"Analyte {analyte!r} has no positive measurements to fit.",
            "no_positive_measurements",
        )

    censoring_limit = _resolve_censoring_limit(analyte_spec, observed)

    return Observations(
        subject_index=subject_index,
        times=times_array,
        values=values_array,
        censored=censored_array,
        censoring_limit=censoring_limit,
        subject_ids=retained_ids,
        n_subjects=len(retained),
        n_excluded_subjects=n_excluded,
        n_dropped_measurements=n_dropped,
    )


import pandas as pd
from scipy import optimize
from scipy.stats import norm

from .shedding_models import (
    half_life_days,
    log10_concentration,
    log10_concentration_pointwise,
    peak_day,
)

# The positivity floor for a natural-scale parameter. Nothing is ever
# *initialized* here — see _DEFAULT_A0 for why that was a bug — but it remains
# the reference point that _DEGENERATE_PARAM is judged against.
_MIN_PARAM = 1e-6
_THETA_BOUNDS = (-25.0, 25.0)
_LOG_SIGMA_BOUNDS = (-10.0, 5.0)

# Parameters are optimized as theta = log(param), so the chain rule gives
# dL/dtheta = param * dL/dparam: the gradient vanishes as a parameter
# approaches zero. Near-zero is therefore an absorbing state — a parameter
# started at _MIN_PARAM has a gradient of order 1e-5 against order 1e+1 for a
# healthy coordinate, and the optimizer can never pull it back. Initialization
# must consequently never place a parameter at or near the floor, so a
# least-squares coefficient that is non-positive, non-finite, or merely
# negligible falls back to one of these data-driven defaults instead.
_DEFAULT_A0 = math.log(2.0) / 7.0  # a one-week half-life
_DEFAULT_B0 = 1.0  # rises then falls, peaking at 1/a0 days

# A fitted parameter at or below this magnitude has collapsed rather than been
# estimated. This sits four orders of magnitude above _MIN_PARAM deliberately:
# the vanishing gradient above means a collapsing parameter stalls somewhere
# near zero rather than reaching the bound exactly, so testing proximity to the
# literal bound does not detect it. Every parameter is already physically
# meaningless well before this point — a0 <= 1e-2 is a half-life beyond 69
# days, longer than any shedding episode in the repository; c0 <= 1e-2 puts the
# log10 intercept below 0.005, orders of magnitude under any assay's detection
# floor; and b0 <= 1e-2 leaves the gamma curve with no rise phase at all.
_DEGENERATE_PARAM = 1e-2

# How close to the upper theta bound counts as pinned against it. The upper
# bound is approached from a region where the gradient is huge rather than
# vanishing, so unlike the floor it really is reached, and the tolerance can be
# tight.
_BOUND_TOLERANCE = 1e-6


def _require_positive_semidefinite(cov: np.ndarray, *, advice: str) -> np.ndarray:
    """
    Raise if ``cov`` is not (numerically) positive semi-definite.

    Shared by ``SheddingFit.sample_params`` and ``SheddingEnsemble``'s
    ``method="moment"`` path, which draws from its own moment-matched
    covariance rather than delegating to a component fit.

    Args:
        cov: Candidate covariance matrix.
        advice: Sentence appended to the error, tailored to how the caller
            produced ``cov`` (e.g. what to try instead).

    Returns:
        ``cov`` as a float array, for convenient chaining at the call site.

    Raises:
        ValueError: If the smallest eigenvalue is negative beyond numerical
            tolerance.
    """
    cov = np.asarray(cov, dtype=float)
    eigenvalues = np.linalg.eigvalsh(cov)
    # Tolerance scaled to the matrix, not exact zero: an all-zeros covariance
    # is legitimate (e.g. a single-subject fit, where no between-subject
    # variance can be estimated) and must still simulate, producing identical
    # individuals. Floating-point noise from np.cov can likewise leave a
    # near-singular but truly PSD matrix with an eigenvalue just below zero;
    # only a *meaningfully* negative eigenvalue indicates a real problem.
    tolerance = (
        max(np.abs(eigenvalues).max(), 1.0) * cov.shape[0] * np.finfo(float).eps * 100
    )
    if eigenvalues.min() < -tolerance:
        raise ValueError(
            "population_cov is not positive semi-definite (smallest eigenvalue "
            f"{eigenvalues.min():.3g}, below numerical tolerance). {advice}"
        )
    return cov


@dataclass
class SheddingFit:
    """
    A fitted shedding model for one analyte of one dataset.

    ``population_mean`` and ``population_cov`` describe the distribution of
    ``theta = log(params)`` across subjects; drawing from that multivariate normal
    and exponentiating produces a new plausible individual, which is what makes
    simulation possible.

    Two-stage estimation does not shrink individual estimates toward the
    population mean the way a hierarchical Bayesian fit does, so
    ``population_cov`` absorbs within-subject estimation error and overestimates
    true between-subject variance. Simulated cohorts are therefore somewhat more
    dispersed than reality, the more so when subjects have few observations.

    For the gamma model specifically, ``population_mean[1]`` (``b0``, the
    rise-rate/shape parameter) is additionally downward-biased at realistic
    sampling densities: roughly -0.15 log units at ~14 observations per
    subject, shrinking toward zero as sampling density rises. Because
    ``peak_day = b0 / a0``, this makes fitted peak-shedding timing somewhat
    early for sparsely-sampled studies. This is a property of two-stage maximum
    likelihood, not a bug, and is the main reason a hierarchical Bayesian
    backend would improve on these estimates.

    That figure was previously about -0.5 log units. Most of the difference was
    not two-stage bias at all but an initialization artifact: seeding the gamma
    fit from its own collinear ``[1, ln(t), -t]`` design produced negative
    coefficients that were clipped onto the parameter floor, which the
    optimizer could not escape. Measured over six seeds after that fix, the
    sparse-sampling bias averages 0.15 log units (range 0.02 to 0.27).
    """

    model: str
    method: str
    population_mean: np.ndarray
    population_cov: np.ndarray
    sigma: float
    subject_params: pd.DataFrame | None
    censoring_limit: float
    dataset_id: str
    analyte: str
    biomarker: str | None
    specimen: str | None
    reference_event: str | None
    unit: str | None
    gene_target: str | None
    dose: int | None
    vaccine_type: str | None
    n_subjects: int
    n_measurements: int
    n_censored: int
    n_excluded_subjects: int
    n_dropped_measurements: int
    converged: bool
    log_likelihood: float
    aic: float
    # Subjects whose fits collapsed onto the parameter bounds. They are present
    # in ``subject_params`` (flagged by its ``degenerate`` column) but excluded
    # from ``population_mean``/``population_cov``, so ``n_subjects`` exceeds the
    # number actually summarized by exactly this count. Defaults to zero so that
    # directly-constructed fits and catalogs written before this field existed
    # both remain loadable.
    n_degenerate_subjects: int = 0

    @property
    def param_names(self) -> tuple[str, ...]:
        return PARAM_NAMES[self.model]

    @property
    def median_params(self) -> np.ndarray:
        """
        Parameters of the median individual.

        Because ``theta`` is normal, the parameters are lognormal, so
        ``exp(population_mean)`` is exactly their median — not their mean. Note
        that the median individual's trajectory is not the population's mean
        trajectory; to aggregate load across a cohort, simulate rather than
        scaling this up.
        """
        return np.exp(self.population_mean)

    @property
    def peak_day(self) -> float:
        return float(peak_day(self.model, self.median_params[None, :])[0])

    @property
    def peak_log10(self) -> float:
        """Log10 concentration of the median individual at its peak."""
        if self.model == "exponential":
            return float(
                log10_concentration(
                    self.model, self.median_params[None, :], np.array([0.0])
                )[0, 0]
            )
        return float(
            log10_concentration(
                self.model, self.median_params[None, :], np.array([self.peak_day])
            )[0, 0]
        )

    @property
    def half_life_days(self) -> float:
        return float(half_life_days(self.model, self.median_params[None, :])[0])

    def sample_params(
        self, rng: np.random.Generator, n: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Draw ``n`` individuals' natural-scale parameters.

        Returns:
            ``(params, sources)`` where ``params`` has shape ``(n, k)`` and
            ``sources`` names the dataset each individual came from — constant
            here, but varying for a mixture ensemble, so both share this
            interface.
        """
        if n < 1:
            raise ValueError("n_individuals must be at least 1")
        cov = _require_positive_semidefinite(
            self.population_cov,
            advice=(
                "This usually means too few subjects survived fitting to "
                "estimate a stable between-subject covariance. Consider "
                "pooling multiple studies into a SheddingEnsemble instead of "
                "simulating from this fit alone."
            ),
        )
        theta = rng.multivariate_normal(self.population_mean, cov, n)
        return np.exp(theta), np.full(n, self.dataset_id, dtype=object)


def _initial_theta(model: str, observations: Observations) -> np.ndarray:
    """
    Initialize per-subject log-parameters by ordinary least squares.

    Fits each subject's uncensored points; subjects with too few uncensored
    points fall back to a pooled fit across all subjects, so an all-censored or
    nearly all-censored subject still starts somewhere sensible.

    Both models are seeded from the same two-column decay design ``[1, -t]``,
    even though the gamma model has three parameters. The gamma model's own
    design ``[1, ln(t), -t]`` is badly conditioned over a realistic sampling
    window — over ``t = 6..23``, ``ln(t)`` and ``t`` correlate at about 0.98 —
    so its least-squares solution is unstable and routinely returns negative
    coefficients, which previously started those parameters at the absorbing
    floor described on ``_DEFAULT_A0``. ``b0`` is instead seeded from a modest
    positive default, and ``a0``/``c0`` from the well-conditioned decay design.

    No coefficient is ever clipped to ``_MIN_PARAM``: one that is non-positive,
    non-finite, or negligible falls back to a data-driven default instead.
    """
    uncensored = ~observations.censored

    def solve(times: np.ndarray, values: np.ndarray) -> np.ndarray:
        design = np.column_stack([np.ones_like(times), -times])
        coefficients, *_ = np.linalg.lstsq(design, values * LN10, rcond=None)
        c0, a0 = coefficients
        if not np.isfinite(a0) or a0 <= _DEGENERATE_PARAM:
            a0 = _DEFAULT_A0
        if not np.isfinite(c0) or c0 <= _DEGENERATE_PARAM:
            # The largest value this subject actually reached, as a natural-log
            # intercept: a curve starting there is consistent with the data even
            # when the regression slope through it was not.
            c0 = LN10 * float(np.max(values))
        if model == "exponential":
            return np.log([a0, c0])
        return np.log([a0, _DEFAULT_B0, c0])

    pooled = solve(observations.times[uncensored], observations.values[uncensored])

    theta = np.tile(pooled, (observations.n_subjects, 1))
    for i in range(observations.n_subjects):
        mask = uncensored & (observations.subject_index == i)
        # Two points determine the decay design, whichever model is being fitted
        # — the gamma model no longer needs three to seed itself.
        if mask.sum() >= 2:
            try:
                theta[i] = solve(observations.times[mask], observations.values[mask])
            except np.linalg.LinAlgError:
                pass
    return np.clip(theta, *_THETA_BOUNDS)


def _degenerate_subjects(theta: np.ndarray) -> np.ndarray:
    """
    Flag subjects whose fit collapsed onto the parameter bounds.

    A collapsed subject's ``theta`` is not an estimate but a boundary solution,
    and averaging it into ``mean(theta_i)`` distorts the population summary out
    of all proportion: one subject pinned at ``log(1e-6) = -13.8`` is enough to
    turn a one-day half-life into a 278-day one.

    Args:
        theta: Fitted log-parameters, shape ``(n_subjects, k)``.

    Returns:
        Boolean array of length ``n_subjects``, True where any coordinate has
        collapsed toward zero or run away to the upper bound.
    """
    collapsed = theta <= math.log(_DEGENERATE_PARAM)
    runaway = theta >= _THETA_BOUNDS[1] - _BOUND_TOLERANCE
    return np.asarray((collapsed | runaway).any(axis=1))


def _negative_log_likelihood(
    x: np.ndarray, model: str, observations: Observations
) -> float:
    """
    Negative log likelihood with left-censored observations.

    Uncensored points contribute a normal density; censored points contribute
    ``Phi((L - mu) / sigma)``, the probability of falling below the limit. This is
    the direct analogue of Stan's ``normal_lcdf`` term.
    """
    k = len(PARAM_NAMES[model])
    n = observations.n_subjects
    theta = x[: n * k].reshape(n, k)
    log_sigma = x[-1]
    sigma = math.exp(log_sigma)

    params = np.exp(theta)[observations.subject_index]
    predicted = log10_concentration_pointwise(model, params, observations.times)
    if not np.all(np.isfinite(predicted)):
        return np.inf

    total = 0.0
    uncensored = ~observations.censored
    if uncensored.any():
        residual = (observations.values[uncensored] - predicted[uncensored]) / sigma
        total += 0.5 * float(np.sum(residual**2))
        total += float(uncensored.sum()) * (log_sigma + 0.5 * math.log(2 * math.pi))
    if observations.censored.any():
        z = (observations.censoring_limit - predicted[observations.censored]) / sigma
        total -= float(np.sum(norm.logcdf(z)))
    return total if np.isfinite(total) else np.inf


def fit_shedding_model(
    dataset: dict,
    *,
    analyte: str,
    model: str = "gamma",
    min_observations: int | None = None,
) -> SheddingFit:
    """
    Fit a shedding model to one analyte by censored maximum likelihood.

    Args:
        dataset: Dataset dictionary from ``load_dataset``.
        analyte: Key into ``dataset["analytes"]``.
        model: ``"exponential"`` or ``"gamma"``.
        min_observations: Minimum usable measurements per subject; defaults to the
            number of per-subject parameters.

    Returns:
        A ``SheddingFit``. Subjects whose parameters collapsed onto the bounds
        are excluded from the population summary but retained in
        ``subject_params`` with ``degenerate`` set; ``n_degenerate_subjects``
        counts them.

    Raises:
        SheddingDataError: The analyte cannot be fitted (see ``reason``). In
            addition to every reason ``prepare_observations`` can raise, this
            adds ``degenerate_fit``: too many subjects' fits collapsed onto the
            parameter bounds for a population covariance to be estimable.

    Note:
        ``aic`` is only comparable between models fitted to the same
        observations. The gamma model drops non-positive times while the
        exponential model keeps them, so compare ``n_measurements`` before
        comparing ``aic`` across models.
    """
    validate_model(model)
    observations = prepare_observations(
        dataset, analyte, model, min_observations=min_observations
    )

    k = len(PARAM_NAMES[model])
    n = observations.n_subjects
    n_parameters = n * k + 1  # every subject's k parameters, plus one shared sigma
    x0 = np.concatenate([_initial_theta(model, observations).ravel(), [math.log(0.5)]])
    bounds = [_THETA_BOUNDS] * (n * k) + [_LOG_SIGMA_BOUNDS]

    # scipy's L-BFGS-B defaults (maxfun=maxiter=15000, ftol=2.22e-9) are tuned
    # for small problems. Jointly optimizing every subject's parameters plus one
    # shared sigma is n_parameters-dimensional, and heavy censoring flattens the
    # likelihood near its optimum, so the default budget/tolerance combination
    # can report spurious non-convergence tens of thousands of evaluations
    # before the fit stops moving in any way that matters (verified against a
    # fully-converged reference run: the recovered parameters agree to 3+
    # decimal places well before the default ftol would be satisfied). Scaling
    # the budget with problem size and relaxing ftol lets these fits report
    # ``converged=True`` honestly instead of exhausting the default cap.
    #
    # The multiplier is 1000 rather than the 500 that sufficed while
    # initialization was collapsing subjects onto the parameter floor. That
    # collapse was flattering these counts: a floor-pinned subject has a
    # vanishing gradient, so the optimizer stopped working on it almost
    # immediately and declared success early on a fit that was simply stuck.
    # With every subject genuinely being optimized, a 60-subject gamma fit
    # needs about 92k evaluations where 500 * n_parameters allowed only 90.5k —
    # it was failing by 1.5%, and the extra work buys a converged fit with
    # identical parameters (verified: nll 644.697 vs 644.709, mean log b0
    # identical to four decimal places).
    max_evaluations = max(15000, 1000 * n_parameters)
    result = optimize.minimize(
        _negative_log_likelihood,
        x0,
        args=(model, observations),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxfun": max_evaluations, "maxiter": max_evaluations, "ftol": 1e-6},
    )
    if not result.success:
        warnings.warn(
            f"Optimizer did not converge for analyte {analyte!r} "
            f"({result.message}). The fit is returned with converged=False.",
            UserWarning,
            stacklevel=2,
        )

    theta = result.x[: n * k].reshape(n, k)
    sigma = float(np.exp(result.x[-1]))

    # Subjects whose fits collapsed onto the bounds stay in subject_params so
    # the raw fit remains inspectable, but are kept out of the population
    # summary, which they would otherwise dominate.
    degenerate = _degenerate_subjects(theta)
    n_degenerate = int(degenerate.sum())
    retained = ~degenerate
    n_retained = int(retained.sum())
    # A single-subject fit legitimately yields a zero covariance and is allowed;
    # what is not allowed is *degeneracy* leaving too little behind. Hence
    # min(2, n) rather than a flat 2: with n == 1 the bar is that the one
    # subject survived, and only with n >= 2 does Sigma become estimable at all.
    if n_retained < min(2, n):
        raise SheddingDataError(
            f"{n_degenerate} of {n} subject(s) for analyte {analyte!r} collapsed to "
            f"the parameter bounds under the {model!r} model, leaving {n_retained} "
            "usable subject(s) — too few to estimate a population covariance. This "
            "usually means the model is not identifiable from this data: for the "
            "gamma model, typically because sampling began after peak shedding, so "
            "there is no rise phase from which to estimate b0.",
            "degenerate_fit",
        )
    if n_degenerate:
        warnings.warn(
            f"{n_degenerate} subject(s) excluded from the population summary of "
            f"analyte {analyte!r}: their fitted parameters collapsed onto the "
            "bounds. They remain in subject_params, flagged by the 'degenerate' "
            "column.",
            UserWarning,
            stacklevel=2,
        )

    kept = theta[retained]
    population_mean = kept.mean(axis=0)
    population_cov = (
        np.cov(kept, rowvar=False, ddof=1) if n_retained > 1 else np.zeros((k, k))
    )
    population_cov = np.atleast_2d(population_cov)

    subject_params = pd.DataFrame(np.exp(theta), columns=list(PARAM_NAMES[model]))
    subject_params.insert(0, "subject_id", observations.subject_ids)
    subject_params["degenerate"] = degenerate

    log_likelihood = -float(result.fun)
    analyte_spec = dataset["analytes"][analyte]
    specimen = analyte_spec.get("specimen")
    if isinstance(specimen, list):
        specimen = "+".join(specimen)

    return SheddingFit(
        model=model,
        method="mle",
        population_mean=population_mean,
        population_cov=population_cov,
        sigma=sigma,
        subject_params=subject_params,
        censoring_limit=observations.censoring_limit,
        dataset_id=dataset.get("dataset_id", "unknown"),
        analyte=analyte,
        biomarker=analyte_spec.get("biomarker"),
        specimen=specimen,
        reference_event=analyte_spec.get("reference_event"),
        unit=analyte_spec.get("unit"),
        gene_target=analyte_spec.get("gene_target"),
        dose=analyte_spec.get("dose"),
        vaccine_type=analyte_spec.get("vaccine_type"),
        n_subjects=n,
        n_measurements=int(observations.times.size),
        n_censored=int(observations.censored.sum()),
        n_excluded_subjects=observations.n_excluded_subjects,
        n_dropped_measurements=observations.n_dropped_measurements,
        converged=bool(result.success),
        log_likelihood=log_likelihood,
        aic=2.0 * n_parameters - 2.0 * log_likelihood,
        n_degenerate_subjects=n_degenerate,
    )
