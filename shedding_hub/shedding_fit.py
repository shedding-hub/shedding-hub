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
            ``no_rise_observed``, ``degenerate_fit``,
            ``too_few_subjects_for_population``. The catalog builder records
            this so a missing study reads as unsuitable, not as a bug.

            Raised in three places, by increasing lateness:

            - ``prepare_observations`` raises the first five, before any fitting
              happens. Note that ``too_few_subjects`` there means no *subject*
              had enough observations — a different failure from
              ``too_few_subjects_for_population`` below, which means there were
              not enough subjects.
            - ``fit_shedding_model`` raises ``no_rise_observed`` once the
              observations are in hand but before optimizing, and
              ``degenerate_fit`` only afterwards, because it describes the fit
              rather than the data.
            - ``require_estimable_population`` raises
              ``too_few_subjects_for_population``. It is applied by the catalog
              builder rather than by ``fit_shedding_model``, so that fitting a
              single subject on purpose stays possible; see its docstring.
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

    Uses the declared limit of quantification, then the limit of detection,
    whenever one is given, and always uses it as-is. Only ``negative``
    measurements are censored — at this limit — while every reported positive is
    kept as observed data, including any that fall below the limit. A number the
    assay reported below its limit of quantification is still a measurement
    (detected, if less precisely), so it is used rather than discarded; the limit
    describes only the value a ``negative`` is known to lie below. The
    observed-value likelihood term never references the limit, so keeping a
    positive below it is well defined.

    Falls back to just below the smallest observed positive only when neither
    limit is declared, so that any ``negative`` still sits below the resolved
    limit.

    Assumes ``observed_log10`` is non-empty; ``prepare_observations`` guarantees
    this by raising ``no_positive_measurements`` itself before ever calling here.
    """
    for key in ("limit_of_quantification", "limit_of_detection"):
        limit = _numeric_limit(analyte_spec.get(key))
        if limit is not None:
            return math.log10(limit)

    smallest = float(observed_log10.min())
    fallback = smallest - CENSORING_MARGIN
    warnings.warn(
        "Falling back to a censoring limit of "
        f"{fallback:.4g} (log10) because no limit of quantification or detection "
        "is declared for this analyte.",
        UserWarning,
        stacklevel=2,
    )
    return fallback


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
            respectively. Note that ``too_few_subjects`` here means no subject
            cleared ``min_observations``, not that the analyte has few
            subjects. The remaining three cannot arise here:
            ``no_rise_observed`` and ``degenerate_fit`` belong to
            ``fit_shedding_model``, and ``too_few_subjects_for_population`` to
            ``require_estimable_population``.
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
    log10_concentration_pointwise,
    log10_concentration_rowwise,
    peak_day,
)

# The historical positivity floor for a natural-scale parameter is 1e-6.
# Nothing is ever *initialized* there — see _DEFAULT_A0 for why that was a
# bug — but it remains the reference point that _DEGENERATE_PARAM below is
# judged against, and is not otherwise given a named constant since nothing
# else in this module reads it directly.
_THETA_BOUNDS = (-25.0, 25.0)
_LOG_SIGMA_BOUNDS = (-10.0, 5.0)

# Parameters are optimized as theta = log(param), so the chain rule gives
# dL/dtheta = param * dL/dparam: the gradient vanishes as a parameter
# approaches zero. Near-zero is therefore an absorbing state — a parameter
# started at the 1e-6 floor above has a gradient of order 1e-5 against order
# 1e+1 for a healthy coordinate, and the optimizer can never pull it back.
# Initialization must consequently never place a parameter at or near the
# floor, so a least-squares coefficient that is non-positive, non-finite, or
# merely negligible falls back to one of these data-driven defaults instead.
_DEFAULT_A0 = math.log(2.0) / 7.0  # a one-week half-life
_DEFAULT_B0 = 1.0  # rises then falls, peaking at 1/a0 days

# A fitted parameter at or below this magnitude has collapsed rather than been
# estimated. This sits four orders of magnitude above the 1e-6 floor
# deliberately, and is a magnitude test rather than a proximity-to-the-bound
# test, because the vanishing gradient described above means a collapsing
# parameter *stalls* near zero instead of ever reaching the bound. Confirmed
# three ways on the fit that motivated it (woelfel stool gamma, subject 3,
# whose c0 is collapsing): with the theta bound at -25 its c0 settles at
# 2.97e-3; moving the bound up to log(1e-6) = -13.82 moves it only to 8.39e-3;
# and restarting L-BFGS-B twelve times from its own output leaves it there,
# improving the likelihood by ~1e-6 per restart. It is a boundary MLE, not
# under-convergence, and no tolerance tight enough to mean "at the bound"
# would ever fire.
#
# The level is set from repository-wide evidence rather than from that one
# subject. Over 7,739 per-subject fits, the number flagged is nearly flat from
# 1e-5 to 1e-2 (3,035 -> 3,613, +19%) and then climbs steeply (4,915 by 1e-1,
# +36% more). 1e-2 is the last point before the threshold starts consuming
# genuinely-estimated parameters.
#
# Each coordinate is independently meaningless by then: a0 <= 1e-2 is a
# half-life beyond 69 days, longer than any shedding episode in the repository;
# c0 <= 1e-2 puts the log10 intercept below 0.005, orders of magnitude under any
# assay's detection floor; and b0 <= 1e-2 leaves the gamma curve with no rise
# phase at all.
_DEGENERATE_PARAM = 1e-2

# How close to the upper theta bound counts as pinned against it. The upper
# bound is approached from a region where the gradient is huge rather than
# vanishing, so unlike the floor it really is reached, and the tolerance can be
# tight.
_BOUND_TOLERANCE = 1e-6

# The runaway counterpart to _DEGENERATE_PARAM, expressed on the derived
# half-life because that is the quantity with physical meaning: a fitted decay
# implying a shedding half-life under 2.4 hours cannot have been estimated from
# data sampled roughly daily. Nothing in the observations distinguishes "fell
# below detection overnight" from "fell a thousandfold overnight", so the
# optimizer is free to run the decay rate up without penalty, and it does — the
# repository contained subjects with a0 of 84 to 142, i.e. half-lives of five to
# eight minutes.
#
# Sited like _DEGENERATE_PARAM, from the repository rather than from the cases
# that prompted it. Half-life alone has no gap to cut in — its distribution is
# smooth, peaking near 0.7 days — but the *pathology* has a sharp edge, because
# a runaway decay drags c0 up with it: extrapolating an implausibly steep slope
# back to t = 0 inflates the intercept. Sweeping the threshold and asking what
# share of subjects in each half-life band carry an implausible c0 (> 100):
#
#     half-life band     < 0.01  0.01-0.05  0.05-0.1  0.1-0.2  0.2-0.5  > 0.5
#     share with c0>100   59.4%      24.2%      6.8%    0.00%    0.00%  0.00%
#
# The pathology vanishes exactly at 0.1 days and never reappears at any slower
# rate. Cutting there takes the worst surviving c0 from 7792.8 to 75.5 and
# leaves no subject above 100; cutting at 0.05 leaves five; cutting at 0.2 buys
# no further improvement and costs 213 more subjects. 0.1 is the knee, and it
# coincides with the independent physical argument above.
#
# It also removes the need for a separate c0 or b0 bound: after this cut the
# worst c0 is 75.5 and the worst b0 is 33.9, both on smooth tails with no gap to
# justify cutting into.
_MIN_HALF_LIFE_DAYS = 0.1

# The gamma model is only fitted where a rise is actually observed. Its ``b0``
# describes the rise to peak, so a study that sampled purely after peak shedding
# carries no information about it: the profile likelihood is then monotone in
# ``b0`` until the ``c0 > 0`` constraint binds, and the MLE is a boundary
# solution no initialization can avoid. A subject counts as observing a rise if
# its highest reading came later than its first, judged only on subjects with
# enough readings for "later" to mean anything.
_MIN_RISE_OBSERVATIONS = 3
_MIN_RISE_FRACTION = 0.5

# Draws used to summarize peak_log10 by its population median. Large enough that
# the Monte Carlo error in the median is well under the third decimal place
# (~1.25 * sd / sqrt(n)), and paired with a fixed seed so a catalog rebuild
# reproduces exactly.
_PEAK_LOG10_DRAWS = 10000
_PEAK_LOG10_SEED = 8601


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
    # Percentage (0-100, matching ``pct_censored``) of adequately-sampled
    # subjects whose highest reading came later than their first. The gamma
    # model is refused below 50%; on exponential fits this is informational, and
    # a low value there is expected rather than alarming — post-peak sampling is
    # exactly what the exponential model is for. NaN when no subject had enough
    # readings to judge.
    pct_subjects_with_rise: float = float("nan")
    # Median, across the subjects that feed the population summary, of each
    # subject's own earliest observation time. Everything the fit says about
    # times before this is extrapolation for most of its subjects: see the
    # warning on ``peak_log10``, which is evaluated at t = 0 for the exponential
    # model however late sampling actually began.
    median_first_observed_day: float = float("nan")

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
        """
        Population median of the log10 concentration at peak.

        Unlike ``peak_day`` and ``half_life_days``, this is estimated by drawing
        from ``MVN(population_mean, population_cov)`` rather than evaluated at
        ``median_params``, because it is the one summary for which those two
        differ.

        ``peak_day = b0/a0 = exp(theta_b - theta_a)`` and
        ``half_life_days = ln(2)/a0`` are monotone transforms of a single
        lognormal quantity, so the value at ``exp(mu)`` *is* their median
        exactly, and no sampling is needed. At the peak, though, ``a0*t = b0``,
        making this quantity ``(c0 + b0*(ln b0 - ln a0) - b0) / ln(10)`` — a
        nonlinear function of all three parameters at once, whose median is not
        the function evaluated at the median parameters.

        The difference is not academic. Evaluating at ``exp(mu)`` reports a
        peak below almost every individual subject's own peak whenever the
        per-subject parameters lie on a correlated ridge, which they do for the
        gamma model: ``b0`` and ``c0`` trade off against each other, so
        averaging their logs coordinate-wise lands on a parameter vector no
        subject actually has.

        For the exponential model the peak is at ``t = 0``, giving
        ``c0 / ln(10)`` — again a monotone transform of one lognormal, so
        sampling is unnecessary. It is done anyway so that both models take one
        code path; the two agree to Monte Carlo error, which is what
        ``_PEAK_LOG10_DRAWS`` is sized to keep negligible.

        .. warning::
            Read this together with ``median_first_observed_day``. The
            exponential model's peak is at ``t = 0`` — the reference event — by
            definition, but many studies only begin sampling well after it. When
            ``median_first_observed_day`` is well above zero, ``peak_log10`` is a
            backward extrapolation to a time most subjects were never observed
            at, not a measured concentration, and it grows roughly as
            ``a0 * median_first_observed_day / ln(10)`` log units beyond the last
            value the study actually saw. A large
            ``median_first_observed_day`` together with a large ``peak_log10``
            means precisely that and should not be read as the study having
            detected ``10 ** peak_log10`` of anything.

            The definition is deliberately left at ``t = 0`` regardless, because
            that is what makes the value comparable across studies that started
            sampling at different times. The extrapolation is surfaced rather
            than hidden.
        """
        rng = np.random.default_rng(_PEAK_LOG10_SEED)
        # check_valid="ignore": an all-zero covariance (a single-subject fit) is
        # legitimate here and simply yields identical draws, and floating-point
        # noise can leave a truly PSD matrix looking marginally indefinite.
        # sample_params is the path that validates, via
        # _require_positive_semidefinite.
        theta = rng.multivariate_normal(
            self.population_mean,
            self.population_cov,
            _PEAK_LOG10_DRAWS,
            check_valid="ignore",
        )
        params = np.exp(theta)
        peaks = peak_day(self.model, params)
        return float(
            np.median(log10_concentration_rowwise(self.model, params, peaks[:, None]))
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

    def to_dict(self) -> dict:
        """
        Serialize this fit to a JSON/YAML-safe dict.

        Everything needed to simulate (``population_mean``, ``population_cov``,
        ``sigma``, ``censoring_limit``) is included. ``subject_params`` is
        deliberately omitted to keep the payload small; a deserialized fit can
        still be passed to ``simulate_shedding``, but has no per-subject table
        to inspect. This is the single serializer for a fit — the catalog's
        on-disk format is exactly one of these per fit, plus ``skipped``.

        Returns:
            A dict of plain Python/numpy-free types.
        """
        return {
            "dataset_id": self.dataset_id,
            "analyte": self.analyte,
            "biomarker": self.biomarker,
            "specimen": self.specimen,
            "reference_event": self.reference_event,
            "unit": self.unit,
            "gene_target": self.gene_target,
            "dose": self.dose,
            "vaccine_type": self.vaccine_type,
            "model": self.model,
            "method": self.method,
            "population_mean": [float(v) for v in self.population_mean],
            "population_cov": [[float(v) for v in row] for row in self.population_cov],
            "sigma": float(self.sigma),
            "censoring_limit": float(self.censoring_limit),
            "n_subjects": int(self.n_subjects),
            "n_measurements": int(self.n_measurements),
            "n_censored": int(self.n_censored),
            "n_excluded_subjects": int(self.n_excluded_subjects),
            "n_degenerate_subjects": int(self.n_degenerate_subjects),
            "pct_subjects_with_rise": float(self.pct_subjects_with_rise),
            "median_first_observed_day": float(self.median_first_observed_day),
            "n_dropped_measurements": int(self.n_dropped_measurements),
            "converged": bool(self.converged),
            "log_likelihood": float(self.log_likelihood),
            "aic": float(self.aic),
        }

    @classmethod
    def from_dict(cls, payload: dict) -> "SheddingFit":
        """
        Reconstruct a fit from ``to_dict``'s output.

        Optional keys default the same way a catalog predating that field
        would need to: ``n_degenerate_subjects`` to 0 (no such fit had any),
        ``pct_subjects_with_rise``/``median_first_observed_day`` to NaN
        ("unknown", not "zero"). ``subject_params`` is always ``None`` — it is
        never serialized, so there is nothing to restore it from.

        Args:
            payload: A dict produced by ``to_dict`` (or a YAML/JSON document
                loaded from one).

        Returns:
            A ``SheddingFit`` with ``subject_params is None``.
        """
        return cls(
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
            # Defaulted, not required: catalogs written before degeneracy
            # detection existed have no such key, and every fit in them
            # predates the concept.
            n_degenerate_subjects=int(payload.get("n_degenerate_subjects", 0)),
            # Likewise defaulted: NaN reads as "this catalog predates the rise
            # gate", which is honest, where 0.0 would assert that no subject
            # rose.
            pct_subjects_with_rise=float(
                payload.get("pct_subjects_with_rise", float("nan"))
            ),
            median_first_observed_day=float(
                payload.get("median_first_observed_day", float("nan"))
            ),
            n_dropped_measurements=int(payload["n_dropped_measurements"]),
            converged=bool(payload["converged"]),
            log_likelihood=float(payload["log_likelihood"]),
            aic=float(payload["aic"]),
        )


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

    No coefficient is ever clipped to the 1e-6 parameter floor: one that is
    non-positive, non-finite, or negligible falls back to a data-driven
    default instead.
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


def _degenerate_subjects(theta: np.ndarray, model: str) -> np.ndarray:
    """
    Flag subjects whose fit is an artifact rather than an estimate.

    Such a ``theta`` is not an estimate but a boundary or runaway solution, and
    averaging it into ``mean(theta_i)`` distorts the population summary out of
    all proportion. One subject pinned at ``log(1e-6) = -13.8`` is enough to
    turn a one-day half-life into a 278-day one; one subject with ``a0 = 142``
    is enough to turn a plausible peak into 10^18 gc/mL.

    Three ways to fail, deliberately symmetric — the check used to catch only
    collapse toward zero, judged against a physically-motivated floor, while
    the top end was left to the raw optimizer bound of ``exp(25) ~ 7.2e10``,
    which is so far past meaninglessness that nothing ever reached it:

    - **Collapsed**: any parameter at or below ``_DEGENERATE_PARAM``.
    - **Runaway decay**: an implied half-life below ``_MIN_HALF_LIFE_DAYS``,
      faster than roughly-daily sampling can resolve.
    - **Pinned**: any parameter at the optimizer's own upper bound. Retained for
      completeness, though the half-life test now fires long before this can.

    Args:
        theta: Fitted log-parameters, shape ``(n_subjects, k)``.
        model: ``"exponential"`` or ``"gamma"``, to locate the decay parameter.

    Returns:
        Boolean array of length ``n_subjects``, True where the subject's fit is
        degenerate by any of the three criteria.
    """
    collapsed = theta <= math.log(_DEGENERATE_PARAM)
    pinned = theta >= _THETA_BOUNDS[1] - _BOUND_TOLERANCE
    # ln(2)/a0 < _MIN_HALF_LIFE_DAYS, rearranged to compare on the log scale
    # theta is already expressed in. Located by name rather than by index so a
    # future change to PARAM_NAMES cannot silently test the wrong coordinate.
    decay = theta[:, PARAM_NAMES[model].index("a0")]
    runaway_decay = decay > math.log(math.log(2.0) / _MIN_HALF_LIFE_DAYS)
    return np.asarray((collapsed | pinned).any(axis=1) | runaway_decay)


def _fraction_observing_a_rise(observations: Observations) -> float:
    """
    Share of adequately-sampled subjects whose highest reading is not their first.

    Only subjects with at least ``_MIN_RISE_OBSERVATIONS`` uncensored positive
    readings at ``t > 0`` are judged: below that, "the maximum came later" is an
    artefact of having almost no readings rather than evidence about the shape
    of the trajectory. Censored points are excluded because a ``negative`` result
    carries no value to be the maximum.

    A subject whose maximum is tied between its first observation and a later
    one counts as *not* observing a rise — the peak is consistent with having
    already passed.

    Args:
        observations: Prepared observations for one analyte.

    Returns:
        The fraction in ``[0, 1]``, or NaN when no subject has enough readings
        to judge. NaN is deliberately not zero: it means "no evidence either
        way", and every comparison against it is False, so a caller gating on
        ``fraction >= threshold`` must spell out how it treats the undecidable
        case rather than silently getting one or the other.
    """
    usable = (~observations.censored) & (observations.times > 0)
    verdicts: list[bool] = []
    for subject in range(observations.n_subjects):
        mask = usable & (observations.subject_index == subject)
        if mask.sum() < _MIN_RISE_OBSERVATIONS:
            continue
        times = observations.times[mask]
        values = observations.values[mask]
        order = np.argsort(times, kind="stable")
        times, values = times[order], values[order]
        # argmax takes the earliest maximum, which is what makes a tie with the
        # first observation read as "no rise".
        verdicts.append(bool(times[int(np.argmax(values))] > times[0]))
    if not verdicts:
        return float("nan")
    return float(np.mean(verdicts))


def require_estimable_population(fit: SheddingFit) -> None:
    """
    Raise unless enough subjects survived for the population to mean anything.

    ``population_cov`` is a ``k x k`` covariance estimated from the retained
    subjects' ``theta``. With one subject it is all zeros, so every simulated
    individual is identical; with ``k`` or fewer it is rank-deficient, and the
    "population" is an artefact of having too little to average. The damage is
    not subtle — a two-subject gamma fit in the repository reported
    ``c_median = 148.6`` and a peak of 10^109 gc/mL. Requiring strictly more
    subjects than parameters gives the covariance a chance of being full rank.

    Applied by ``fit_shedding_models`` when building the catalog, deliberately
    **not** by ``fit_shedding_model`` itself. Fitting a single subject is a
    legitimate thing to ask for directly — validating the port against the
    published Rstan tutorial does exactly that — and such a fit is honest about
    itself, carrying a zero covariance. What must not happen is a fit like that
    being published in the catalog as though it described a population.

    The count is taken after degenerate-subject exclusion, because that is what
    actually feeds ``mu``/``Sigma``.

    Args:
        fit: A fit returned by ``fit_shedding_model``.

    Raises:
        SheddingDataError: With reason ``too_few_subjects_for_population``.
    """
    k = len(PARAM_NAMES[fit.model])
    retained = fit.n_subjects - fit.n_degenerate_subjects
    if retained > k:
        return
    excluded = (
        ""
        if not fit.n_degenerate_subjects
        else (
            f" ({fit.n_subjects} fitted, {fit.n_degenerate_subjects} excluded as "
            "degenerate)"
        )
    )
    raise SheddingDataError(
        f"Analyte {fit.analyte!r} has {retained} subject(s){excluded} feeding the "
        f"population summary of the {fit.model!r} model, which has {k} parameters. "
        f"At least {k + 1} are needed for the between-subject covariance to be "
        "estimable rather than rank-deficient; below that the 'population' is an "
        "artefact of averaging too few individuals, and simulating from it produces "
        "near-identical or wildly extrapolated individuals.",
        "too_few_subjects_for_population",
    )


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
        A ``SheddingFit``. Subjects whose fits are degenerate — collapsed onto
        the parameter bounds, or decaying faster than the sampling can resolve —
        are excluded from the population summary but retained in
        ``subject_params`` with ``degenerate`` set; ``n_degenerate_subjects``
        counts them.

    Raises:
        SheddingDataError: The analyte cannot be fitted (see ``reason``). In
            addition to every reason ``prepare_observations`` can raise, this
            adds two of its own: ``no_rise_observed``, when the gamma model is
            asked of data in which fewer than half the adequately-sampled
            subjects show a rise, and ``degenerate_fit``, when too many
            subjects' fits collapsed onto the parameter bounds for a population
            covariance to be estimable.

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

    rise_fraction = _fraction_observing_a_rise(observations)
    # `not (x >= t)` rather than `x < t` so that a NaN fraction — no subject had
    # enough readings to judge — refuses too. Absence of evidence that a rise
    # was ever observed is not evidence that fitting a rise is justified.
    if model == "gamma" and not (rise_fraction >= _MIN_RISE_FRACTION):
        observed = (
            "no subject had enough readings to judge"
            if math.isnan(rise_fraction)
            else f"only {100 * rise_fraction:.0f}% of subjects did"
        )
        raise SheddingDataError(
            f"The gamma model needs a rise to fit, and analyte {analyte!r} does "
            f"not show one: {observed}, against a required "
            f"{100 * _MIN_RISE_FRACTION:.0f}%. Sampling here appears to begin at "
            "or after peak shedding, which leaves b0 unidentifiable — the "
            "likelihood is monotone in it until the c0 > 0 constraint binds. Fit "
            "the exponential model instead; it is the appropriate model for "
            "post-peak sampling.",
            "no_rise_observed",
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

    # Subjects whose fits are artifacts — collapsed, pinned, or decaying faster
    # than the sampling can resolve — stay in subject_params so the raw fit
    # remains inspectable, but are kept out of the population summary, which
    # they would otherwise dominate.
    degenerate = _degenerate_subjects(theta, model)
    n_degenerate = int(degenerate.sum())
    retained = ~degenerate
    n_retained = int(retained.sum())
    # A single-subject fit legitimately yields a zero covariance and is allowed;
    # what is not allowed is *degeneracy* leaving too little behind. Hence
    # min(2, n) rather than a flat 2: with n == 1 the bar is that the one
    # subject survived, and only with n >= 2 does Sigma become estimable at all.
    if n_retained < min(2, n):
        raise SheddingDataError(
            f"{n_degenerate} of {n} subject(s) for analyte {analyte!r} produced "
            f"degenerate fits under the {model!r} model — collapsed onto the "
            "parameter bounds, or decaying faster than the sampling can resolve — "
            f"leaving {n_retained} usable subject(s), too few to estimate a "
            "population covariance. This usually means the model is not "
            "identifiable from this data: for the gamma model, typically because "
            "sampling began after peak shedding, so there is no rise phase from "
            "which to estimate b0.",
            "degenerate_fit",
        )
    if n_degenerate:
        warnings.warn(
            f"{n_degenerate} subject(s) excluded from the population summary of "
            f"analyte {analyte!r}: their fitted parameters are degenerate "
            "(collapsed onto the bounds, or an implied half-life below "
            f"{_MIN_HALF_LIFE_DAYS} days). They remain in subject_params, flagged "
            "by the 'degenerate' column.",
            UserWarning,
            stacklevel=2,
        )

    # When the typical retained subject was first sampled: each subject's own
    # earliest time, then the median across subjects. The median rather than the
    # overall minimum, because the minimum is the most generous statistic
    # available and hides how late most subjects started —
    # zuo2020alterations:stool_SARSCoV2_SymptomOnset has retained first days of
    # [5, 8, 16, 16, 18, 19, 22], where a minimum of 5 badly understates the
    # 16-day backward extrapolation behind its peak_log10.
    #
    # Censored observations count: a `negative` at day 2 still says the study
    # looked at day 2, and it still constrains the curve there. Taken over
    # retained subjects only, so it describes exactly the data behind
    # population_mean/population_cov.
    median_first_observed_day = float(
        np.median(
            [
                observations.times[observations.subject_index == subject].min()
                for subject in np.flatnonzero(retained)
            ]
        )
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
        pct_subjects_with_rise=100.0 * rise_fraction,
        median_first_observed_day=median_first_observed_day,
    )
