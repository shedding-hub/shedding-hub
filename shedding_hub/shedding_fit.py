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

from .shedding_models import LN10, PARAM_NAMES, theta_to_params, validate_model

CENSORING_MARGIN = 0.01
NEGATIVE_VALUE = "negative"

# Ct is affine in log10 concentration (Ct = alpha - beta * log10 C), so the
# shedding models describe it unchanged once it is negated -- Ct falls as
# shedding rises -- and offset to keep fitted levels positive.
#
# The offset is one constant for every analyte, not each study's own cutoff.
# Recorded cutoffs run 37 to 41, so anchoring per study would put two studies
# measuring identical samples up to 4 cycles apart on height alone. 40 sits
# above the observed Ct median of 31 and above 95% of all readings, so fitted
# peak heights -- which occur at LOW Ct -- stay comfortably positive.
CT_REFERENCE = 40.0


def _to_response(value: float, value_type: str) -> float:
    """
    Map a reported measurement onto the scale the models are fitted on.

    Concentrations are fitted on log10. Cycle thresholds are fitted as cycles
    below ``CT_REFERENCE``, which is decreasing in Ct and therefore increasing
    in viral load, exactly like a log10 concentration.
    """
    if value_type == "ct":
        return CT_REFERENCE - float(value)
    return math.log10(float(value))


# Fecal-strength / normalization indicators, not pathogens shed by infected
# people. They have no time-since-infection trajectory, so fitting a shedding
# curve to them is meaningless regardless of how much data is available.
NON_PATHOGEN_BIOMARKERS = frozenset({"crAssphage", "PMMoV", "mtDNA"})

# Readings earlier than this many days before the reference event are discarded.
#
# Every measurement in the repository earlier than about day -3 is reported
# `negative`: 5,101 of them, across the bins [-60,-30), [-30,-14), [-14,-7),
# [-7,-5) and [-5,-3), all 100% censored. There are only 151 detected
# measurements at negative times at all, and the two earliest sit at exactly
# day -5.
#
# Those censored points are not neutral to a decay-only model. The exponential
# curve peaks at t = 0 by construction, so at day -53 it predicts
# c0 * exp(53 * a0) -- astronomically high -- and a reading "below the limit"
# there is near-impossible under the model. The censored likelihood then pushes
# hard on a0 and c0 to accommodate it: tsang2016individual NPSOPS had its
# pre-reference-event readings sitting a median 4.01 log10 below its own fitted
# curve, against -0.26 for its 1,100-odd later ones.
#
# -5 rather than -3 so that no detected measurement in the repository is ever
# discarded. It is a parameter, so a future dataset with earlier positives does
# not silently lose them.
#
# The gamma model is unaffected: it already drops everything at t <= 0, which is
# stricter. Were a shifted gamma ever added, this cutoff should NOT apply to it
# -- a censored reading at day -53 means "shedding had not started yet", which is
# precisely what would identify the shift.
_MIN_TIME_DAYS = -5.0

# How far below a subject's earliest retained reading its shedding onset must
# sit, under gamma_shifted. The curve is undefined at t <= t0 and dives steeply
# just above it, so the onset needs clearance from the first reading it has to
# explain. Half a day is small enough not to distort a real onset and large
# enough to keep ln(t - t0) finite in the optimizer's arithmetic.
_ONSET_MARGIN_DAYS = 0.5


class SheddingDataError(ValueError):
    """
    Raised when an analyte cannot be fitted.

    Attributes:
        reason: Machine-readable cause, one of
            ``non_pathogen_biomarker``, ``too_few_subjects``,
            ``no_positive_measurements``, ``no_data_after_reference_event``,
            ``no_pre_event_readings``, ``unknown_analyte``,
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

    Examples:
        >>> import shedding_hub as sh
        >>> error = sh.SheddingDataError(
        ...     "no positive measurements", "no_positive_measurements"
        ... )
        >>> error.reason
        'no_positive_measurements'
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
    # Which scale ``values`` and ``censoring_limit`` are on: log10 concentration,
    # or cycles below CT_REFERENCE. Carried on the observations rather than
    # re-derived downstream, so a plot or a fit cannot disagree with the fitter
    # about what its own numbers mean.
    value_type: str = "concentration"
    subject_ids: list = field(default_factory=list)
    n_subjects: int = 0
    n_excluded_subjects: int = 0
    n_dropped_measurements: int = 0
    # The plottable subset of what was dropped: readings with a usable time and
    # value that a model-specific rule discarded. Recorded so a diagnostic plot
    # can mark them without re-deriving the rules, which would drift.
    # ``dropped_values`` is log10, NaN where the reading was censored, matching
    # ``values``. Readings with no usable time or value are counted in
    # ``n_dropped_measurements`` but cannot be placed on a plot and are absent.
    dropped_times: np.ndarray = field(default_factory=lambda: np.empty(0))
    dropped_values: np.ndarray = field(default_factory=lambda: np.empty(0))


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


def _resolve_censoring_limit(
    analyte_spec: dict,
    observed: np.ndarray,
    value_type: str = "concentration",
) -> float:
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

    Assumes ``observed`` is non-empty; ``prepare_observations`` guarantees
    this by raising ``no_positive_measurements`` itself before ever calling here.

    For a cycle-threshold analyte the declared limit is itself a Ct, so it is
    transformed the same way the observations are and the resolved limit is
    ``CT_REFERENCE - cutoff`` — zero at a cutoff of 40, negative where an assay
    runs to 41, positive where it stops at 37. The fallback branch needs no
    special case: ``observed`` has already been transformed, so "just below the
    smallest response" means the same thing on either scale.
    """
    for key in ("limit_of_quantification", "limit_of_detection"):
        limit = _numeric_limit(analyte_spec.get(key))
        if limit is not None:
            return _to_response(limit, value_type)

    smallest = float(observed.min())
    fallback = smallest - CENSORING_MARGIN
    scale = "cycles below reference" if value_type == "ct" else "log10"
    warnings.warn(
        "Falling back to a censoring limit of "
        f"{fallback:.4g} ({scale}) because no limit of quantification or "
        "detection is declared for this analyte.",
        UserWarning,
        stacklevel=2,
    )
    return fallback


def _record_dropped(
    measurement: dict,
    time: float,
    times: list,
    values: list,
    value_type: str = "concentration",
) -> None:
    """Note a discarded reading, if it can be placed on a plot at all."""
    value = measurement.get("value")
    if isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0:
        times.append(time)
        values.append(_to_response(float(value), value_type))
    elif value == NEGATIVE_VALUE:
        times.append(time)
        values.append(float("nan"))


def prepare_observations(
    dataset: dict,
    analyte: str,
    model: str,
    *,
    min_observations: int | None = None,
    min_time: float = _MIN_TIME_DAYS,
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
        model: ``"exponential"``, ``"gamma"`` or ``"gamma_shifted"``.
        min_observations: Minimum usable measurements a subject must have to be
            retained. Defaults to the number of per-subject parameters (3 for
            gamma, 2 for exponential). ``sigma`` is shared across subjects, so a
            subject does not need residual degrees of freedom of its own.
            Independently of this count, a subject with no positive measurement
            at all is excluded: its readings locate no curve, only an upper
            bound, so its fitted parameters are arbitrary.
        min_time: Earliest time, in days from the reference event, that a
            reading may carry and still be used. See ``_MIN_TIME_DAYS`` for why
            the default is -5 and why it does not affect the gamma model.

    Returns:
        An ``Observations`` instance with subject indices renumbered
        contiguously. For a cycle-threshold analyte, ``values`` are cycles
        below ``CT_REFERENCE`` rather than log10 concentrations, and
        ``value_type`` says which.

    Raises:
        SheddingDataError: The analyte is unknown, is a non-pathogen
            indicator biomarker, has no positive measurements, or leaves no
            subject with enough data — reasons ``unknown_analyte``,
            ``non_pathogen_biomarker``, ``no_positive_measurements`` and
            ``too_few_subjects`` respectively, plus
            ``no_data_after_reference_event`` when every usable reading
            falls at or before day 0. Note that ``too_few_subjects`` here means no subject
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

    # Cycle thresholds are affine in log10 concentration, so both models
    # describe them once the response is transformed. See ``_to_response``.
    value_type = "ct" if _is_ct_unit(analyte_spec.get("unit")) else "concentration"

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
    dropped_times: list[float] = []
    dropped_values: list[float] = []

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
            gamma_drops_it = model == "gamma" and time <= 0
            if gamma_drops_it or time < min_time:
                n_dropped += 1
                _record_dropped(
                    measurement, time, dropped_times, dropped_values, value_type
                )
                continue
            value = measurement.get("value")
            if isinstance(value, str):
                if value == NEGATIVE_VALUE and not (
                    model == "gamma_shifted" and time <= 0
                ):
                    times.append(time)
                    values.append(np.nan)
                    censored.append(True)
                else:
                    # A censored reading at or before the reference event is
                    # dropped under gamma_shifted. Its curve dives toward minus
                    # infinity as t approaches t0, so "below the limit" there is
                    # explained for free and t0 becomes a support parameter
                    # pulled onto its own bound. A *detected* reading at the same
                    # time is kept, and repels t0 instead: a diving curve
                    # mispredicts a measured value badly.
                    n_dropped += 1
                    if value == NEGATIVE_VALUE:
                        dropped_times.append(time)
                        dropped_values.append(np.nan)
                continue
            if not isinstance(value, (int, float)) or value <= 0:
                n_dropped += 1
                continue
            times.append(time)
            values.append(_to_response(float(value), value_type))
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

    # Checked before the retention filter below, not after it. Every retained
    # subject now has a positive reading by construction, so asking afterwards
    # could only ever report the less specific 'too_few_subjects' for an analyte
    # whose real problem is that nothing was ever detected in it.
    if all(censored for subject in per_subject for censored in subject["censored"]):
        raise SheddingDataError(
            f"Analyte {analyte!r} has no positive measurements to fit.",
            "no_positive_measurements",
        )

    # Filter subject_ids and per_subject together, from a single zipped list,
    # so the retention predicate appears exactly once: applying it separately
    # to each list could silently desync subject_ids from the arrays below.
    # The two exclusion reasons are counted apart so each can say what it
    # actually means, but both land in n_excluded_subjects.
    retained_pairs = []
    n_too_few = 0
    n_no_positive = 0
    for subject_id, subject in zip(subject_ids, per_subject):
        if len(subject["times"]) < min_observations:
            n_too_few += 1
        elif all(subject["censored"]):
            n_no_positive += 1
        else:
            retained_pairs.append((subject_id, subject))
    n_excluded = n_too_few + n_no_positive
    retained_ids = [subject_id for subject_id, _ in retained_pairs]
    retained = [subject for _, subject in retained_pairs]

    if n_too_few:
        warnings.warn(
            f"{n_too_few} subject(s) excluded from the {analyte!r} fit for having "
            f"fewer than {min_observations} usable measurements.",
            UserWarning,
            stacklevel=2,
        )
    if n_no_positive:
        warnings.warn(
            f"{n_no_positive} subject(s) excluded from the {analyte!r} fit for "
            "having no positive measurement: every reading was 'negative'. Such a "
            "subject constrains the curve only to stay below the limit, which "
            "every curve that does so satisfies equally, so the optimizer returns "
            "an arbitrary point estimate that this two-stage estimator would then "
            "average into the population summary at full weight.",
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

    # gamma_shifted exists to use readings at or before the reference event, and
    # is only defensible where there are some. Without one, t0 has nothing to
    # locate and merely absorbs curve shape: on woelfel2020virological stool the
    # onset came back at +0.97 days -- after the reference event -- and AIC
    # preferred the plain gamma model, 223.8 against 239.8, on the identical 79
    # observations.
    #
    # This is also what keeps the two rise-and-fall models comparable at all.
    # Where gamma_shifted is admitted it is fitted to more observations than
    # gamma, so their AICs are not comparable and the choice between them is
    # made by data availability -- this gate -- rather than by fit statistic.
    if model == "gamma_shifted" and not any(
        time <= 0 for subject in retained for time in subject["times"]
    ):
        raise SheddingDataError(
            f"Analyte {analyte!r} has no detected reading at or before its "
            "reference event, so a shifted onset has nothing to locate and "
            "would only absorb curve shape. Fit the 'gamma' model instead.",
            "no_pre_event_readings",
        )

    # Both models describe shedding from the reference event onwards -- the
    # exponential decays from it, the gamma rises and falls after it -- so an
    # analyte sampled only up to it cannot constrain either.
    #
    # Checked after the retention filter, unlike no_positive_measurements above.
    # A cross-sectional study sampled once per subject at day 0 trips both this
    # and too_few_subjects, and "no subject has enough measurements" is the more
    # useful diagnosis there; an analyte that reaches here has subjects with
    # ample readings and genuinely stops at the reference event.
    #
    # jones2021estimating swab_SARSCoV2_confirmation is that case: sampled days
    # -7 to 0, it optimized to convergence and was published with 1990 of its
    # 2075 subjects degenerate, a sigma of 5.41 against a catalog median of
    # 0.84, and a median individual 1.26 log10 below its own censoring limit.
    if not any(time > 0 for subject in retained for time in subject["times"]):
        raise SheddingDataError(
            f"Analyte {analyte!r} has no measurement after its reference event "
            "(every usable reading is at or before day 0), so there is no "
            "post-event trajectory for either model to describe.",
            "no_data_after_reference_event",
        )

    subject_index = np.concatenate(
        [np.full(len(s["times"]), i, dtype=int) for i, s in enumerate(retained)]
    )
    times_array = np.concatenate([np.asarray(s["times"], float) for s in retained])
    values_array = np.concatenate([np.asarray(s["values"], float) for s in retained])
    censored_array = np.concatenate([np.asarray(s["censored"], bool) for s in retained])

    # Non-empty: retention requires at least one positive per subject, and the
    # analyte-wide check above already rejected the case where there are none.
    observed = values_array[~censored_array]
    censoring_limit = _resolve_censoring_limit(analyte_spec, observed, value_type)

    return Observations(
        subject_index=subject_index,
        times=times_array,
        values=values_array,
        censored=censored_array,
        censoring_limit=censoring_limit,
        value_type=value_type,
        subject_ids=retained_ids,
        n_subjects=len(retained),
        n_excluded_subjects=n_excluded,
        n_dropped_measurements=n_dropped,
        dropped_times=np.asarray(dropped_times, dtype=float),
        dropped_values=np.asarray(dropped_values, dtype=float),
    )


import pandas as pd
from scipy import optimize
from scipy.stats import norm

from .shedding_models import (
    POPULATION_COORDS,
    from_population_coords,
    half_life_days,
    log10_concentration_pointwise,
    log10_concentration_rowwise,
    peak_day,
    to_population_coords,
)

# The historical positivity floor for a natural-scale parameter is 1e-6.
# Nothing is ever *initialized* there — see _DEFAULT_A0 for why that was a
# bug — but it remains the reference point that _DEGENERATE_PARAM below is
# judged against, and is not otherwise given a named constant since nothing
# else in this module reads it directly.
_THETA_BOUNDS = (-25.0, 25.0)
_LOG_SIGMA_BOUNDS = (-10.0, 5.0)

# How many times the optimizer may be handed its own result and told to keep
# going. Six is well past what any fit in the repository needs -- the worst
# observed is a second round -- and exists so that a pathological surface
# terminates instead of running forever.
_MAX_OPTIMIZER_ROUNDS = 6

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


# How far above the highest observed concentration a subject's implied peak may
# sit before it is treated as extrapolation rather than estimate, in log10 units.
# Three is a thousandfold: enough headroom that a study sampling a few days after
# the reference event still summarizes every real subject, tight enough to
# exclude the hundred-half-life backward extrapolations that late-starting
# studies produce. See ``_over_extrapolated_subjects``.
_MAX_PEAK_ABOVE_OBSERVED = 3.0

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


def validate_dispersion(dispersion: float) -> float:
    """Check a dispersion factor and return it as a float.

    Shared by ``SheddingFit`` and ``SheddingEnsemble`` so both reject the same
    values with the same message.
    """
    dispersion = float(dispersion)
    if not np.isfinite(dispersion) or dispersion < 0:
        raise ValueError(
            f"dispersion must be a finite, non-negative number; got {dispersion!r}. "
            "It scales the between-subject covariance by dispersion**2, so 1.0 "
            "leaves the fitted population alone and 0.0 makes every individual "
            "the median one."
        )
    return dispersion


def _scaled(cov: np.ndarray, dispersion: float) -> np.ndarray:
    """The covariance to draw from, after applying ``dispersion``."""
    return cov * validate_dispersion(dispersion) ** 2


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

    Examples:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> fit = catalog.select(
        ...     dataset_id='woelfel2020virological', analyte='stool', model='gamma'
        ... )
        >>> fit.param_names
        ('a0', 'b0', 'c0')
        >>> round(fit.peak_day, 2)
        1.18
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
    def population_coords(self) -> tuple[str, ...]:
        """Names of the coordinates ``population_mean``/``population_cov`` use."""
        return POPULATION_COORDS[self.model]

    @property
    def median_params(self) -> np.ndarray:
        """
        Parameters of the median individual.

        ``population_mean`` is the mean of a normal, hence also its median, so
        mapping it back through ``from_population_coords`` gives the parameters
        of the median individual in each summarized coordinate. Note that the
        median individual's trajectory is not the population's mean trajectory;
        to aggregate load across a cohort, simulate rather than scaling this up.
        """
        return from_population_coords(self.model, self.population_mean[None, :])[0]

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
        params = from_population_coords(self.model, theta)
        peaks = peak_day(self.model, params)
        return float(
            np.median(log10_concentration_rowwise(self.model, params, peaks[:, None]))
        )

    @property
    def half_life_days(self) -> float:
        return float(half_life_days(self.model, self.median_params[None, :])[0])

    def sample_params(
        self, rng: np.random.Generator, n: int, dispersion: float = 1.0
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Draw ``n`` individuals' natural-scale parameters.

        Args:
            rng: Generator to draw from.
            n: Number of individuals.
            dispersion: Scales the between-subject covariance by
                ``dispersion ** 2``, so the population's spread scales by
                ``dispersion`` while its centre and correlation structure are
                untouched. See ``simulate_shedding`` for why you might want it
                below 1.

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
        theta = rng.multivariate_normal(
            self.population_mean, _scaled(cov, dispersion), n
        )
        return (
            from_population_coords(self.model, theta),
            np.full(n, self.dataset_id, dtype=object),
        )

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
            # Recorded so a catalog cannot be silently misread: the coordinates
            # are the same *length* under any convention, so without this an
            # older file loads cleanly and yields wrong curves.
            "population_coords": list(self.population_coords),
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
        model = payload["model"]
        expected = list(POPULATION_COORDS[model])
        if payload.get("population_coords") != expected:
            raise ValueError(
                f"This {model!r} fit records its population summary in "
                f"{payload.get('population_coords')!r}, but this version of "
                f"shedding_hub reads {expected!r}. Both are the same length, so "
                "loading it anyway would silently produce wrong curves rather "
                "than fail. Rebuild the catalog with `make catalog`."
            )
        return cls(
            model=model,
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
        if model == "gamma_shifted":
            # Start the onset one margin below the earliest reading it must
            # explain, which is both feasible and the least presumptuous guess:
            # shedding began shortly before this subject was first sampled.
            onset = float(np.min(times)) - _ONSET_MARGIN_DAYS
            return np.array([np.log(a0), np.log(_DEFAULT_B0), np.log(c0), onset])
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


def _over_extrapolated_subjects(
    theta: np.ndarray,
    model: str,
    observations: Observations,
    margin: float = _MAX_PEAK_ABOVE_OBSERVED,
) -> np.ndarray:
    """
    Flag subjects whose implied peak sits far above anything the study observed.

    A peak is only measured when sampling reaches it. Many studies begin days or
    weeks after the reference event, and the exponential model's peak is at
    ``t = 0`` by definition, so a subject with a steep fitted decay is
    extrapolated backwards through many half-lives. In
    ``fajnzylber2020sars`` nasopharyngeal, one subject with a 0.30-day half-life
    first sampled on day 30 implies ``10**33`` gc/mL — a hundred half-lives of
    extrapolation — while its own highest reading was ``10**2.7`` and the whole
    analyte never exceeded ``10**5.5``.

    Such a value is not an estimate of a concentration; it is the functional form
    continued past every observation that constrains it. Averaging it into the
    population summary is what left 87% of that analyte's observations *below*
    its median-individual curve.

    The threshold is referenced to the data rather than fixed, since what counts
    as absurd depends on what the assay can see: more than
    ``_MAX_PEAK_ABOVE_OBSERVED`` log10 — a thousandfold — above the highest
    concentration the analyte ever recorded. Measured over the repository this
    flags 66 subjects across 34 of 92 fits, while leaving well-sampled fits
    untouched (``tsang2016individual`` NPSOPS loses 4 of 440 and does not move).

    Why this could not be left to the log scale: while the exponential model was
    summarized in ``log c0``, the logarithm compressed these subjects enough to
    hide them. Summarizing in ``peak_log10`` — which is what keeps simulated
    concentrations from being double-exponential — makes the coordinate linear in
    log10, so one such subject dominates the mean outright. The two changes
    belong together.

    Args:
        theta: Fitted log-parameters, shape ``(n_subjects, k)``.
        model: ``"exponential"``, ``"gamma"`` or ``"gamma_shifted"``.
        observations: The observations the subjects were fitted to.
        margin: Log10 units of headroom above the highest observed
            concentration. Defaults to ``_MAX_PEAK_ABOVE_OBSERVED``.

    Returns:
        Boolean array of length ``n_subjects``, True where the subject's implied
        peak is unsupportable. All False when the analyte has no positive
        reading to reference, which leaves the judgement to the other checks.
    """
    positive = ~observations.censored
    if not positive.any():
        return np.zeros(theta.shape[0], dtype=bool)
    ceiling = float(np.nanmax(observations.values[positive])) + margin
    # Indexed by name, not by position. `peak_log10` is last for `exponential`
    # and `gamma`, so `[:, -1]` read correctly for those two -- but
    # `gamma_shifted` ends ('...', 'peak_log10', 't0'), so the same expression
    # returned an onset in days. Comparing a t0 of about -3 against a ceiling of
    # about 9 log10 is never true, and the gate did nothing at all for every
    # gamma_shifted fit.
    peak = POPULATION_COORDS[model].index("peak_log10")
    heights = to_population_coords(model, theta_to_params(model, theta))[:, peak]
    return np.asarray(heights > ceiling)


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
        model: ``"exponential"``, ``"gamma"`` or ``"gamma_shifted"``, to
            locate the decay parameter.

    Returns:
        Boolean array of length ``n_subjects``, True where the subject's fit is
        degenerate by any of the three criteria.
    """
    # The floor and ceiling tests apply to the log-parameters only. Under
    # gamma_shifted the last coordinate is t0, an absolute time: an onset of
    # -6.0 days is entirely ordinary yet sits below log(1e-2) = -4.6, so
    # including it here would condemn every subject whose shedding began more
    # than 4.6 days before the reference event. Its own boundary -- landing on
    # the per-subject bound below its first reading -- is judged in
    # ``_onset_on_its_bound``, which needs the observations to know where that
    # bound was.
    scales = theta[:, :3] if model == "gamma_shifted" else theta
    collapsed = scales <= math.log(_DEGENERATE_PARAM)
    pinned = scales >= _THETA_BOUNDS[1] - _BOUND_TOLERANCE
    # ln(2)/a0 < _MIN_HALF_LIFE_DAYS, rearranged to compare on the log scale
    # theta is already expressed in. Located by name rather than by index so a
    # future change to PARAM_NAMES cannot silently test the wrong coordinate.
    decay = theta[:, PARAM_NAMES[model].index("a0")]
    runaway_decay = decay > math.log(math.log(2.0) / _MIN_HALF_LIFE_DAYS)
    return np.asarray((collapsed | pinned).any(axis=1) | runaway_decay)


def _fraction_observing_a_rise(
    observations: Observations, model: str = "gamma"
) -> float:
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
    # Judged over the readings the model will use. gamma_shifted retains
    # detected readings at t <= 0, and those are exactly where a rise crossing
    # the reference event shows itself, so excluding them would hide the very
    # shape the gate is looking for.
    usable = ~observations.censored
    if model != "gamma_shifted":
        usable = usable & (observations.times > 0)
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

    params = theta_to_params(model, theta)[observations.subject_index]
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
    min_time: float = _MIN_TIME_DAYS,
    max_peak_above_observed: float = _MAX_PEAK_ABOVE_OBSERVED,
) -> SheddingFit:
    """
    Fit a shedding model to one analyte by censored maximum likelihood.

    Args:
        dataset: Dataset dictionary from ``load_dataset``.
        analyte: Key into ``dataset["analytes"]``.
        model: ``"exponential"``, ``"gamma"`` or ``"gamma_shifted"``.
        min_time: Earliest time, in days from the reference event, a reading may
            carry and still be used. Passed to ``prepare_observations``.
        max_peak_above_observed: How far above the analyte's highest observed
            concentration a subject's implied peak may sit, in log10 units,
            before it is excluded from the population summary as extrapolation.
            See ``_over_extrapolated_subjects``. Lower it to judge the catalog
            under a stricter reading of what counts as supportable.
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

    Examples:
        Fitting is a joint optimization over every subject's parameters. One
        analyte of one study, as below, takes about a second; a
        repository-wide build over every analyte of every study is what takes
        minutes.

        >>> import shedding_hub as sh
        >>> data = sh.load_dataset(
        ...     'woelfel2020virological', local='./data'
        ... )
        >>> fit = sh.fit_shedding_model(data, analyte='stool', model='gamma')
        >>> fit.model
        'gamma'
    """
    validate_model(model)
    observations = prepare_observations(
        dataset, analyte, model, min_observations=min_observations, min_time=min_time
    )

    rise_fraction = _fraction_observing_a_rise(observations, model)
    # `not (x >= t)` rather than `x < t` so that a NaN fraction — no subject had
    # enough readings to judge — refuses too. Absence of evidence that a rise
    # was ever observed is not evidence that fitting a rise is justified.
    if model in ("gamma", "gamma_shifted") and not (
        rise_fraction >= _MIN_RISE_FRACTION
    ):
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
    if model == "gamma_shifted":
        # t0 must sit below every reading the subject has, or its own curve is
        # undefined at its own observations. Imposed as a per-subject upper
        # bound rather than by reparameterizing, so t0 stays interpretable as an
        # absolute time and the population summary can average it directly.
        onset_index = PARAM_NAMES[model].index("t0")
        for subject in range(n):
            mine = observations.times[observations.subject_index == subject]
            bounds[subject * k + onset_index] = (
                _THETA_BOUNDS[0],
                float(mine.min()) - _ONSET_MARGIN_DAYS,
            )

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
    # gamma_shifted needs more than the others: its fourth coordinate is
    # bounded per subject, and the constrained surface takes longer to settle.
    # Measured on a 25-subject synthetic fit (101 parameters), it converged at
    # 148,716 evaluations -- past the 101,000 that a multiplier of 1000 allows.
    #
    # That measurement was taken on one machine, and a count measured on one
    # machine is not a property of the problem. The same synthetic fit exhausts
    # a multiplier of 2000 on Linux, where it is still descending: it reaches
    # log-likelihood 9.828 against Windows' 9.749 -- a *better* optimum -- and
    # then reports failure for having run out of allowance rather than for
    # having stopped improving. Heavy censoring flattens this surface near its
    # optimum, so exactly where L-BFGS-B satisfies ftol depends on the BLAS
    # underneath, and any fixed cap will be generous on one platform and mean on
    # another.
    #
    # So the cap is no longer asked to be right. It is a chunk size, and when a
    # round ends *only* because the chunk ran out (status 1), the optimizer is
    # handed its own result and told to keep going. Convergence then depends on
    # the problem rather than on the machine, and the reported flag means what
    # it says. Restarting also resets the limited-memory Hessian approximation,
    # which is often what a stalled L-BFGS-B run needs.
    #
    # status 2 -- a genuine breakdown, such as a line-search failure -- is not
    # retried: repeating it would only burn evaluations to fail identically.
    #
    # The per-model multiplier stays as it was, deliberately. Every fit that
    # already converged still runs one round with the same budget it always
    # had, so its result is bit-identical and the shipped catalog can only move
    # where it was previously reporting non-convergence.
    multiplier = 2000 if model == "gamma_shifted" else 1000
    max_evaluations = max(15000, multiplier * n_parameters)
    options = {
        "maxfun": max_evaluations,
        "maxiter": max_evaluations,
        "ftol": 1e-6,
    }
    result = optimize.minimize(
        _negative_log_likelihood,
        x0,
        args=(model, observations),
        method="L-BFGS-B",
        bounds=bounds,
        options=options,
    )
    rounds = 1
    while not result.success and result.status == 1 and rounds < _MAX_OPTIMIZER_ROUNDS:
        result = optimize.minimize(
            _negative_log_likelihood,
            result.x,
            args=(model, observations),
            method="L-BFGS-B",
            bounds=bounds,
            options=options,
        )
        rounds += 1

    if not result.success:
        warnings.warn(
            f"Optimizer did not converge for analyte {analyte!r} "
            f"({result.message}). The fit is returned with converged=False.",
            UserWarning,
            stacklevel=2,
        )

    theta = result.x[: n * k].reshape(n, k)
    sigma = float(np.exp(result.x[-1]))

    # Subjects whose fits are artifacts — collapsed, pinned, decaying faster than
    # the sampling can resolve, or implying a peak far above anything the study
    # observed — stay in subject_params so the raw fit remains inspectable, but
    # are kept out of the population summary, which they would otherwise
    # dominate. The last of those is judged against the data rather than against
    # the parameter bounds, which is why it needs the observations.
    degenerate = _degenerate_subjects(theta, model) | _over_extrapolated_subjects(
        theta, model, observations, max_peak_above_observed
    )
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

    # Summarized in population coordinates, not in the log-parameters the
    # optimizer works in. For the gamma model those are not the same space, and
    # averaging the log-parameters coordinate-wise lands off the ridge the
    # subjects lie on — see ``to_population_coords``. The exponential model's
    # coordinates are its log-parameters, so this is an identity for it.
    kept = to_population_coords(model, theta_to_params(model, theta[retained]))
    population_mean = kept.mean(axis=0)
    population_cov = (
        np.cov(kept, rowvar=False, ddof=1) if n_retained > 1 else np.zeros((k, k))
    )
    population_cov = np.atleast_2d(population_cov)

    subject_params = pd.DataFrame(
        theta_to_params(model, theta), columns=list(PARAM_NAMES[model])
    )
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
