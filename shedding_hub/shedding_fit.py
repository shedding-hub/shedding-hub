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

from .shedding_models import PARAM_NAMES, validate_model

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
            ``no_positive_measurements``, ``unknown_analyte``. The catalog
            builder records this so a missing study reads as unsuitable, not
            as a bug.
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
            measurements, or leaves no subject with enough data.
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
