"""
Simulate shedding trajectories for synthetic infected individuals.

Intended for agent-based models: draw a cohort of individuals from a fitted
population distribution, then evaluate each one's shedding curve at whatever
times the simulation needs.
"""

from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd

from .shedding_models import log10_concentration_rowwise


def _resolve_incubation(
    incubation_period: Any, rng: np.random.Generator, n: int
) -> tuple[np.ndarray, bool]:
    """Return per-individual offsets in days, and whether a shift was applied."""
    if incubation_period is None:
        return np.zeros(n), False
    if callable(incubation_period):
        offsets = np.asarray(incubation_period(rng, n), dtype=float)
    elif np.isscalar(incubation_period):
        offsets = np.full(n, float(incubation_period))
    else:
        offsets = np.asarray(incubation_period, dtype=float)
    if offsets.shape != (n,):
        raise ValueError(
            "incubation_period must be None, a scalar, a callable, or an array of "
            f"length n_individuals ({n}); got shape {offsets.shape}."
        )
    return offsets, True


def simulate_shedding(
    source,
    *,
    n_individuals: int,
    times: Sequence[float] | np.ndarray,
    incubation_period: float | np.ndarray | Callable | None = None,
    include_measurement_error: bool = False,
    seed: int | None = None,
) -> pd.DataFrame:
    """
    Simulate shedding trajectories for synthetic individuals.

    Args:
        source: A ``SheddingFit`` or ``SheddingEnsemble``. Each individual's
            parameters are drawn from its population distribution.
        n_individuals: Number of individuals to simulate.
        times: Times at which to evaluate each trajectory. Measured from the
            fit's reference event, or from infection when ``incubation_period``
            is given.
        incubation_period: Days from infection to the reference event. ``None``
            leaves output in reference-event time. A scalar applies to everyone;
            an array of length ``n_individuals`` or a callable
            ``(rng, n) -> array`` gives each individual its own, which adds
            realistic timing variability across the cohort.
        include_measurement_error: Add ``N(0, sigma)`` assay noise on the log10
            scale. Off by default: an agent-based model wants the true shed
            concentration, and assay noise is a property of sampling rather than
            of the host.
        seed: Seed for a ``numpy`` generator, making runs reproducible.

    Returns:
        A tidy DataFrame with columns ``individual_id``, ``time``,
        ``log10_value``, ``value``, ``detected``, ``source_dataset_id``.
        ``detected`` is whether the value reaches the censoring limit. Values
        below the limit are reported as-is rather than clipped, so downstream
        mass-balance calculations stay correct. Under the gamma model,
        non-positive times yield NaN.

        ``result.attrs`` records ``time_origin``, ``incubation_applied``,
        ``model``, ``unit``, ``biomarker``, and ``specimen``.

    Note:
        The exponential model is defined for negative times but grows without
        bound going backwards, so simulating before the reference event
        extrapolates into implausible concentrations. Prefer ``times >= 0``.
    """
    if n_individuals < 1:
        raise ValueError("n_individuals must be at least 1")

    rng = np.random.default_rng(seed)
    params, sources = source.sample_params(rng, n_individuals)
    times = np.asarray(times, dtype=float)
    offsets, incubation_applied = _resolve_incubation(
        incubation_period, rng, n_individuals
    )

    shifted = times[None, :] - offsets[:, None]
    log10_values = log10_concentration_rowwise(source.model, params, shifted)

    if include_measurement_error:
        log10_values = log10_values + rng.normal(
            0.0, source.sigma, size=log10_values.shape
        )

    n_times = times.size
    frame = pd.DataFrame(
        {
            "individual_id": np.repeat(np.arange(n_individuals), n_times),
            "time": np.tile(times, n_individuals),
            "log10_value": log10_values.ravel(),
            "source_dataset_id": np.repeat(sources, n_times),
        }
    )
    frame["value"] = np.power(10.0, frame["log10_value"])
    detected = frame["log10_value"] >= source.censoring_limit
    frame["detected"] = detected.fillna(False).astype(bool)
    frame = frame[
        [
            "individual_id",
            "time",
            "log10_value",
            "value",
            "detected",
            "source_dataset_id",
        ]
    ]
    frame.attrs = {
        "time_origin": "infection" if incubation_applied else source.reference_event,
        "incubation_applied": incubation_applied,
        "model": source.model,
        "unit": source.unit,
        "biomarker": getattr(source, "biomarker", None),
        "specimen": getattr(source, "specimen", None),
    }
    return frame
