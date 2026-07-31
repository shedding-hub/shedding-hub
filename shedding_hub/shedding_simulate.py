"""
Simulate shedding trajectories for synthetic infected individuals.

Intended for agent-based models: draw a cohort of individuals from a fitted
population distribution, then evaluate each one's shedding curve at whatever
times the simulation needs.
"""

import warnings
from typing import Any, Callable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from .shedding_models import log10_concentration_rowwise
from .shedding_select import classify_reference_event


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
    dispersion: float = 1.0,
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

            Only meaningful for a reference event that is a natural-history
            landmark (symptom onset). Applying it to an administrative event
            (enrollment, confirmation date, hospital admission, treatment) or to
            the exposure itself (inoculation, vaccination) warns, and
            ``attrs["time_origin"]`` records ``"<event>_shifted"`` rather than
            ``"infection"``.
        include_measurement_error: Add ``N(0, sigma)`` assay noise on the log10
            scale. Off by default: an agent-based model wants the true shed
            concentration, and assay noise is a property of sampling rather than
            of the host.
        dispersion: Scales the between-subject covariance by ``dispersion ** 2``,
            so the cohort's spread scales by ``dispersion`` while its centre and
            correlation structure stay put. ``1.0`` (default) simulates the
            fitted population as estimated.

            Reach for a value below 1 when a handful of agents dominate the
            cohort's total shed load. Two-stage estimation does not shrink
            individual estimates toward the population mean, so
            ``population_cov`` carries within-subject estimation error on top of
            true between-subject variance and every simulated cohort is somewhat
            over-dispersed — the more so the fewer observations each subject
            has. The bias only ever runs one way, which is what makes a
            shrinkage factor a defensible correction rather than a fudge, but
            there is no automatic way to choose it: it is a judgement about how
            much of the fitted spread is real.
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
    params, sources = source.sample_params(rng, n_individuals, dispersion)
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
    event = source.reference_event
    event_class = classify_reference_event(event)
    time_origin = event
    if incubation_applied:
        if event_class == "landmark":
            time_origin = "infection"
        else:
            time_origin = f"{event}_shifted" if event else "shifted reference event"
            if event_class == "exposure":
                warnings.warn(
                    f"{event!r} is already the exposure, so there is no "
                    "incubation period to bridge: shifting moves the origin to "
                    f"before the exposure itself. time_origin is recorded as "
                    f"{time_origin!r}, not 'infection'.",
                    UserWarning,
                    stacklevel=2,
                )
            else:
                warnings.warn(
                    f"{event!r} is an administrative reference event, which has "
                    "no fixed offset from infection -- it reflects testing "
                    "behaviour and health-system access. time_origin is recorded "
                    f"as {time_origin!r}, not 'infection'.",
                    UserWarning,
                    stacklevel=2,
                )

    frame.attrs = {
        "time_origin": time_origin,
        "reference_event_class": event_class,
        "incubation_applied": incubation_applied,
        "model": source.model,
        "unit": source.unit,
        "biomarker": getattr(source, "biomarker", None),
        "specimen": getattr(source, "specimen", None),
    }
    return frame


def plot_simulated_shedding(
    traj: pd.DataFrame,
    *,
    source=None,
    observed: dict | None = None,
    quantiles: tuple[float, float, float] = (0.05, 0.5, 0.95),
    figsize: tuple[float, float] = (8, 6),
) -> Figure:
    """
    Plot the median and a credible band of simulated trajectories.

    Args:
        traj: Output of ``simulate_shedding``.
        source: Optional fit or ensemble, used to draw the censoring limit and,
            when ``observed`` is given, to determine which analyte(s) to
            overlay.
        observed: Optional dataset dictionary; its measurements are overlaid as
            points so simulated and real trajectories can be compared, filtered
            to the analyte(s) contributed by ``source`` (a ``SheddingFit``'s
            own analyte, or a ``SheddingEnsemble``'s component analytes).
            Requires ``source``.
        quantiles: Lower, middle, and upper quantiles for the band.
        figsize: Figure size in inches.

    Returns:
        The figure. It is closed in the pyplot state so notebooks do not display
        it twice, matching the convention in ``shedding_peak.py``.
    """
    if traj.empty:
        raise ValueError("Simulation result is empty, cannot create plot")
    if traj["log10_value"].isna().all():
        raise ValueError(
            "Every simulated value is NaN, so there is nothing to plot. Under "
            "the gamma model this happens when every requested time falls at "
            "or before the reference event (or, with incubation_period set, "
            "within the incubation window). Request later times or reduce "
            "incubation_period."
        )

    lower, middle, upper = quantiles
    summary = (
        traj.dropna(subset=["log10_value"])
        .groupby("time")["log10_value"]
        .quantile([lower, middle, upper])
        .unstack()
    )

    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(
        summary.index,
        summary[lower],
        summary[upper],
        alpha=0.25,
        color="tab:blue",
        label=f"{int((upper - lower) * 100)}% of simulated individuals",
    )
    ax.plot(summary.index, summary[middle], color="tab:blue", lw=2, label="Median")

    if source is not None:
        ax.axhline(
            source.censoring_limit,
            ls=":",
            color="gray",
            label="Limit of quantification",
        )

    if observed is not None:
        if source is None:
            raise ValueError(
                "observed requires a source (a SheddingFit or SheddingEnsemble) "
                "to identify which analyte(s) to overlay; pass source=... or "
                "drop observed."
            )
        # Duck-typed rather than isinstance-checked against SheddingEnsemble, to
        # avoid this plotting module importing shedding_ensemble.py. A
        # SheddingFit contributes its own analyte; a SheddingEnsemble (which has
        # no analyte of its own) contributes the analytes of its component fits.
        if hasattr(source, "fits"):
            analytes = {fit.analyte for fit in source.fits}
        else:
            analytes = {source.analyte}
        times, values = [], []
        for participant in observed.get("participants", []):
            for measurement in participant.get("measurements") or []:
                if measurement.get("analyte") not in analytes:
                    continue
                time = measurement.get("time")
                value = measurement.get("value")
                if isinstance(time, (int, float)) and isinstance(value, (int, float)):
                    times.append(float(time))
                    values.append(np.log10(float(value)))
        if times:
            ax.scatter(times, values, s=18, color="black", alpha=0.5, label="Observed")

    origin = traj.attrs.get("time_origin", "reference event")
    unit = traj.attrs.get("unit", "")
    ax.set_xlabel(f"Days after {origin}")
    ax.set_ylabel(f"log10 concentration ({unit})" if unit else "log10 concentration")
    ax.set_title("Simulated shedding trajectories")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.close(fig)
    return fig
