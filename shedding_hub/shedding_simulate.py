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


# How far below zero the y axis will follow a band that sets it. A simulated
# range descends past 10**-100 gc/mL -- measured at -137 log10 on a
# gamma_shifted cohort from the shipped catalog -- which is not a concentration
# anyone needs plotted, and following it leaves the data a ribbon at the top of
# the panel. -3 sits below every censoring limit in the repository (the lowest
# is -2.37) while keeping the data legible. It bounds the band and never the
# observations: the axis takes the lower of this and the data. Same value and
# same rule as ``viz.FIT_DIAGNOSTIC_YLIM_FLOOR``, so the two agree.
SIMULATION_YLIM_FLOOR = -3.0


def plot_simulated_shedding(
    traj: pd.DataFrame,
    *,
    source=None,
    observed: dict | None = None,
    band_quantiles: tuple[float, float] = (0.0, 1.0),
    band_inner_quantiles: tuple[float, float] | None = (0.025, 0.975),
    ylim_floor: float = SIMULATION_YLIM_FLOOR,
    figsize: tuple[float, float] = (8, 6),
) -> Figure:
    """
    Plot the median, an outer band, and an inner interval of a simulated cohort.

    The layout follows the review pages ``make review_range`` produces: the wide
    interval is shaded, the narrow one is drawn inside it as a dashed pair, and
    the median is solid. Shading the range shows what the fit considers
    possible; the dashed interval shows where the mass actually sits, which the
    range alone hides.

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
        band_quantiles: Lower and upper quantiles of the shaded band. Defaults
            to ``(0.0, 1.0)``, the full simulated range.
        band_inner_quantiles: Lower and upper quantiles of the dashed interval
            drawn inside the band, or ``None`` to omit it. Defaults to the
            central 95%.
        ylim_floor: How far below zero the axis will follow the band. Never
            clips the observations.
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

    defined = traj.dropna(subset=["log10_value"])
    lower, upper = band_quantiles
    summary = (
        defined.groupby("time")["log10_value"].quantile([lower, 0.5, upper]).unstack()
    )

    fig, ax = plt.subplots(figsize=figsize)

    if band_inner_quantiles is not None:
        # Drawn before the fill so the shading does not sit on top of it. Only
        # the first edge carries a label, so the legend gets one entry for the
        # pair rather than two identical ones.
        inner = (
            defined.groupby("time")["log10_value"]
            .quantile(list(band_inner_quantiles))
            .unstack()
        )
        width = int(round((band_inner_quantiles[1] - band_inner_quantiles[0]) * 100))
        for index, quantile in enumerate(band_inner_quantiles):
            ax.plot(
                inner.index,
                inner[quantile],
                ls="--",
                color="tab:blue",
                lw=1.1,
                alpha=0.75,
                zorder=3,
                label=f"{width}% of individuals" if index == 0 else "_nolegend_",
            )

    ax.fill_between(
        summary.index,
        summary[lower],
        summary[upper],
        alpha=0.18,
        color="tab:blue",
        lw=0,
        zorder=1,
        # A full-range band is labelled by its draw count, not "100% of
        # individuals": how far the extremes reach is a property of how many
        # individuals were drawn as much as of the population, and the figure
        # should not imply otherwise.
        label=(
            f"Simulated range, {traj['individual_id'].nunique()} individuals"
            if lower <= 0.0 and upper >= 1.0
            else f"Simulated {int(round((upper - lower) * 100))}% of individuals"
        ),
    )
    ax.plot(
        summary.index, summary[0.5], color="tab:blue", lw=2, zorder=4, label="Median"
    )

    # The values the axis must keep whatever the band does: real readings and
    # the censoring limit, never the simulated band itself.
    protected_low = None
    protected_high = None

    if source is not None:
        ax.axhline(
            source.censoring_limit,
            ls=":",
            color="gray",
            label="Limit of quantification",
        )
        protected_low = float(source.censoring_limit)

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
            finite = [value for value in values if np.isfinite(value)]
            if finite:
                seen_low, seen_high = min(finite), max(finite)
                protected_low = (
                    seen_low if protected_low is None else min(protected_low, seen_low)
                )
                protected_high = (
                    seen_high
                    if protected_high is None
                    else max(protected_high, seen_high)
                )

    # The axis follows the inner interval, not the band. A range over a few
    # hundred draws reaches both directions of nonsense -- 10**-137 and 10**76
    # gc/mL were both measured on the shipped catalog -- and letting it set the
    # limits hands the panel to one agent and leaves everything real a ribbon in
    # the middle. The band is still drawn; it is simply allowed to clip.
    focus = inner if band_inner_quantiles is not None else summary[[lower, upper]]
    low, high = float(np.nanmin(focus.values)), float(np.nanmax(focus.values))
    pad = 0.05 * (high - low) if high > low else 1.0
    bottom, top = low - pad, high + pad

    # Anything the study actually recorded, and the limit it was read against,
    # outrank that: the axis bounds the simulation, never the data.
    if protected_low is not None:
        bottom = min(bottom, protected_low - pad)
    if protected_high is not None:
        top = max(top, protected_high + pad)
    if protected_low is None or protected_low >= ylim_floor:
        bottom = max(bottom, ylim_floor)
    ax.set_ylim(bottom, top)

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
