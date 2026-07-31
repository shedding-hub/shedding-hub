"""
Parametric shedding-curve models.

Two models are supported, both taken from the Shedding Hub Rstan tutorial and
expressed on the log10 scale, which is the scale likelihoods are evaluated on:

- ``exponential``: ``c(t) = c0 * exp(-a0 * t)``, a pure decay appropriate when
  sampling begins at or after peak shedding.
- ``gamma``: ``c(t) = c0 * t**b0 * exp(-a0 * t)``, which rises and falls, peaking
  at ``t = b0 / a0``.

``c0`` is on the natural-log scale, so the log10 concentration at ``t = 0`` is
``c0 / ln(10)``. All parameters are strictly positive.

This module is pure math: no dataset handling, no I/O.

Example:
    >>> import shedding_hub as sh
    >>> sh.MODELS
    ('exponential', 'gamma', 'gamma_shifted')
    >>> sh.PARAM_NAMES['gamma']
    ('a0', 'b0', 'c0')
"""

import numpy as np

MODELS = ("exponential", "gamma", "gamma_shifted")

PARAM_NAMES = {
    "exponential": ("a0", "c0"),
    "gamma": ("a0", "b0", "c0"),
    # t0 is the onset of shedding: the gamma model's support starts there rather
    # than at the reference event, so readings before the reference event become
    # evaluable and the onset itself is estimated instead of assumed.
    "gamma_shifted": ("a0", "b0", "c0", "t0"),
}

# The coordinates each model's population summary is averaged in; see
# ``to_population_coords`` for why they differ between the two models. Named
# rather than merely implied, so that a serialized fit records which space its
# ``population_mean``/``population_cov`` live in and a catalog written under the
# old convention fails loudly instead of being silently misread.
POPULATION_COORDS = {
    "exponential": ("log_a0", "peak_log10"),
    "gamma": ("log_a0", "log_peak_day", "peak_log10"),
    # log_rise_days is the interval from onset to peak, b0/a0, which is the
    # shape quantity; the absolute peak is t0 + b0/a0. t0 rides along
    # untransformed, being a time on the whole real line rather than a positive
    # scale.
    "gamma_shifted": ("log_a0", "log_rise_days", "peak_log10", "t0"),
}

LN10 = float(np.log(10.0))


def validate_model(model: str) -> None:
    """Raise ``ValueError`` unless ``model`` is a supported model name."""
    if model not in MODELS:
        raise ValueError(f"Unknown model {model!r}. Choose one of {list(MODELS)}.")


def _safe_log(times: np.ndarray) -> np.ndarray:
    """Natural log of ``times``, NaN where ``times <= 0``."""
    positive = times > 0
    return np.where(positive, np.log(np.where(positive, times, 1.0)), np.nan)


def log10_concentration(
    model: str, params: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """
    Evaluate the model for every combination of individual and time.

    Args:
        model: ``"exponential"`` or ``"gamma"``.
        params: Natural-scale parameters, shape ``(n, k)``, ordered as
            ``PARAM_NAMES[model]``.
        times: Times since the reference event, shape ``(m,)``.

    Returns:
        Log10 concentrations, shape ``(n, m)``. Under the gamma model,
        non-positive times yield NaN because ``ln(t)`` is undefined there.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    times = np.asarray(times, dtype=float)
    return log10_concentration_rowwise(
        model, params, np.broadcast_to(times, (params.shape[0], times.size))
    )


def log10_concentration_rowwise(
    model: str, params: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """
    Evaluate the model giving each individual its own row of times.

    Args:
        model: ``"exponential"`` or ``"gamma"``.
        params: Natural-scale parameters, shape ``(n, k)``.
        times: Times, shape ``(n, m)`` — row ``i`` belongs to individual ``i``.

    Returns:
        Log10 concentrations, shape ``(n, m)``.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    times = np.asarray(times, dtype=float)
    if model == "exponential":
        a0 = params[:, 0:1]
        c0 = params[:, 1:2]
        return (c0 - a0 * times) / LN10
    a0 = params[:, 0:1]
    b0 = params[:, 1:2]
    c0 = params[:, 2:3]
    if model == "gamma_shifted":
        times = times - params[:, 3:4]
    return (c0 + b0 * _safe_log(times) - a0 * times) / LN10


def log10_concentration_pointwise(
    model: str, params: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """
    Evaluate the model once per observation.

    Args:
        model: ``"exponential"`` or ``"gamma"``.
        params: Natural-scale parameters, shape ``(n_obs, k)`` — row ``j`` holds
            the parameters of the subject that observation ``j`` belongs to.
        times: Observation times, shape ``(n_obs,)``.

    Returns:
        Log10 concentrations, shape ``(n_obs,)``.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    times = np.asarray(times, dtype=float)
    if model == "exponential":
        return (params[:, 1] - params[:, 0] * times) / LN10
    if model == "gamma_shifted":
        times = times - params[:, 3]
    return (
        params[:, 2] + params[:, 1] * _safe_log(times) - params[:, 0] * times
    ) / LN10


def theta_to_params(model: str, theta: np.ndarray) -> np.ndarray:
    """
    Map the optimizer's unconstrained coordinates to natural parameters.

    The fitter works in ``theta`` so that positivity is automatic: every
    parameter is a positive scale, so ``exp`` both enforces the constraint and
    puts the optimizer on a well-behaved scale.

    ``gamma_shifted``'s ``t0`` is the exception. It is a *time*, on the whole
    real line, and exponentiating it would forbid an onset before the reference
    event — the one thing the model exists to allow. It is therefore carried
    untransformed, and its constraint (``t0`` below every reading it must
    explain) is imposed by the optimizer's bounds instead.
    """
    validate_model(model)
    theta = np.atleast_2d(np.asarray(theta, dtype=float))
    if model == "gamma_shifted":
        return np.column_stack([np.exp(theta[:, :3]), theta[:, 3]])
    return np.exp(theta)


def params_to_theta(model: str, params: np.ndarray) -> np.ndarray:
    """Invert ``theta_to_params``."""
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    if model == "gamma_shifted":
        return np.column_stack([np.log(params[:, :3]), params[:, 3]])
    return np.log(params)


def to_population_coords(model: str, params: np.ndarray) -> np.ndarray:
    """
    Map natural-scale parameters into the space the population is summarized in.

    A population summary averages subjects and treats the result as a Gaussian.
    That is only defensible in coordinates where the subjects actually form a
    compact, roughly-elliptical cloud, and for the gamma model the natural
    log-parameters are not such coordinates: ``c0`` is the concentration at
    ``t = 1``, so what counts as a plausible value depends entirely on ``b0``.
    Real subjects therefore lie along a *curved* ridge — in
    ``woelfel2020virological`` stool, ``b0`` spans 0.03 to 13.7 while ``c0``
    counter-varies from 19.7 down to 0.04 — and the coordinate-wise mean of
    their logs lands off that ridge, describing a curve no subject resembles.
    Measured across the repository's gamma fits, that mean sat a median of 2.14
    log10 away from the subjects' own median curve, and simulating from it
    produced 95th-percentile concentrations above 10^20 gc/mL in 16 of 23 fits.

    The gamma model is therefore summarized as ``(log a0, log t_peak, y_peak)``:
    the decay rate, the day of peak shedding, and the log10 concentration
    reached at that peak. These are the quantities the model is actually
    identified in — each is separately interpretable, and they are close to
    uncorrelated in practice (``corr(log t_peak, y_peak)`` is 0.07 for woelfel
    stool against -0.55 for ``corr(log b0, log c0)``), so averaging them
    coordinate-wise is meaningful.

    The exponential model takes the same treatment for a different reason. Its
    ``c0`` *is* the level, with no shape parameter to trade against, so its
    subjects do form a compact cloud in ``log(params)`` and the *ridge* argument
    above does not apply to it. What does apply is the tail. Because
    ``log10 c(0) = c0 / ln(10)``, modelling ``log c0`` as normal makes the log10
    concentration itself lognormal, hence the concentration a *double*
    exponential of the draw: on ``woelfel2020virological`` stool that put the
    99.9th percentile of simulated day-0 concentrations at ``10**16.6`` and the
    worst draw at ``10**25``, against eight real subjects topping out at
    ``10**9``, and left the top 0.1% of a simulated cohort carrying
    essentially all of its shed load. Taking ``peak_log10`` — the model peaks at
    ``t = 0`` — as the coordinate makes a draw ``k`` units above the mean land
    exactly ``k`` log10 above it, so concentration is merely lognormal, as it is
    for the gamma model.

    This is a genuine trade, and it is why the exponential model was originally
    left in ``log(params)``. The level summary moves from the geometric mean of
    ``c0`` to the arithmetic mean of the subjects' log10 heights, and measured
    against the subjects' own pointwise median curve over a 10-fit sample, the
    new coordinate was closer on only 3 of 10 (median RMS gap 1.10 against 1.16
    log10, mean 2.36 against 1.85). It is accepted because the two costs are not
    remotely the same size: the median individual moves by a fraction of a log10,
    while the simulated maximum falls from ``10**27`` to ``10**15``. Simulation
    is what the catalog is for.

    Args:
        model: ``"exponential"`` or ``"gamma"``.
        params: Natural-scale parameters, shape ``(n, k)``, ordered as
            ``PARAM_NAMES[model]``.

    Returns:
        Population coordinates, shape ``(n, k)``. Exactly invertible by
        ``from_population_coords``.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    if model == "exponential":
        # The peak is at t = 0, where the log10 concentration is c0 / ln(10).
        return np.column_stack([np.log(params[:, 0]), params[:, 1] / LN10])
    a0 = params[:, 0]
    b0 = params[:, 1]
    c0 = params[:, 2]
    rise = b0 / a0
    # At the peak a0 * rise == b0, so the log10 height collapses to this. The
    # shift cancels: the height depends on the interval from onset, not on when
    # the onset was.
    height = (c0 + b0 * np.log(rise) - b0) / LN10
    if model == "gamma_shifted":
        return np.column_stack([np.log(a0), np.log(rise), height, params[:, 3]])
    return np.column_stack([np.log(a0), np.log(rise), height])


def from_population_coords(model: str, coords: np.ndarray) -> np.ndarray:
    """
    Invert ``to_population_coords``, back to natural-scale parameters.

    Args:
        model: ``"exponential"`` or ``"gamma"``.
        coords: Population coordinates, shape ``(n, k)``.

    Returns:
        Natural-scale parameters, shape ``(n, k)``, ordered as
        ``PARAM_NAMES[model]``.

        ``a0`` and ``b0`` are always strictly positive, being exponentials of
        the coordinates. ``c0`` may come back non-positive, because ``y_peak``
        is a log10 concentration and is modelled on the whole real line: a
        drawn individual whose peak sits below 1 gc/mL is meaningful, and every
        function in this module evaluates such a curve correctly. Only the
        *fitted* parameters are constrained positive, by the optimizer.
    """
    validate_model(model)
    coords = np.atleast_2d(np.asarray(coords, dtype=float))
    if model == "exponential":
        return np.column_stack([np.exp(coords[:, 0]), coords[:, 1] * LN10])
    a0 = np.exp(coords[:, 0])
    rise = np.exp(coords[:, 1])
    height = coords[:, 2]
    b0 = a0 * rise
    c0 = LN10 * height - b0 * np.log(rise) + b0
    if model == "gamma_shifted":
        return np.column_stack([a0, b0, c0, coords[:, 3]])
    return np.column_stack([a0, b0, c0])


def peak_day(model: str, params: np.ndarray) -> np.ndarray:
    """
    Time of peak shedding, in days after the reference event.

    The gamma model peaks at ``b0 / a0``. The exponential model is monotonically
    decreasing, so its maximum is at the reference event and this returns 0.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    if model == "exponential":
        return np.zeros(params.shape[0])
    peak = params[:, 1] / params[:, 0]
    if model == "gamma_shifted":
        # Absolute, so it stays comparable with the other models' peak times.
        return params[:, 3] + peak
    return peak


def half_life_days(model: str, params: np.ndarray) -> np.ndarray:
    """
    Half-life of the late-phase decline, ``ln(2) / a0``.

    Exact for the exponential model; asymptotic for the gamma model, whose
    decline approaches rate ``a0`` once ``t`` is well past the peak.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    return np.log(2.0) / params[:, 0]
