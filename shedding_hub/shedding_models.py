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
"""

import numpy as np

MODELS = ("exponential", "gamma")

PARAM_NAMES = {
    "exponential": ("a0", "c0"),
    "gamma": ("a0", "b0", "c0"),
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
    return (
        params[:, 2] + params[:, 1] * _safe_log(times) - params[:, 0] * times
    ) / LN10


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
    return params[:, 1] / params[:, 0]


def half_life_days(model: str, params: np.ndarray) -> np.ndarray:
    """
    Half-life of the late-phase decline, ``ln(2) / a0``.

    Exact for the exponential model; asymptotic for the gamma model, whose
    decline approaches rate ``a0`` once ``t`` is well past the peak.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    return np.log(2.0) / params[:, 0]
