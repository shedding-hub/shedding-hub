import numpy as np
import pytest

from shedding_hub.shedding_models import (
    LN10,
    MODELS,
    PARAM_NAMES,
    half_life_days,
    log10_concentration,
    log10_concentration_pointwise,
    log10_concentration_rowwise,
    from_population_coords,
    peak_day,
    to_population_coords,
    validate_model,
)


def test_models_and_param_names():
    assert MODELS == ("exponential", "gamma")
    assert PARAM_NAMES["exponential"] == ("a0", "c0")
    assert PARAM_NAMES["gamma"] == ("a0", "b0", "c0")


def test_exponential_matches_closed_form():
    params = np.array([[0.5, 20.0]])  # a0, c0
    times = np.array([0.0, 2.0, 10.0])
    got = log10_concentration("exponential", params, times)
    expected = (20.0 - 0.5 * times) / LN10
    np.testing.assert_allclose(got[0], expected)


def test_gamma_matches_closed_form():
    params = np.array([[0.5, 2.0, 10.0]])  # a0, b0, c0
    times = np.array([1.0, 4.0])
    got = log10_concentration("gamma", params, times)
    expected = (10.0 + 2.0 * np.log(times) - 0.5 * times) / LN10
    np.testing.assert_allclose(got[0], expected)


def test_gamma_is_nan_for_non_positive_times():
    params = np.array([[0.5, 2.0, 10.0]])
    got = log10_concentration("gamma", params, np.array([-1.0, 0.0, 1.0]))
    assert np.isnan(got[0, 0])
    assert np.isnan(got[0, 1])
    assert np.isfinite(got[0, 2])


def test_exponential_is_finite_for_negative_times():
    params = np.array([[0.5, 20.0]])
    got = log10_concentration("exponential", params, np.array([-3.0]))
    assert np.isfinite(got[0, 0])


def test_gamma_peaks_at_b_over_a():
    params = np.array([[0.5, 2.0, 10.0]])
    assert peak_day("gamma", params)[0] == pytest.approx(4.0)
    # confirm the curve really is maximal there
    times = np.linspace(0.1, 20, 500)
    values = log10_concentration("gamma", params, times)[0]
    assert times[np.argmax(values)] == pytest.approx(4.0, abs=0.1)


def test_exponential_peak_day_is_zero():
    params = np.array([[0.5, 20.0]])
    assert peak_day("exponential", params)[0] == 0.0


def test_half_life():
    params = np.array([[np.log(2.0), 20.0]])
    assert half_life_days("exponential", params)[0] == pytest.approx(1.0)


def test_rowwise_gives_each_individual_its_own_times():
    params = np.array([[0.5, 20.0], [1.0, 20.0]])
    times = np.array([[0.0, 1.0], [0.0, 1.0]])
    got = log10_concentration_rowwise("exponential", params, times)
    assert got.shape == (2, 2)
    np.testing.assert_allclose(got[0, 1], (20.0 - 0.5) / LN10)
    np.testing.assert_allclose(got[1, 1], (20.0 - 1.0) / LN10)


def test_pointwise_pairs_each_observation_with_its_params():
    params = np.array([[0.5, 20.0], [1.0, 20.0]])
    times = np.array([2.0, 2.0])
    got = log10_concentration_pointwise("exponential", params, times)
    assert got.shape == (2,)
    np.testing.assert_allclose(got[0], (20.0 - 1.0) / LN10)
    np.testing.assert_allclose(got[1], (20.0 - 2.0) / LN10)


def test_unknown_model_raises():
    with pytest.raises(ValueError, match="Unknown model"):
        validate_model("weibull")


# --- population coordinates -------------------------------------------------
#
# The population summary averages subjects in a coordinate system where a
# Gaussian is defensible. For the gamma model the natural log-parameters are
# not such a system: c0 is the concentration at t = 1, so its meaning depends
# on b0, and the two trade off along a curved ridge.


def test_population_coords_round_trip_exactly():
    params = np.array([[0.5, 2.0, 12.0], [1.81, 13.74, 0.04], [0.24, 0.03, 9.1]])
    coords = to_population_coords("gamma", params)
    np.testing.assert_allclose(from_population_coords("gamma", coords), params)


def test_gamma_population_coords_are_decay_peak_day_and_peak_height():
    params = np.array([[0.5, 2.0, 12.0]])
    a0, b0, c0 = params[0]
    coords = to_population_coords("gamma", params)[0]
    peak = b0 / a0
    np.testing.assert_allclose(coords[0], np.log(a0))
    np.testing.assert_allclose(coords[1], np.log(peak))
    # the third coordinate is the log10 concentration the curve reaches at peak
    np.testing.assert_allclose(
        coords[2], log10_concentration("gamma", params, np.array([peak]))[0, 0]
    )


def test_exponential_population_coords_are_decay_and_peak_height():
    """The exponential model peaks at t = 0, where log10 c = c0 / ln(10)."""
    params = np.array([[0.6, 18.0]])
    a0, _ = params[0]
    coords = to_population_coords("exponential", params)[0]
    np.testing.assert_allclose(coords[0], np.log(a0))
    np.testing.assert_allclose(
        coords[1], log10_concentration("exponential", params, np.array([0.0]))[0, 0]
    )


def test_exponential_population_coords_round_trip_exactly():
    params = np.array([[0.6, 18.0], [0.24, 9.1], [1.9, 22.5]])
    coords = to_population_coords("exponential", params)
    np.testing.assert_allclose(from_population_coords("exponential", coords), params)


def test_exponential_height_coordinate_is_linear_in_log10_concentration():
    """The property that keeps simulated tails finite.

    Modelling ``log c0`` as normal would make the *log10* concentration
    lognormal, hence the concentration a double exponential of the draw, so a
    3.5-sigma individual reached 10**18.8 where the median reached 10**6.1.
    Making the log10 height itself a coordinate means a draw k units above the
    mean lands exactly k log10 above it.
    """
    for height in (6.0, 9.0, 12.0):
        params = from_population_coords(
            "exponential", np.array([[np.log(0.6), height]])
        )
        reached = log10_concentration("exponential", params, np.array([0.0]))[0, 0]
        np.testing.assert_allclose(reached, height)


def test_population_coords_reject_unknown_model():
    with pytest.raises(ValueError, match="Unknown model"):
        to_population_coords("weibull", np.array([[1.0, 1.0]]))
