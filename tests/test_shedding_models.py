import numpy as np
import pytest

from shedding_hub.shedding_models import (
    LN10,
    MODELS,
    PARAM_NAMES,
    POPULATION_COORDS,
    half_life_days,
    log10_concentration,
    log10_concentration_pointwise,
    log10_concentration_rowwise,
    from_population_coords,
    peak_day,
    population_coord_names,
    to_population_coords,
    validate_model,
)


def test_models_and_param_names():
    assert MODELS == ("exponential", "gamma", "gamma_shifted")
    assert PARAM_NAMES["exponential"] == ("a0", "c0")
    assert PARAM_NAMES["gamma"] == ("a0", "b0", "c0")
    assert PARAM_NAMES["gamma_shifted"] == ("a0", "b0", "c0", "t0")


# ---------------------------------------------------------------------------
# gamma_shifted: c(t) = c0 * (t - t0)**b0 * exp(-a0 * (t - t0)),  t > t0
#
# The gamma model's support starts at the reference event, so every reading at
# t <= 0 is discarded -- 1,026 of them across 14 catalog fits, and 29% of
# hakki2022onset's. Freeing the onset lets those readings inform the fit, and
# makes the onset itself estimable rather than assumed, which matters because
# the catalog's five reference events do not mean the same thing.
# ---------------------------------------------------------------------------


def test_theta_maps_positive_parameters_through_exp():
    from shedding_hub.shedding_models import params_to_theta, theta_to_params

    for model, params in [
        ("exponential", np.array([[0.6, 18.0]])),
        ("gamma", np.array([[0.5, 2.0, 12.0]])),
    ]:
        theta = params_to_theta(model, params)
        np.testing.assert_allclose(theta, np.log(params))
        np.testing.assert_allclose(theta_to_params(model, theta), params)


def test_theta_leaves_the_onset_untransformed():
    """t0 is a time, not a positive scale: exponentiating it would forbid an
    onset before the reference event, which is the whole purpose of the model."""
    from shedding_hub.shedding_models import params_to_theta, theta_to_params

    params = np.array([[0.5, 2.0, 12.0, -4.0], [1.0, 3.0, 9.0, 0.5]])
    theta = params_to_theta("gamma_shifted", params)
    np.testing.assert_allclose(theta[:, :3], np.log(params[:, :3]))
    np.testing.assert_allclose(theta[:, 3], params[:, 3])
    np.testing.assert_allclose(theta_to_params("gamma_shifted", theta), params)


def test_gamma_shifted_reduces_to_gamma_when_the_shift_is_zero():
    times = np.array([0.5, 1.0, 4.0, 9.0])
    plain = log10_concentration("gamma", np.array([[0.5, 2.0, 12.0]]), times)
    shifted = log10_concentration(
        "gamma_shifted", np.array([[0.5, 2.0, 12.0, 0.0]]), times
    )
    np.testing.assert_allclose(shifted, plain)


def test_gamma_shifted_uses_times_before_the_reference_event():
    """The whole point: readings at t <= 0 become evaluable once t0 < 0."""
    params = np.array([[0.5, 2.0, 12.0, -4.0]])
    values = log10_concentration("gamma_shifted", params, np.array([-3.0, -1.0, 2.0]))
    assert np.isfinite(values).all()


def test_gamma_shifted_is_undefined_at_or_before_its_own_onset():
    params = np.array([[0.5, 2.0, 12.0, -4.0]])
    values = log10_concentration("gamma_shifted", params, np.array([-5.0, -4.0, -3.9]))
    assert np.isnan(values[0, 0]) and np.isnan(values[0, 1])
    assert np.isfinite(values[0, 2])


def test_gamma_shifted_peaks_a_rise_after_its_onset():
    a0, b0, c0, t0 = 0.5, 2.0, 12.0, -4.0
    params = np.array([[a0, b0, c0, t0]])
    expected = t0 + b0 / a0
    np.testing.assert_allclose(peak_day("gamma_shifted", params)[0], expected)

    grid = np.linspace(t0 + 1e-3, t0 + 30.0, 4000)
    values = log10_concentration("gamma_shifted", params, grid)[0]
    np.testing.assert_allclose(grid[values.argmax()], expected, atol=0.02)


def test_gamma_shifted_half_life_is_the_same_asymptotic_decay():
    params = np.array([[0.5, 2.0, 12.0, -4.0]])
    np.testing.assert_allclose(
        half_life_days("gamma_shifted", params)[0], np.log(2) / 0.5
    )


def test_gamma_shifted_population_coords_round_trip_exactly():
    params = np.array(
        [[0.5, 2.0, 12.0, -4.0], [1.2, 0.3, 6.0, 1.5], [0.2, 5.0, 3.0, 0.0]]
    )
    coords = to_population_coords("gamma_shifted", params)
    np.testing.assert_allclose(from_population_coords("gamma_shifted", coords), params)


def test_gamma_shifted_height_coordinate_is_the_log10_at_its_peak():
    params = np.array([[0.5, 2.0, 12.0, -4.0]])
    coords = to_population_coords("gamma_shifted", params)[0]
    at_peak = log10_concentration(
        "gamma_shifted", params, np.array([peak_day("gamma_shifted", params)[0]])
    )[0, 0]
    np.testing.assert_allclose(coords[2], at_peak)
    # the onset rides along untransformed, being a time on the whole real line
    np.testing.assert_allclose(coords[3], -4.0)


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


def test_population_coord_names_default_to_concentration():
    assert population_coord_names("gamma") == POPULATION_COORDS["gamma"]


def test_population_coord_names_rename_the_height_for_ct():
    assert population_coord_names("gamma", "ct") == (
        "log_a0",
        "log_peak_day",
        "peak_cycles",
    )


def test_population_coord_names_leave_temporal_coordinates_alone():
    # t0 is a time on either scale and must not be renamed.
    assert population_coord_names("gamma_shifted", "ct")[-1] == "t0"
