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
    peak_day,
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
