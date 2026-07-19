import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from shedding_hub.shedding_fit import fit_shedding_model
from shedding_hub.shedding_simulate import simulate_shedding


@pytest.fixture
def exponential_fit(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=25
    )
    return fit_shedding_model(dataset, analyte="stool", model="exponential")


@pytest.fixture
def gamma_fit(make_synthetic_dataset):
    mu = np.array([np.log(0.5), np.log(2.0), np.log(12.0)])
    dataset = make_synthetic_dataset(
        "gamma", mu, np.diag([0.04, 0.04, 0.04]), n_subjects=25
    )
    return fit_shedding_model(dataset, analyte="stool", model="gamma")


def test_returns_tidy_frame_of_expected_shape(exponential_fit):
    times = np.arange(0.0, 10.0)
    traj = simulate_shedding(exponential_fit, n_individuals=20, times=times, seed=0)
    assert isinstance(traj, pd.DataFrame)
    assert list(traj.columns) == [
        "individual_id",
        "time",
        "log10_value",
        "value",
        "detected",
        "source_dataset_id",
    ]
    assert len(traj) == 20 * len(times)
    assert traj["individual_id"].nunique() == 20


def test_is_reproducible_with_a_seed(exponential_fit):
    a = simulate_shedding(exponential_fit, n_individuals=10, times=[1, 2], seed=7)
    b = simulate_shedding(exponential_fit, n_individuals=10, times=[1, 2], seed=7)
    c = simulate_shedding(exponential_fit, n_individuals=10, times=[1, 2], seed=8)
    pd.testing.assert_frame_equal(a, b)
    assert not np.allclose(a["log10_value"], c["log10_value"])


def test_measurement_error_is_off_by_default(exponential_fit):
    clean = simulate_shedding(exponential_fit, n_individuals=200, times=[3.0], seed=3)
    noisy = simulate_shedding(
        exponential_fit,
        n_individuals=200,
        times=[3.0],
        seed=3,
        include_measurement_error=True,
    )
    assert noisy["log10_value"].std() > clean["log10_value"].std()


def test_detected_reflects_the_censoring_limit(exponential_fit):
    traj = simulate_shedding(
        exponential_fit, n_individuals=50, times=np.arange(0.0, 40.0), seed=1
    )
    finite = traj.dropna(subset=["log10_value"])
    expected = finite["log10_value"] >= exponential_fit.censoring_limit
    pd.testing.assert_series_equal(finite["detected"], expected, check_names=False)


def test_values_below_the_limit_are_not_clipped(exponential_fit):
    traj = simulate_shedding(
        exponential_fit, n_individuals=50, times=np.arange(0.0, 60.0), seed=2
    )
    below = traj[~traj["detected"]].dropna(subset=["log10_value"])
    assert (below["log10_value"] < exponential_fit.censoring_limit).all()
    assert below["log10_value"].min() < exponential_fit.censoring_limit - 1


def test_gamma_non_positive_times_are_nan_and_undetected(gamma_fit):
    traj = simulate_shedding(gamma_fit, n_individuals=5, times=[-1.0, 0.0, 5.0], seed=0)
    non_positive = traj[traj["time"] <= 0]
    assert non_positive["log10_value"].isna().all()
    assert not non_positive["detected"].any()


def test_scalar_incubation_shifts_the_curve(exponential_fit):
    unshifted = simulate_shedding(
        exponential_fit, n_individuals=30, times=[0.0], seed=5
    )
    shifted = simulate_shedding(
        exponential_fit,
        n_individuals=30,
        times=[5.0],
        incubation_period=5.0,
        seed=5,
    )
    np.testing.assert_allclose(
        unshifted["log10_value"].to_numpy(),
        shifted["log10_value"].to_numpy(),
    )


def test_array_incubation_must_match_n_individuals(exponential_fit):
    with pytest.raises(ValueError, match="incubation_period"):
        simulate_shedding(
            exponential_fit,
            n_individuals=5,
            times=[1.0],
            incubation_period=np.array([1.0, 2.0]),
        )


def test_callable_incubation_is_drawn_per_individual(exponential_fit):
    def incubation(rng, n):
        return rng.uniform(2.0, 8.0, size=n)

    traj = simulate_shedding(
        exponential_fit,
        n_individuals=40,
        times=[6.0],
        incubation_period=incubation,
        seed=4,
    )
    # Different offsets mean different effective times, hence spread in values.
    assert traj["log10_value"].std() > 0


def test_result_attrs_record_the_time_origin(exponential_fit):
    native = simulate_shedding(exponential_fit, n_individuals=3, times=[1.0], seed=0)
    assert native.attrs["time_origin"] == "symptom onset"
    assert native.attrs["incubation_applied"] is False

    infection = simulate_shedding(
        exponential_fit,
        n_individuals=3,
        times=[1.0],
        incubation_period=5.0,
        seed=0,
    )
    assert infection.attrs["time_origin"] == "infection"
    assert infection.attrs["incubation_applied"] is True
    assert infection.attrs["unit"] == "gc/mL"


def test_n_individuals_must_be_positive(exponential_fit):
    with pytest.raises(ValueError):
        simulate_shedding(exponential_fit, n_individuals=0, times=[1.0])
