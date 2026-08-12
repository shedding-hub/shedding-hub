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


from matplotlib.figure import Figure

from shedding_hub.shedding_simulate import plot_simulated_shedding


def test_plot_returns_a_figure(exponential_fit):
    traj = simulate_shedding(
        exponential_fit, n_individuals=40, times=np.arange(0.0, 20.0), seed=0
    )
    fig = plot_simulated_shedding(traj, source=exponential_fit)
    assert isinstance(fig, Figure)


def test_plot_rejects_an_empty_frame():
    with pytest.raises(ValueError, match="empty"):
        plot_simulated_shedding(pd.DataFrame())


def test_plot_rejects_an_all_nan_trajectory(gamma_fit):
    """A non-empty but all-NaN trajectory must fail with a helpful message.

    Reproduction from the review: requesting only times at or before the
    reference event under the gamma model, further shifted earlier by
    incubation_period, leaves every row NaN. traj.empty is False (500 rows),
    so the old guard let execution reach groupby.quantile on an empty frame,
    which raised a bare KeyError naming a quantile -- useless to the user.
    """
    traj = simulate_shedding(
        gamma_fit,
        n_individuals=100,
        times=np.arange(0.0, 5.0),
        incubation_period=5.0,
        seed=42,
    )
    assert not traj.empty
    assert traj["log10_value"].isna().all()
    with pytest.raises(ValueError, match="NaN"):
        plot_simulated_shedding(traj, source=gamma_fit)


def test_plot_requires_a_source_to_overlay_observed_points(exponential_fit):
    """observed with source=None must raise rather than plot everything.

    Without a source there is no analyte to filter observed measurements to,
    so silently plotting the whole dataset would be exactly the bug this
    guards against.
    """
    traj = simulate_shedding(
        exponential_fit, n_individuals=10, times=[1.0, 2.0], seed=0
    )
    observed = {
        "participants": [
            {"measurements": [{"analyte": "stool", "time": 1.0, "value": 100.0}]}
        ]
    }
    with pytest.raises(ValueError, match="source"):
        plot_simulated_shedding(traj, source=None, observed=observed)


def _stub_fit_for_plotting(analyte, dataset_id):
    from shedding_hub.shedding_fit import SheddingFit

    return SheddingFit(
        model="exponential",
        method="mle",
        population_mean=np.array([np.log(0.6), np.log(18.0)]),
        population_cov=np.diag([0.04, 0.04]),
        sigma=0.3,
        subject_params=None,
        censoring_limit=2.0,
        dataset_id=dataset_id,
        analyte=analyte,
        biomarker="SARS-CoV-2",
        specimen="stool",
        reference_event="symptom onset",
        unit="gc/mL",
        gene_target=None,
        dose=None,
        vaccine_type=None,
        n_subjects=10,
        n_measurements=100,
        n_censored=10,
        n_excluded_subjects=0,
        n_dropped_measurements=0,
        converged=True,
        log_likelihood=-1.0,
        aic=2.0,
    )


def test_plot_filters_observed_points_to_the_ensembles_own_analytes():
    """An ensemble source must overlay only its components' analytes.

    Before the fix, getattr(source, "analyte", None) was None for any
    SheddingEnsemble (which has no analyte of its own), so the filter was
    skipped entirely and every measurement in `observed` was plotted --
    mixing unrelated analytes on one axis. A SheddingEnsemble's analyte set is
    the union of its component fits' analytes.
    """
    from matplotlib.collections import PathCollection

    from shedding_hub.shedding_ensemble import make_ensemble

    ensemble = make_ensemble(
        [
            _stub_fit_for_plotting("stool_a", "study_a"),
            _stub_fit_for_plotting("stool_b", "study_b"),
        ],
        weights="equal",
    )
    traj = simulate_shedding(
        ensemble, n_individuals=20, times=np.arange(0.0, 20.0), seed=0
    )
    observed = {
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool_a", "time": 1.0, "value": 1e3},
                    {"analyte": "stool_b", "time": 2.0, "value": 1e3},
                    # A third analyte not in this ensemble -- must be excluded.
                    {"analyte": "sputum_other", "time": 3.0, "value": 1e3},
                ]
            }
        ]
    }

    fig = plot_simulated_shedding(traj, source=ensemble, observed=observed)
    ax = fig.axes[0]
    scatter = next(c for c in ax.collections if isinstance(c, PathCollection))
    assert scatter.get_offsets().shape[0] == 2


def test_public_exports_are_available():
    import shedding_hub as sh

    for name in [
        "fit_shedding_model",
        "fit_shedding_models",
        "load_shedding_catalog",
        "make_ensemble",
        "simulate_shedding",
        "plot_simulated_shedding",
        "SheddingFit",
        "SheddingEnsemble",
        "SheddingCatalog",
    ]:
        assert hasattr(sh, name), name


# ---------------------------------------------------------------------------
# dispersion
# ---------------------------------------------------------------------------


def _spread(fit, dispersion, day=0.0, n=4000):
    traj = simulate_shedding(
        fit,
        n_individuals=n,
        times=np.array([day]),
        dispersion=dispersion,
        seed=11,
    )
    return traj["log10_value"].std()


def test_dispersion_defaults_to_leaving_the_fitted_covariance_alone(exponential_fit):
    times = np.arange(0.0, 6.0)
    plain = simulate_shedding(exponential_fit, n_individuals=50, times=times, seed=3)
    explicit = simulate_shedding(
        exponential_fit, n_individuals=50, times=times, dispersion=1.0, seed=3
    )
    pd.testing.assert_frame_equal(plain, explicit)


def test_dispersion_below_one_narrows_the_cohort(exponential_fit):
    """The knob users reach for when a few agents dominate total shed load."""
    assert _spread(exponential_fit, 0.5) < _spread(exponential_fit, 1.0)


def test_dispersion_scales_the_standard_deviation_it_multiplies(exponential_fit):
    """Sigma is scaled by dispersion**2, so the spread scales by dispersion."""
    full = _spread(exponential_fit, 1.0)
    half = _spread(exponential_fit, 0.5)
    assert half == pytest.approx(0.5 * full, rel=0.1)


def test_dispersion_of_zero_gives_every_agent_the_median_individual(exponential_fit):
    from shedding_hub.shedding_models import log10_concentration

    times = np.array([0.0, 5.0])
    traj = simulate_shedding(
        exponential_fit, n_individuals=25, times=times, dispersion=0.0, seed=3
    )
    # Every agent identical *at each time*; the two times of course differ.
    assert (traj.groupby("time")["log10_value"].std() < 1e-9).all()
    expected = log10_concentration(
        exponential_fit.model, exponential_fit.median_params[None, :], times
    )[0]
    np.testing.assert_allclose(
        traj.groupby("time")["log10_value"].first().to_numpy(), expected
    )


def test_dispersion_preserves_the_median(exponential_fit):
    """Shrinking spread must not move the centre of the cohort."""
    wide = simulate_shedding(
        exponential_fit, n_individuals=4000, times=np.array([0.0]), seed=11
    )["log10_value"].median()
    narrow = simulate_shedding(
        exponential_fit,
        n_individuals=4000,
        times=np.array([0.0]),
        dispersion=0.4,
        seed=11,
    )["log10_value"].median()
    assert narrow == pytest.approx(wide, abs=0.15)


def test_dispersion_rejects_a_negative_value(exponential_fit):
    with pytest.raises(ValueError, match="dispersion"):
        simulate_shedding(
            exponential_fit,
            n_individuals=5,
            times=np.array([0.0]),
            dispersion=-0.5,
        )


def _fit_with_reference_event(make_synthetic_dataset, event):
    import numpy as np

    from shedding_hub.shedding_fit import fit_shedding_model

    dataset = make_synthetic_dataset(
        "exponential",
        np.array([np.log(0.6), np.log(18.0)]),
        np.diag([0.04, 0.04]),
        n_subjects=12,
        seed=5,
    )
    for analyte in dataset["analytes"].values():
        analyte["reference_event"] = event
    return fit_shedding_model(dataset, analyte="stool", model="exponential")


def test_symptom_onset_shifts_to_infection(make_synthetic_dataset):
    import numpy as np

    from shedding_hub import simulate_shedding

    fit = _fit_with_reference_event(make_synthetic_dataset, "symptom onset")
    traj = simulate_shedding(
        fit, n_individuals=5, times=np.arange(0, 6), incubation_period=5.0, seed=1
    )
    assert traj.attrs["time_origin"] == "infection"
    assert traj.attrs["reference_event_class"] == "landmark"


def test_administrative_event_warns_and_does_not_claim_infection(
    make_synthetic_dataset,
):
    import numpy as np
    import pytest

    from shedding_hub import simulate_shedding

    fit = _fit_with_reference_event(make_synthetic_dataset, "enrollment")
    with pytest.warns(UserWarning, match="administrative"):
        traj = simulate_shedding(
            fit, n_individuals=5, times=np.arange(0, 6), incubation_period=5.0, seed=1
        )
    assert traj.attrs["time_origin"] == "enrollment_shifted"
    assert traj.attrs["reference_event_class"] == "administrative"


def test_exposure_event_warns_because_there_is_nothing_to_bridge(
    make_synthetic_dataset,
):
    import numpy as np
    import pytest

    from shedding_hub import simulate_shedding

    fit = _fit_with_reference_event(make_synthetic_dataset, "inoculation")
    with pytest.warns(UserWarning, match="already the exposure"):
        traj = simulate_shedding(
            fit, n_individuals=5, times=np.arange(0, 6), incubation_period=5.0, seed=1
        )
    assert traj.attrs["time_origin"] == "inoculation_shifted"


def test_no_incubation_period_leaves_the_origin_alone(make_synthetic_dataset):
    import numpy as np

    from shedding_hub import simulate_shedding

    fit = _fit_with_reference_event(make_synthetic_dataset, "enrollment")
    traj = simulate_shedding(fit, n_individuals=5, times=np.arange(0, 6), seed=1)
    assert traj.attrs["time_origin"] == "enrollment"
    assert traj.attrs["incubation_applied"] is False


def _wide_traj(exponential_fit, low_tail=False):
    """A cohort whose band descends far below any plottable concentration."""
    return simulate_shedding(
        exponential_fit,
        n_individuals=400,
        times=np.arange(1.0, 25.0),
        dispersion=3.0 if low_tail else 1.0,
        seed=5,
    )


def test_plot_floors_the_axis_rather_than_following_the_band_down(exponential_fit):
    """
    A range band descends past concentrations nobody needs plotted.

    Left to autoscale it squashes every real value into a ribbon at the top of
    the panel: measured on the shipped catalog, a gamma_shifted cohort dragged
    the axis to -137 log10.
    """
    from shedding_hub.shedding_simulate import SIMULATION_YLIM_FLOOR

    traj = _wide_traj(exponential_fit, low_tail=True)
    fig = plot_simulated_shedding(traj, source=exponential_fit)
    bottom, _ = fig.axes[0].get_ylim()
    assert bottom >= SIMULATION_YLIM_FLOOR
    # and the floor is genuinely binding on this cohort
    assert traj["log10_value"].min() < SIMULATION_YLIM_FLOOR


def test_the_floor_never_hides_a_real_observation(exponential_fit):
    """It bounds the simulated band, never the data the study actually saw."""
    from shedding_hub.shedding_simulate import SIMULATION_YLIM_FLOOR

    traj = _wide_traj(exponential_fit, low_tail=True)
    observed = {
        "participants": [
            {
                "measurements": [
                    {"analyte": exponential_fit.analyte, "time": 3.0, "value": 10.0**-6}
                ]
            }
        ]
    }
    fig = plot_simulated_shedding(traj, source=exponential_fit, observed=observed)
    bottom, _ = fig.axes[0].get_ylim()
    assert bottom <= -6.0, "an observation below the floor must keep its place"
    assert bottom < SIMULATION_YLIM_FLOOR


def test_plot_leaves_an_ordinary_band_alone(exponential_fit):
    """The floor must not drag a well-behaved axis down to meet it."""
    from shedding_hub.shedding_simulate import SIMULATION_YLIM_FLOOR

    traj = _wide_traj(exponential_fit)
    fig = plot_simulated_shedding(traj, source=exponential_fit)
    bottom, _ = fig.axes[0].get_ylim()
    assert bottom > SIMULATION_YLIM_FLOOR


def test_plot_draws_the_inner_interval_as_dashed_lines(exponential_fit):
    """Two dashed edges, one legend entry between them."""
    traj = _wide_traj(exponential_fit)
    fig = plot_simulated_shedding(
        traj, source=exponential_fit, band_inner_quantiles=(0.025, 0.975)
    )
    ax = fig.axes[0]
    dashed = [line for line in ax.get_lines() if line.get_linestyle() == "--"]
    assert len(dashed) == 2
    labelled = [l for l in dashed if not l.get_label().startswith("_")]
    assert len(labelled) == 1
    assert "95%" in labelled[0].get_label()


def test_inner_interval_can_be_switched_off(exponential_fit):
    traj = _wide_traj(exponential_fit)
    fig = plot_simulated_shedding(
        traj, source=exponential_fit, band_inner_quantiles=None
    )
    dashed = [l for l in fig.axes[0].get_lines() if l.get_linestyle() == "--"]
    assert dashed == []


def test_a_full_range_band_is_labelled_by_draw_count(exponential_fit):
    """
    The extremes a range reaches depend on how many individuals were drawn.

    Calling it "100% of individuals" would imply the population reaches there,
    so the label names the cohort size instead -- the same wording the review
    pages use.
    """
    traj = _wide_traj(exponential_fit)
    fig = plot_simulated_shedding(
        traj, source=exponential_fit, band_quantiles=(0.0, 1.0)
    )
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert any("range" in l and "400" in l for l in labels), labels
    assert not any("100% of individuals" in l for l in labels), labels


def test_the_bands_extremes_do_not_drive_the_axis(exponential_fit):
    """
    Shading the range must not hand the axis to one absurd agent.

    Over a few hundred draws the range reaches concentrations no biology
    supports -- 10**76 gc/mL was measured on the shipped catalog -- so the axis
    follows the inner interval and the band is allowed to clip.
    """
    traj = _wide_traj(exponential_fit, low_tail=True)
    fig = plot_simulated_shedding(
        traj,
        source=exponential_fit,
        band_quantiles=(0.0, 1.0),
        band_inner_quantiles=(0.025, 0.975),
    )
    bottom, top = fig.axes[0].get_ylim()
    inner_high = traj["log10_value"].quantile(0.975)
    band_high = traj["log10_value"].max()
    assert top < band_high, "the band's extreme must not set the top"
    assert top >= inner_high, "the inner interval must remain visible"


def test_an_observation_above_the_interval_still_fits(exponential_fit):
    """The axis bounds the band, never the data — at the top as at the bottom."""
    traj = _wide_traj(exponential_fit)
    high = float(traj["log10_value"].quantile(0.975)) + 5.0
    observed = {
        "participants": [
            {
                "measurements": [
                    {
                        "analyte": exponential_fit.analyte,
                        "time": 3.0,
                        "value": 10.0**high,
                    }
                ]
            }
        ]
    }
    fig = plot_simulated_shedding(traj, source=exponential_fit, observed=observed)
    _, top = fig.axes[0].get_ylim()
    assert top >= high
