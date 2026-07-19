import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from shedding_hub.shedding_catalog import fit_shedding_models
from shedding_hub.shedding_ensemble import SheddingEnsemble, make_ensemble
from shedding_hub.shedding_fit import SheddingFit


def _stub_fit(dataset_id, mean, cov, n_subjects, sigma=0.3):
    """A directly-constructed SheddingFit for hand-computable ensemble tests.

    Only ``population_mean``/``population_cov``/``sigma``/``n_subjects`` vary
    across the tests that use this; every other field is a plausible
    placeholder so the dataclass can be built at all.
    """
    return SheddingFit(
        model="exponential",
        method="mle",
        population_mean=np.asarray(mean, float),
        population_cov=np.asarray(cov, float),
        sigma=sigma,
        subject_params=pd.DataFrame(),
        censoring_limit=2.0,
        dataset_id=dataset_id,
        analyte="stool",
        biomarker="SARS-CoV-2",
        specimen="stool",
        reference_event="symptom onset",
        unit="gc/mL",
        gene_target=None,
        dose=None,
        vaccine_type=None,
        n_subjects=n_subjects,
        n_measurements=100,
        n_censored=10,
        n_excluded_subjects=0,
        n_dropped_measurements=0,
        converged=True,
        log_likelihood=-1.0,
        aic=2.0,
    )


@pytest.fixture
def catalog(make_synthetic_dataset):
    cov = np.diag([0.04, 0.04])

    def study(dataset_id, mu, seed):
        return make_synthetic_dataset(
            "exponential",
            np.asarray(mu),
            cov,
            n_subjects=20,
            seed=seed,
            dataset_id=dataset_id,
        )

    return fit_shedding_models(
        [
            study("study_a", [np.log(0.6), np.log(18.0)], 1),
            study("study_b", [np.log(0.4), np.log(20.0)], 2),
            study("study_c", [np.log(0.8), np.log(16.0)], 3),
        ],
        models=("exponential",),
    )


def test_ensemble_over_all_matching_studies(catalog):
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2", specimen="stool")
    assert len(ensemble.fits) == 3
    assert len(ensemble.components) == 3
    assert set(ensemble.components["dataset_id"]) == {
        "study_a",
        "study_b",
        "study_c",
    }


def test_ensemble_restricted_to_named_studies(catalog):
    ensemble = catalog.ensemble(
        biomarker="SARS-CoV-2", dataset_ids=["study_a", "study_c"]
    )
    assert set(ensemble.components["dataset_id"]) == {"study_a", "study_c"}


def test_unmatched_dataset_id_raises(catalog):
    with pytest.raises(ValueError, match="study_z"):
        catalog.ensemble(biomarker="SARS-CoV-2", dataset_ids=["study_a", "study_z"])


def test_single_component_ensemble_matches_the_underlying_fit(catalog):
    fit = catalog.select(dataset_id="study_a")
    ensemble = make_ensemble([fit])
    a = fit.sample_params(np.random.default_rng(0), 50)[0]
    b = ensemble.sample_params(np.random.default_rng(0), 50)[0]
    np.testing.assert_allclose(a, b)


def test_mixture_draws_from_every_study_and_records_the_source(catalog):
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2", weights="equal")
    _, sources = ensemble.sample_params(np.random.default_rng(0), 600)
    assert set(sources.tolist()) == {"study_a", "study_b", "study_c"}
    counts = np.array(
        [(sources == name).sum() for name in ["study_a", "study_b", "study_c"]]
    )
    # Equal weights: each study should take roughly a third.
    assert (counts > 120).all()


def test_weights_default_to_subject_counts(catalog):
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2")
    expected = np.array([fit.n_subjects for fit in ensemble.fits], dtype=float)
    np.testing.assert_allclose(ensemble.weights, expected / expected.sum())


def test_negative_weight_entries_raise(catalog):
    fits = list(catalog.fits)
    with pytest.raises(ValueError, match="non-negative"):
        make_ensemble(fits, weights=[1.0, -2.0, 3.0])


def test_weights_accepts_a_numpy_array(catalog):
    """A numpy array must reach the length/negativity validation, not crash.

    ``if weights == "equal"`` on a numpy array is an elementwise comparison
    that raises "The truth value of an array ... is ambiguous" before the
    validation below it ever runs -- even though the docstring documents "an
    explicit array" as a valid input for a numpy-centric package.
    """
    fits = list(catalog.fits)
    array_weights = np.array([1.0, 2.0, 3.0])
    ensemble = make_ensemble(fits, weights=array_weights)
    np.testing.assert_allclose(ensemble.weights, array_weights / array_weights.sum())


def test_weights_array_still_validates_length_and_negativity(catalog):
    fits = list(catalog.fits)
    with pytest.raises(ValueError, match="length 3"):
        make_ensemble(fits, weights=np.array([1.0, 2.0]))
    with pytest.raises(ValueError, match="non-negative"):
        make_ensemble(fits, weights=np.array([1.0, -2.0, 3.0]))


def test_unrecognized_weights_string_raises(catalog):
    fits = list(catalog.fits)
    with pytest.raises(ValueError, match="n_subjects.*equal"):
        make_ensemble(fits, weights="bogus")


def test_moment_covariance_is_within_plus_between():
    """Hand-computable two-study example."""
    a = _stub_fit("a", [0.0, 0.0], np.eye(2), 10)
    b = _stub_fit("b", [2.0, 0.0], np.eye(2), 10)
    ensemble = make_ensemble([a, b], weights="equal", method="moment")

    np.testing.assert_allclose(ensemble.population_mean, [1.0, 0.0])
    # within = I; between = weighted cov of means = [[1, 0], [0, 0]]
    expected = np.eye(2) + np.array([[1.0, 0.0], [0.0, 0.0]])
    np.testing.assert_allclose(ensemble.population_cov, expected)
    np.testing.assert_allclose(ensemble.median_params, np.exp([1.0, 0.0]))


def test_moment_sample_params_rejects_non_positive_semi_definite_covariance():
    """The moment path computes its own covariance and must guard it too.

    Two components with means far enough apart along an axis where each
    component's own covariance is (numerically) near zero can leave the
    moment-matched covariance non-PSD-clean only through floating point noise,
    so instead we directly force a non-PSD case by constructing components
    whose within/between combination is unambiguously indefinite.
    """
    # Both components share a mean (so "between" contributes nothing) and an
    # indefinite own-covariance (eigenvalues 3 and -1), so within + between is
    # exactly that indefinite matrix -- unambiguously not PSD, not numerical
    # noise near zero.
    non_psd_cov = np.array([[1.0, 2.0], [2.0, 1.0]])
    a = _stub_fit("a", [0.0, 0.0], non_psd_cov, 10)
    b = _stub_fit("b", [0.0, 0.0], non_psd_cov, 10)
    ensemble = make_ensemble([a, b], weights="equal", method="moment")

    with pytest.raises(ValueError, match="positive semi-definite"):
        ensemble.sample_params(np.random.default_rng(0), 5)


def test_mixture_sample_params_rejects_a_components_non_positive_semi_definite_covariance():
    """The mixture path must validate each component's own covariance too.

    Previously only the single-component shortcut (which delegates to
    SheddingFit.sample_params) validated; a genuine multi-component mixture
    called rng.multivariate_normal directly on each component's covariance
    inside the loop, skipping _require_positive_semidefinite entirely. This
    also means the moment path's advice to "consider method='mixture' instead,
    which uses ... (already validated) covariances" was false for exactly the
    case it targets -- this test pins that the mixture path now validates, so
    the advice is true.
    """
    non_psd_cov = np.array([[1.0, 2.0], [2.0, 1.0]])
    a = _stub_fit("a", [0.0, 0.0], non_psd_cov, 10)
    b = _stub_fit("b", [5.0, 0.0], np.eye(2), 10)
    ensemble = make_ensemble([a, b], weights="equal", method="mixture")

    with pytest.raises(ValueError, match="positive semi-definite"):
        # n large enough that both components are certain to be drawn from
        # with a fixed seed, so the bad component's covariance is actually used.
        ensemble.sample_params(np.random.default_rng(0), 200)


def test_mixture_sigma_is_root_mean_square_not_arithmetic_mean():
    """Variances add across a mixture; standard deviations do not.

    With equal weights and component sigmas 0.2 and 0.8, the correct combined
    sigma is sqrt(mean(sigma**2)) = sqrt(0.5 * (0.04 + 0.64)) ~= 0.583, which
    is strictly greater than the (wrong) arithmetic mean of 0.5.
    """
    a = _stub_fit("a", [0.0, 0.0], np.eye(2), 10, sigma=0.2)
    b = _stub_fit("b", [0.0, 0.0], np.eye(2), 10, sigma=0.8)
    ensemble = make_ensemble([a, b], weights="equal")

    expected = np.sqrt(np.mean(np.array([0.2, 0.8]) ** 2))
    np.testing.assert_allclose(ensemble.sigma, expected)
    assert ensemble.sigma > np.mean([0.2, 0.8])


def test_mixed_units_raise(catalog):
    fits = list(catalog.fits)
    fits[1].unit = "gc/dry gram"
    with pytest.raises(ValueError, match="unit"):
        make_ensemble(fits)


def test_mixed_reference_events_raise(catalog):
    fits = list(catalog.fits)
    fits[1].reference_event = "enrollment"
    with pytest.raises(ValueError, match="reference_event"):
        make_ensemble(fits)


def test_two_analytes_from_one_study_raise(catalog):
    fits = list(catalog.fits)
    fits[1].dataset_id = "study_a"
    fits[1].analyte = "stool_orf1a"
    with pytest.raises(ValueError, match="more than one"):
        make_ensemble(fits)


def test_ensemble_can_be_simulated(catalog):
    from shedding_hub.shedding_simulate import simulate_shedding

    ensemble = catalog.ensemble(biomarker="SARS-CoV-2")
    traj = simulate_shedding(ensemble, n_individuals=30, times=[1.0, 5.0], seed=0)
    assert len(traj) == 60
    assert traj["source_dataset_id"].nunique() > 1


def test_mixture_has_no_median_params(catalog):
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2", method="mixture")
    with pytest.raises(ValueError, match="mixture"):
        _ = ensemble.median_params


def test_ensemble_to_dict_from_dict_round_trip(catalog):
    """SheddingEnsemble.to_dict()/from_dict() completes the persistence story:

    a fit can be saved and reloaded (SheddingFit.to_dict/from_dict), and so
    can an ensemble built from several such fits -- its component fits,
    resolved weights, and combination method all round-trip.
    """
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2", weights="equal")

    payload = ensemble.to_dict()
    restored = SheddingEnsemble.from_dict(payload)

    assert len(restored.fits) == len(ensemble.fits)
    assert [fit.dataset_id for fit in restored.fits] == [
        fit.dataset_id for fit in ensemble.fits
    ]
    assert all(fit.subject_params is None for fit in restored.fits)
    np.testing.assert_allclose(restored.weights, ensemble.weights)
    assert restored.method == ensemble.method
    np.testing.assert_allclose(restored.population_mean, ensemble.population_mean)


def test_deserialized_ensemble_can_still_simulate(catalog):
    """The property that lets an ABM fit once and simulate across many runs."""
    from shedding_hub.shedding_simulate import simulate_shedding

    ensemble = catalog.ensemble(biomarker="SARS-CoV-2")
    restored = SheddingEnsemble.from_dict(ensemble.to_dict())
    traj = simulate_shedding(restored, n_individuals=30, times=[1.0, 5.0], seed=0)
    assert len(traj) == 60
    assert traj["source_dataset_id"].nunique() > 1
