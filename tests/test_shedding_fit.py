import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from shedding_hub.shedding_fit import (
    SheddingDataError,
    prepare_observations,
)


@pytest.fixture
def simple_dataset():
    """Two subjects, one analyte, with a censored point and a qualitative one."""
    return {
        "dataset_id": "test_study",
        "analytes": {
            "stool": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": 1e6},
                    {"analyte": "stool", "time": 2, "value": 1e5},
                    {"analyte": "stool", "time": 3, "value": "negative"},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": 1e7},
                    {"analyte": "stool", "time": 2, "value": 1e6},
                    {"analyte": "stool", "time": 3, "value": "positive"},
                    {"analyte": "stool", "time": "unknown", "value": 1e4},
                ]
            },
        ],
    }


def test_extracts_values_on_log10_scale(simple_dataset):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.n_subjects == 2
    np.testing.assert_allclose(obs.values[0], 6.0)
    np.testing.assert_allclose(obs.values[1], 5.0)


def test_negative_becomes_censored(simple_dataset):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censored.sum() == 1
    assert np.isnan(obs.values[obs.censored][0])


def test_censoring_limit_uses_declared_loq_when_below_smallest_positive(
    simple_dataset,
):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    # LOQ 100 -> log10 = 2, which is below the smallest observed positive (5.0)
    assert obs.censoring_limit == pytest.approx(2.0)


def test_censoring_limit_falls_back_below_smallest_positive(simple_dataset):
    # Declare a limit above every observed value; the fallback must kick in.
    simple_dataset["analytes"]["stool"]["limit_of_quantification"] = 1e8
    with pytest.warns(UserWarning, match="censoring limit"):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censoring_limit == pytest.approx(5.0 - 0.01)


def test_censoring_limit_uses_lod_when_loq_is_unusable(simple_dataset):
    # LOQ is "unknown" (unusable); LOD is a valid number below the smallest
    # observed positive (5.0), so the declared LOD should resolve cleanly,
    # with no fallback-below-smallest-positive warning.
    simple_dataset["analytes"]["stool"]["limit_of_quantification"] = "unknown"
    simple_dataset["analytes"]["stool"]["limit_of_detection"] = 10
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censoring_limit == pytest.approx(1.0)


def test_qualitative_and_unknown_time_are_dropped_with_warning(simple_dataset):
    with pytest.warns(UserWarning):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    # 7 measurements total, minus 1 qualitative "positive", minus 1 unknown time
    assert obs.times.size == 5
    assert obs.n_dropped_measurements == 2


def test_ct_analyte_is_rejected():
    dataset = {
        "dataset_id": "ct_study",
        "analytes": {
            "swab": {
                "specimen": "saliva",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "cycle threshold",
                "limit_of_quantification": "unknown",
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "swab", "time": 1, "value": 20.0},
                    {"analyte": "swab", "time": 2, "value": 25.0},
                ]
            }
        ],
    }
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(dataset, "swab", "exponential")
    assert excinfo.value.reason == "ct_units"


def test_non_pathogen_biomarker_is_rejected():
    # Plenty of clean, numeric, positive-time data for both subjects, so a
    # rejection here can only be the biomarker check, not a data shortage.
    dataset = {
        "dataset_id": "crassphage_study",
        "analytes": {
            "stool_crAssphage": {
                "specimen": "stool",
                "biomarker": "crAssphage",
                "reference_event": "symptom onset",
                "unit": "gc/dry gram",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool_crAssphage", "time": 1, "value": 1e6},
                    {"analyte": "stool_crAssphage", "time": 2, "value": 1e5},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool_crAssphage", "time": 1, "value": 1e7},
                    {"analyte": "stool_crAssphage", "time": 2, "value": 1e6},
                ]
            },
        ],
    }
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(dataset, "stool_crAssphage", "exponential")
    assert excinfo.value.reason == "non_pathogen_biomarker"


def test_vaccine_strain_biomarker_is_not_rejected():
    # Live-attenuated vaccine shedding has a real trajectory (reference event
    # "vaccination"), unlike the fecal-strength/normalization indicators in
    # NON_PATHOGEN_BIOMARKERS. Guards against that set later being "tidied"
    # into a broader substring/category match that would catch this too.
    dataset = {
        "dataset_id": "vaccine_study",
        "analytes": {
            "stool_rotavirus": {
                "specimen": "stool",
                "biomarker": "rotavirus vaccine",
                "reference_event": "vaccination",
                "unit": "gc/mL",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool_rotavirus", "time": 1, "value": 1e6},
                    {"analyte": "stool_rotavirus", "time": 2, "value": 1e5},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool_rotavirus", "time": 1, "value": 1e7},
                    {"analyte": "stool_rotavirus", "time": 2, "value": 1e6},
                ]
            },
        ],
    }
    obs = prepare_observations(dataset, "stool_rotavirus", "exponential")
    assert obs.n_subjects == 2


def test_no_positive_measurements_raises():
    # Both subjects have enough censored measurements to clear
    # min_observations on their own, so they are retained rather than
    # excluded — the failure must come from having no positives at all.
    dataset = {
        "dataset_id": "all_negative_study",
        "analytes": {
            "stool": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": "negative"},
                    {"analyte": "stool", "time": 2, "value": "negative"},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": "negative"},
                    {"analyte": "stool", "time": 2, "value": "negative"},
                ]
            },
        ],
    }
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(dataset, "stool", "exponential")
    assert excinfo.value.reason == "no_positive_measurements"


def test_gamma_drops_non_positive_times(simple_dataset):
    simple_dataset["participants"][0]["measurements"].append(
        {"analyte": "stool", "time": 0, "value": 1e6}
    )
    simple_dataset["participants"][0]["measurements"].append(
        {"analyte": "stool", "time": -2, "value": 1e6}
    )
    with pytest.warns(UserWarning):
        gamma_obs = prepare_observations(simple_dataset, "stool", "gamma")
    exponential_obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert (gamma_obs.times > 0).all()
    # Two effects separate exponential from gamma here: the two appended
    # non-positive times are dropped outright under gamma, and subject 2 (only
    # 2 usable measurements) falls below gamma's 3-observation minimum and is
    # excluded entirely, while it clears exponential's 2-observation minimum.
    assert exponential_obs.times.size == gamma_obs.times.size + 4


def test_subjects_below_min_observations_are_excluded(simple_dataset):
    simple_dataset["participants"].append(
        {"measurements": [{"analyte": "stool", "time": 1, "value": 1e5}]}
    )
    with pytest.warns(UserWarning, match="excluded"):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.n_subjects == 2
    assert obs.n_excluded_subjects == 1


def test_no_usable_subject_raises(simple_dataset):
    for participant in simple_dataset["participants"]:
        participant["measurements"] = participant["measurements"][:1]
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(simple_dataset, "stool", "exponential")
    assert excinfo.value.reason == "too_few_subjects"


def test_unknown_analyte_raises(simple_dataset):
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(simple_dataset, "sputum", "exponential")
    assert excinfo.value.reason == "unknown_analyte"


def test_subject_index_is_contiguous_after_exclusions(simple_dataset):
    simple_dataset["participants"].insert(
        0, {"measurements": [{"analyte": "stool", "time": 1, "value": 1e5}]}
    )
    with pytest.warns(UserWarning):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert set(obs.subject_index.tolist()) == {0, 1}


from shedding_hub.shedding_fit import SheddingFit, fit_shedding_model


def test_recovers_known_exponential_population(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.09, 0.04]), sigma=0.3, n_subjects=60
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    np.testing.assert_allclose(fit.population_mean, mu, atol=0.15)
    assert fit.sigma == pytest.approx(0.3, abs=0.15)
    assert fit.converged


def test_recovers_known_gamma_population(make_synthetic_dataset):
    mu = np.array([np.log(0.5), np.log(1.5), np.log(12.0)])
    # A dense sampling grid is deliberate: the gamma model's b0 (rise-rate/
    # shape) MLE has an O(1/observations-per-subject) finite-sample bias, so
    # the default ~14-point grid would mostly test that bias rather than
    # whether the estimator is correct. 56 points/subject is dense enough for
    # all three parameters to land within atol=0.25; see
    # test_gamma_b0_is_downward_biased_at_sparse_sampling for the sparse-grid
    # behavior, which is pinned separately rather than papered over here.
    dataset = make_synthetic_dataset(
        "gamma",
        mu,
        np.diag([0.04, 0.04, 0.04]),
        sigma=0.3,
        n_subjects=60,
        times=np.linspace(1.0, 14.0, 56),
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="gamma")
    np.testing.assert_allclose(fit.population_mean, mu, atol=0.25)
    assert fit.converged


def test_gamma_b0_is_downward_biased_at_sparse_sampling(make_synthetic_dataset):
    """
    Documents a known limitation rather than a defect awaiting a fix.

    The gamma model's b0 MLE is finite-sample biased downward, roughly
    O(1/observations-per-subject) (confirmed by refitting the same truth at
    14/28/56/112 observations per subject: the bias shrinks by about half
    each time the density doubles, then vanishes). At the ~14-point sampling
    density typical of real shedding studies, this shows up as population_mean
    understating b0 by roughly 0.5 log units, which -- because
    peak_day = b0 / a0 -- makes the fitted peak-shedding day systematically
    early. This test pins that direction and mechanism (not just its
    existence) so that a future change to the aggregation that makes the bias
    worse, or removes the effect entirely without a corresponding fix
    elsewhere, is caught rather than silently accepted.
    """
    mu = np.array([np.log(0.5), np.log(1.5), np.log(12.0)])
    cov = np.diag([0.04, 0.04, 0.04])

    sparse = make_synthetic_dataset("gamma", mu, cov, sigma=0.3, n_subjects=60)
    sparse_fit = fit_shedding_model(sparse, analyte="stool", model="gamma")

    dense = make_synthetic_dataset(
        "gamma", mu, cov, sigma=0.3, n_subjects=60, times=np.linspace(1.0, 14.0, 56)
    )
    dense_fit = fit_shedding_model(dense, analyte="stool", model="gamma")

    true_b0 = mu[1]
    sparse_error = true_b0 - sparse_fit.population_mean[1]
    dense_error = abs(dense_fit.population_mean[1] - true_b0)

    # (a) the sparse fit understates b0 by a clear margin.
    assert sparse_error > 0.2
    # (b) denser sampling recovers b0 closer to the truth than sparse sampling
    # -- the direction and mechanism (more data shrinks this specific bias),
    # not just that some bias exists.
    assert dense_error < sparse_error


def test_censored_fit_beats_dropping_negatives(make_synthetic_dataset):
    """The property the whole design rests on."""
    mu = np.array([np.log(0.6), np.log(14.0)])
    # A high limit censors much of the decay phase.
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), sigma=0.2, n_subjects=60, loq=1e4
    )

    censored_fit = fit_shedding_model(dataset, analyte="stool", model="exponential")

    dropped = {
        **dataset,
        "participants": [
            {"measurements": [m for m in p["measurements"] if m["value"] != "negative"]}
            for p in dataset["participants"]
        ],
    }
    naive_fit = fit_shedding_model(dropped, analyte="stool", model="exponential")

    true_decay = mu[0]
    censored_error = abs(censored_fit.population_mean[0] - true_decay)
    naive_error = abs(naive_fit.population_mean[0] - true_decay)
    assert censored_error < naive_error
    # Dropping negatives should understate the decay rate.
    assert naive_fit.population_mean[0] < censored_fit.population_mean[0]


def test_fit_carries_metadata_and_counts(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    # loq=1e2 (the fixture default) never censors this particular truth curve
    # over t=1..14 -- its noisy values stay in [2.9, 9.6] on the log10 scale,
    # nowhere near log10(1e2)=2. A higher loq is needed so this dataset
    # actually exercises the censored branch it claims to check.
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=10, loq=1e5
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert isinstance(fit, SheddingFit)
    assert fit.dataset_id == "synthetic"
    assert fit.analyte == "stool"
    assert fit.biomarker == "SARS-CoV-2"
    assert fit.specimen == "stool"
    assert fit.reference_event == "symptom onset"
    assert fit.unit == "gc/mL"
    assert fit.method == "mle"
    assert fit.n_subjects == 10
    # 10 subjects x 14 time points, none dropped by the exponential model.
    assert fit.n_measurements == 140
    assert 0 < fit.n_censored < fit.n_measurements
    assert fit.param_names == ("a0", "c0")


def test_median_params_are_exp_of_population_mean(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=10
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    np.testing.assert_allclose(fit.median_params, np.exp(fit.population_mean))


def test_gamma_peak_day_is_b_over_a(make_synthetic_dataset):
    mu = np.array([np.log(0.5), np.log(2.0), np.log(12.0)])
    dataset = make_synthetic_dataset(
        "gamma", mu, np.diag([0.04, 0.04, 0.04]), n_subjects=20
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="gamma")
    assert fit.peak_day == pytest.approx(fit.median_params[1] / fit.median_params[0])


def test_sample_params_shape_and_source(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=20
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    rng = np.random.default_rng(1)
    params, sources = fit.sample_params(rng, 25)
    assert params.shape == (25, 2)
    assert (params > 0).all()
    assert sources.shape == (25,)
    assert set(sources.tolist()) == {"synthetic"}


def test_subject_params_has_one_row_per_subject(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=12
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert len(fit.subject_params) == 12
    assert list(fit.subject_params.columns) == ["subject_id", "a0", "c0"]


def _minimal_fit(population_mean, population_cov, model="exponential"):
    """A directly-constructed SheddingFit, bypassing fit_shedding_model entirely.

    Only ``population_mean``/``population_cov`` matter for the tests that use
    this; every other field is a plausible placeholder so the dataclass can be
    built at all.
    """
    return SheddingFit(
        model=model,
        method="mle",
        population_mean=np.asarray(population_mean, dtype=float),
        population_cov=np.asarray(population_cov, dtype=float),
        sigma=0.3,
        subject_params=None,
        censoring_limit=2.0,
        dataset_id="synthetic",
        analyte="stool",
        biomarker="SARS-CoV-2",
        specimen="stool",
        reference_event="symptom onset",
        unit="gc/mL",
        gene_target=None,
        dose=None,
        vaccine_type=None,
        n_subjects=1,
        n_measurements=1,
        n_censored=0,
        n_excluded_subjects=0,
        n_dropped_measurements=0,
        converged=True,
        log_likelihood=0.0,
        aic=0.0,
    )


def test_sample_params_rejects_non_positive_semi_definite_covariance():
    # Eigenvalues 3 and -1: unambiguously not positive semi-definite, not just
    # numerical noise near zero.
    non_psd_cov = np.array([[1.0, 2.0], [2.0, 1.0]])
    fit = _minimal_fit([np.log(0.6), np.log(18.0)], non_psd_cov)
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="positive semi-definite"):
        fit.sample_params(rng, 5)


def test_sample_params_accepts_zero_covariance_producing_identical_individuals():
    # A single-subject fit legitimately has an all-zero population_cov (no
    # between-subject variance can be estimated from one subject). This must
    # still simulate, giving every individual the same parameters.
    mean = np.array([np.log(0.6), np.log(18.0)])
    fit = _minimal_fit(mean, np.zeros((2, 2)))
    rng = np.random.default_rng(0)
    params, sources = fit.sample_params(rng, 4)
    np.testing.assert_allclose(params, np.tile(np.exp(mean), (4, 1)))
    assert set(sources.tolist()) == {"synthetic"}
