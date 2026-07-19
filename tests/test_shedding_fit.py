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
    understating b0, which -- because peak_day = b0 / a0 -- makes the fitted
    peak-shedding day somewhat early. This test pins that direction and
    mechanism (not just its existence) so that a future change to the
    aggregation that makes the bias worse, or removes the effect entirely
    without a corresponding fix elsewhere, is caught rather than silently
    accepted.

    The margin below was 0.2, calibrated when the bias measured roughly 0.5
    log units. Most of that turned out to be an initialization artifact rather
    than two-stage bias: the gamma fit was seeded from its own collinear
    [1, ln(t), -t] design, whose negative coefficients were clipped onto the
    parameter floor that the optimizer cannot escape. With the fit seeded from
    the well-conditioned decay design instead, the remaining genuine bias
    measures 0.15 log units on average over six seeds (range 0.02 to 0.27), so
    0.2 no longer holds at every seed and the margin is 0.05. Both quantities
    are seed-noisy at 60 subjects; a failure here means checking the
    distribution across seeds before concluding anything has regressed.
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
    assert sparse_error > 0.05
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
    assert list(fit.subject_params.columns) == [
        "subject_id",
        "a0",
        "c0",
        "degenerate",
    ]


def _flat_dataset(n_flat, n_decaying=0):
    """Subjects whose concentration never changes, optionally plus real decays.

    A perfectly flat trajectory has no decay to estimate, so its a0 is driven
    to the positivity floor -- the collapse this fix detects, reproduced
    deliberately rather than waited for.
    """
    participants = [
        {
            "measurements": [
                {"analyte": "stool", "time": t, "value": 1e5} for t in range(1, 11)
            ]
        }
        for _ in range(n_flat)
    ]
    participants += [
        {
            "measurements": [
                {"analyte": "stool", "time": t, "value": 1e7 * np.exp(-rate * t)}
                for t in range(1, 11)
            ]
        }
        for rate in np.linspace(0.4, 0.7, n_decaying)
    ]
    return {
        "dataset_id": "degenerate_study",
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
        "participants": participants,
    }


def test_all_subjects_collapsing_raises_degenerate_fit():
    with pytest.raises(SheddingDataError) as excinfo:
        fit_shedding_model(_flat_dataset(3), analyte="stool", model="exponential")
    assert excinfo.value.reason == "degenerate_fit"
    assert "not identifiable" in str(excinfo.value)


def test_degenerate_subject_is_flagged_but_excluded_from_the_population():
    """The collapsed subject stays inspectable; it just stops voting."""
    dataset = _flat_dataset(n_flat=1, n_decaying=3)
    with pytest.warns(UserWarning, match="collapsed onto the bounds"):
        fit = fit_shedding_model(dataset, analyte="stool", model="exponential")

    assert fit.n_degenerate_subjects == 1
    # Retained in subject_params, flagged, and still the 4 subjects fitted.
    assert len(fit.subject_params) == 4
    assert fit.n_subjects == 4
    assert fit.subject_params["degenerate"].sum() == 1

    # The population summary reflects only the three real decays, whose rates
    # span 0.4 to 0.7 -- the collapsed subject's a0 of ~1e-6 would have dragged
    # the mean of log(a0) down by several units had it been included.
    assert 0.4 <= fit.median_params[0] <= 0.7
    assert fit.half_life_days == pytest.approx(np.log(2) / fit.median_params[0])


def test_population_covariance_ignores_degenerate_subjects():
    """Sigma must be estimated from the retained subjects alone."""
    dataset = _flat_dataset(n_flat=1, n_decaying=3)
    with pytest.warns(UserWarning, match="collapsed onto the bounds"):
        fit = fit_shedding_model(dataset, analyte="stool", model="exponential")

    theta = np.log(
        fit.subject_params.loc[
            ~fit.subject_params["degenerate"], list(fit.param_names)
        ].to_numpy(dtype=float)
    )
    np.testing.assert_allclose(fit.population_mean, theta.mean(axis=0))
    np.testing.assert_allclose(fit.population_cov, np.cov(theta, rowvar=False, ddof=1))


def test_single_subject_fit_still_allowed_when_not_degenerate():
    """A one-subject fit has no estimable covariance, but is not degenerate.

    The 'fewer than two subjects' guard must key off collapse, not off the
    subject count -- test_shedding_tutorial_agreement.py fits exactly one
    subject and must keep working.
    """
    dataset = _flat_dataset(n_flat=0, n_decaying=1)
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert fit.n_subjects == 1
    assert fit.n_degenerate_subjects == 0
    np.testing.assert_allclose(fit.population_cov, np.zeros((2, 2)))


def test_single_degenerate_subject_raises():
    with pytest.raises(SheddingDataError) as excinfo:
        fit_shedding_model(
            _flat_dataset(n_flat=1), analyte="stool", model="exponential"
        )
    assert excinfo.value.reason == "degenerate_fit"


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


# ---------------------------------------------------------------------------
# Real-data regression guards.
#
# Everything above this line is synthetic, and every one of those tests passed
# while the fitter was silently broken. Per-subject fits were initialized at the
# parameter floor _MIN_PARAM and could never escape it -- because the optimizer
# works on theta = log(param), the gradient carries a factor of param and
# vanishes as the parameter approaches zero -- and those collapsed subjects were
# then averaged into the population mean. On woelfel2020virological stool under
# the gamma model, four of nine subjects pinned at a0 = b0 = 1e-6, which
# produced a reported SARS-CoV-2 half-life of 278 days and a peak of 0.31 log10
# (about 2 gc/mL) for data whose concentrations run 10^2 to 10^7. Repository-wide
# it produced 37 fits with half-lives beyond a year.
#
# Synthetic fixtures never reproduced it because they are generated from
# well-behaved parameters over dense, well-conditioned sampling grids. Only real
# data has the sparse, post-peak sampling windows that make the gamma design
# matrix collinear. These tests are therefore the guard that the synthetic suite
# structurally cannot provide, and they assert physical plausibility rather than
# exact numbers.
# ---------------------------------------------------------------------------

import pathlib

import shedding_hub as sh
from shedding_hub.shedding_fit import _MIN_PARAM

DATA = pathlib.Path(__file__).parent.parent / "data"


@pytest.fixture(scope="module")
def woelfel():
    """The real SARS-CoV-2 study the collapse was first observed on."""
    return sh.load_dataset("woelfel2020virological", local=str(DATA))


@pytest.fixture(scope="module")
def woelfel_fits(woelfel):
    """Every analyte/model combination this bug affected, fitted once."""
    return {
        (analyte, model): fit_shedding_model(woelfel, analyte=analyte, model=model)
        for analyte in ("stool", "sputum")
        for model in ("exponential", "gamma")
    }


@pytest.mark.parametrize("analyte", ["stool", "sputum"])
@pytest.mark.parametrize("model", ["exponential", "gamma"])
def test_woelfel_fit_is_physically_plausible(woelfel_fits, analyte, model):
    """
    Real-data guard: a synthetic-only suite missed a 278-day half-life.

    Both analytes are SARS-CoV-2 in gc/mL, so the median individual's peak must
    land inside the assay's actual dynamic range and its decay must be measured
    in days. The collapsed fit this guards against reported peak_log10 = 0.31
    and half_life_days = 278 for stool under the gamma model; some repository
    fits reached ln(2)/1e-6 = 693147 days exactly.

    The bounds are deliberately loose -- these are plausibility limits, not
    pinned values, so that legitimate refits are free to move within them while
    a collapse of the kind described above still fails loudly.
    """
    fit = woelfel_fits[(analyte, model)]
    assert 2.0 < fit.peak_log10 < 10.0, (
        f"{analyte}/{model} peaks at {fit.peak_log10:.3g} log10 gc/mL, outside "
        "the plausible range for a SARS-CoV-2 concentration assay"
    )
    assert 0.1 < fit.half_life_days < 60.0, (
        f"{analyte}/{model} has a half-life of {fit.half_life_days:.4g} days; "
        "SARS-CoV-2 shedding does not decay that slowly"
    )


@pytest.mark.parametrize("analyte", ["stool", "sputum"])
@pytest.mark.parametrize("model", ["exponential", "gamma"])
def test_woelfel_no_subject_sits_at_the_parameter_floor(woelfel_fits, analyte, model):
    """
    Real-data guard: no per-subject fit may collapse to _MIN_PARAM.

    This is the mechanism behind the 278-day half-life rather than its symptom.
    A factor of ten above the floor is the margin: the floor is an absorbing
    state, so a parameter that is genuinely being estimated has no reason to
    linger anywhere near it, and anything that does was placed there by
    initialization rather than found by the optimizer.
    """
    fit = woelfel_fits[(analyte, model)]
    params = fit.subject_params.drop(columns=["subject_id", "degenerate"])
    smallest = params.min().min()
    assert smallest > 10 * _MIN_PARAM, (
        f"{analyte}/{model} has a subject parameter at {smallest:.3g}, within a "
        f"factor of ten of the floor {_MIN_PARAM:.0e}:\n{fit.subject_params}"
    )
