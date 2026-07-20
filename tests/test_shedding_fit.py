import matplotlib

matplotlib.use("Agg")

import math

import numpy as np
import pytest

from shedding_hub.shedding_fit import (
    _MIN_HALF_LIFE_DAYS,
    SheddingDataError,
    _degenerate_subjects,
    _fraction_observing_a_rise,
    prepare_observations,
    require_estimable_population,
)
from shedding_hub.shedding_models import PARAM_NAMES


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


def test_censoring_limit_falls_back_when_no_limit_declared(simple_dataset):
    # With neither a limit of quantification nor detection declared, fall back to
    # just below the smallest observed positive (5.0) so any `negative` still
    # sits below the resolved limit.
    simple_dataset["analytes"]["stool"]["limit_of_quantification"] = "unknown"
    simple_dataset["analytes"]["stool"]["limit_of_detection"] = "unknown"
    with pytest.warns(UserWarning, match="censoring limit"):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censoring_limit == pytest.approx(5.0 - 0.01)


def test_declared_limit_above_all_positives_is_still_used(simple_dataset):
    # A declared limit is always used as-is, even when it sits above every
    # observed positive: positives stay as data, only `negative`s are censored
    # (here at log10 1e8 = 8.0). No data-derived fallback.
    simple_dataset["analytes"]["stool"]["limit_of_quantification"] = 1e8
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censoring_limit == pytest.approx(8.0)


def test_censoring_limit_uses_lod_when_loq_is_unusable(simple_dataset):
    # LOQ is "unknown" (unusable); LOD is a valid number below the smallest
    # observed positive (5.0), so the declared LOD should resolve cleanly,
    # with no fallback-below-smallest-positive warning.
    simple_dataset["analytes"]["stool"]["limit_of_quantification"] = "unknown"
    simple_dataset["analytes"]["stool"]["limit_of_detection"] = 10
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censoring_limit == pytest.approx(1.0)


def test_positive_below_loq_is_kept_as_observed(simple_dataset):
    # A positive reading below the declared LOQ of 100 (log10 2.0): 50 gc/mL is
    # still a real measurement, so it is kept as observed data rather than
    # censored, and the declared limit (2.0) is used as-is for the `negative`s —
    # it is not lowered to accommodate the sub-LOQ positive.
    simple_dataset["participants"][0]["measurements"].append(
        {"analyte": "stool", "time": 4, "value": 50.0}  # log10 ~1.70, below LOQ
    )
    obs = prepare_observations(simple_dataset, "stool", "exponential")

    # Declared LOQ is used as-is, not lowered to a data-derived fallback.
    assert obs.censoring_limit == pytest.approx(2.0)

    # The 50 gc/mL reading is kept as an observed (non-censored) value, and an
    # observed value below the censoring limit is allowed.
    subj0 = obs.subject_index == 0
    observed_subj0 = obs.values[subj0 & ~obs.censored]
    assert np.any(np.isclose(observed_subj0, np.log10(50.0)))
    assert observed_subj0.min() < obs.censoring_limit

    # Only the original `negative` (t=3) is censored for subject 1.
    assert obs.censored[subj0].sum() == 1


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


def _cliff_dataset(n_cliff, n_normal):
    """Subjects that vanish between two daily samples, plus ordinary decays.

    A subject at the limit on one day and far below it the next constrains the
    decay only from above: nothing in daily-resolution data distinguishes "fell
    tenfold overnight" from "fell a millionfold overnight", so the optimizer runs
    a0 up without penalty. This is the runaway that _MIN_HALF_LIFE_DAYS catches.
    """
    cliff = [(1, 1e8), (2, 1e7), (3, "negative"), (4, "negative"), (5, "negative")]
    normal = [(t, 1e7 * np.exp(-0.5 * t)) for t in range(1, 11)]

    def subject(points):
        return {
            "measurements": [
                {"analyte": "stool", "time": t, "value": v} for t, v in points
            ]
        }

    return {
        "dataset_id": "cliff_study",
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
        "participants": (
            [subject(cliff) for _ in range(n_cliff)]
            + [subject(normal) for _ in range(n_normal)]
        ),
    }


def test_runaway_decay_is_flagged_degenerate():
    """The upper end of the check, symmetric with the collapse floor.

    Before this existed the check was asymmetric: a physically-sited floor of
    1e-2 against a raw optimizer bound of exp(25), so subjects with a0 of 84 to
    142 -- half-lives of five to eight minutes -- passed as ordinary estimates
    and pushed population peaks to 10^18 gc/mL.
    """
    dataset = _cliff_dataset(n_cliff=1, n_normal=4)
    with pytest.warns(UserWarning, match="implied half-life"):
        fit = fit_shedding_model(dataset, analyte="stool", model="exponential")

    assert fit.n_degenerate_subjects == 1
    flagged = fit.subject_params[fit.subject_params["degenerate"]]
    assert (np.log(2.0) / flagged["a0"] < _MIN_HALF_LIFE_DAYS).all()
    # And the runaway is kept out of the population it would otherwise dominate.
    assert fit.half_life_days > _MIN_HALF_LIFE_DAYS


def test_runaway_subject_is_retained_for_inspection():
    """Flagged, not deleted -- same contract as a collapsed subject."""
    dataset = _cliff_dataset(n_cliff=1, n_normal=4)
    with pytest.warns(UserWarning):
        fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert len(fit.subject_params) == 5
    assert fit.n_subjects == 5


def test_ordinary_fast_decay_is_not_flagged():
    """The threshold must not reach into legitimately fast but resolvable decay.

    A half-life of 0.5 days is four times _MIN_HALF_LIFE_DAYS and is perfectly
    estimable from daily samples, so nothing here may be flagged.
    """
    dataset = _flat_dataset(n_flat=0, n_decaying=3)
    for participant in dataset["participants"]:
        for measurement in participant["measurements"]:
            measurement["value"] = 1e9 * np.exp(
                -np.log(2.0) / 0.5 * measurement["time"]
            )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert fit.n_degenerate_subjects == 0
    assert fit.half_life_days == pytest.approx(0.5, abs=0.05)


def test_degenerate_check_locates_the_decay_parameter_by_name():
    """Both models must apply the half-life test to a0, wherever it sits."""
    for model, k in (("exponential", 2), ("gamma", 3)):
        index = PARAM_NAMES[model].index("a0")
        theta = np.zeros((1, k))
        theta[0, index] = np.log(np.log(2.0) / (_MIN_HALF_LIFE_DAYS / 2))
        assert _degenerate_subjects(theta, model)[0]
        # A decay comfortably inside the threshold is not flagged.
        theta[0, index] = np.log(np.log(2.0) / (_MIN_HALF_LIFE_DAYS * 10))
        assert not _degenerate_subjects(theta, model)[0]


def _rise_dataset(n_rising, n_falling, n_short=0):
    """Subjects that peak mid-window, subjects that peak at their first reading.

    The rise gate counts subjects, not measurements, so these are built one
    subject at a time with an unambiguous shape each.

    Each subject gets freshly-built measurement dicts. Handing every subject a
    copy of one shared list would let a caller that edits one subject silently
    edit them all.
    """

    def subject(times, values):
        return {
            "measurements": [
                {"analyte": "stool", "time": t, "value": v}
                for t, v in zip(times, values)
            ]
        }

    def rising():
        return subject([1, 2, 3, 4, 5], [1e3, 1e5, 1e7, 1e5, 1e3])

    def falling():
        return subject([1, 2, 3, 4, 5], [1e7, 1e6, 1e5, 1e4, 1e3])

    def short():
        # Three usable measurements, so prepare_observations keeps this subject,
        # but only two are uncensored -- below what the gate needs to judge a
        # shape. Deliberately built to *look* like a rise, so a gate that forgot
        # its minimum would count it.
        return subject([1, 2, 3], [1e3, 1e7, "negative"])

    return {
        "dataset_id": "rise_study",
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
        "participants": (
            [rising() for _ in range(n_rising)]
            + [falling() for _ in range(n_falling)]
            + [short() for _ in range(n_short)]
        ),
    }


def test_gamma_refused_when_most_subjects_never_rise():
    with pytest.raises(SheddingDataError) as excinfo:
        fit_shedding_model(_rise_dataset(1, 3), analyte="stool", model="gamma")
    assert excinfo.value.reason == "no_rise_observed"
    assert "exponential" in str(excinfo.value)


def test_gamma_allowed_at_exactly_the_threshold():
    """The gate is >= 0.5, so a half-and-half analyte fits rather than refuses."""
    fit = fit_shedding_model(_rise_dataset(2, 2), analyte="stool", model="gamma")
    assert fit.pct_subjects_with_rise == pytest.approx(50.0)


def test_exponential_is_never_gated_on_rise():
    """Post-peak sampling is what the exponential model is for."""
    fit = fit_shedding_model(_rise_dataset(0, 4), analyte="stool", model="exponential")
    assert fit.pct_subjects_with_rise == pytest.approx(0.0)


def test_rise_gate_ignores_subjects_with_too_few_readings():
    """Two rising readings are not evidence of a rise; they are just two points.

    Without the minimum, the three short subjects here would each read as a
    rise and drag a 1-in-4 analyte to 4-in-7, over the threshold.
    """
    with pytest.raises(SheddingDataError) as excinfo:
        fit_shedding_model(
            _rise_dataset(1, 3, n_short=3), analyte="stool", model="gamma"
        )
    assert excinfo.value.reason == "no_rise_observed"


def test_rise_fraction_is_nan_when_no_subject_can_be_judged():
    observations = prepare_observations(
        _rise_dataset(0, 0, n_short=4), "stool", "gamma"
    )
    assert math.isnan(_fraction_observing_a_rise(observations))


def test_gamma_refused_when_no_subject_can_be_judged():
    """NaN must refuse, not slip through the comparison as neither < nor >=."""
    with pytest.raises(SheddingDataError) as excinfo:
        fit_shedding_model(
            _rise_dataset(0, 0, n_short=4), analyte="stool", model="gamma"
        )
    assert excinfo.value.reason == "no_rise_observed"
    assert "no subject had enough readings" in str(excinfo.value)


def test_two_subject_exponential_fit_is_refused_as_a_population():
    """Three subjects are the minimum for a 2-parameter between-subject covariance.

    With two, the covariance is rank-deficient and the "population" is an
    artefact of averaging too little. The repository's two-subject gamma fit
    reported a peak of 10^109 gc/mL before this gate existed.
    """
    dataset = _flat_dataset(n_flat=0, n_decaying=2)
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert fit.n_subjects == 2  # the fit itself is allowed
    with pytest.raises(SheddingDataError) as excinfo:
        require_estimable_population(fit)
    assert excinfo.value.reason == "too_few_subjects_for_population"


def test_three_subject_exponential_fit_is_accepted_as_a_population():
    """One more subject than parameters is enough; the gate is not stricter."""
    dataset = _flat_dataset(n_flat=0, n_decaying=3)
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert fit.n_subjects == 3
    require_estimable_population(fit)  # must not raise


def test_population_gate_counts_after_degenerate_exclusion():
    """What matters is the subjects feeding mu/Sigma, not the subjects fitted."""
    dataset = _flat_dataset(n_flat=1, n_decaying=3)
    with pytest.warns(UserWarning, match="collapsed onto the bounds"):
        fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    # 4 fitted, 1 degenerate -> 3 retained, which is exactly enough.
    assert (fit.n_subjects, fit.n_degenerate_subjects) == (4, 1)
    require_estimable_population(fit)

    # The same arithmetic one subject lower must refuse, and say so.
    dataset = _flat_dataset(n_flat=1, n_decaying=2)
    with pytest.warns(UserWarning, match="collapsed onto the bounds"):
        fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert (fit.n_subjects, fit.n_degenerate_subjects) == (3, 1)
    with pytest.raises(SheddingDataError, match="excluded as degenerate"):
        require_estimable_population(fit)


def test_single_subject_fit_is_still_allowed_by_the_fitter():
    """The gate must not live in fit_shedding_model.

    tests/test_shedding_tutorial_agreement.py validates the port against the
    published Rstan posterior by fitting woelfel subject 3 alone, through
    fit_shedding_model directly. A one-subject fit is honest about itself -- it
    carries a zero covariance -- so it is the *catalog* that must refuse it, not
    the fitter.
    """
    dataset = _flat_dataset(n_flat=0, n_decaying=1)
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert fit.n_subjects == 1
    with pytest.raises(SheddingDataError):
        require_estimable_population(fit)


def test_median_first_observed_day_is_the_median_not_the_minimum():
    """One early-enrolled subject must not make a late study look well-observed.

    This is the distinction the column's name now carries: with subjects
    starting on days 1, 21 and 31, the minimum would report 1.0 and imply the
    study saw the reference event, while the median reports 21.0 and correctly
    says the typical subject was first sampled three weeks in.
    """
    dataset = _flat_dataset(n_flat=0, n_decaying=3)
    # _flat_dataset samples days 1..10 for every subject; stagger two of them.
    for offset, participant in zip((20, 30), dataset["participants"][1:]):
        for measurement in participant["measurements"]:
            measurement["time"] += offset

    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert fit.median_first_observed_day == pytest.approx(21.0)

    for participant in dataset["participants"]:
        for measurement in participant["measurements"]:
            measurement["time"] += 5
    later = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert later.median_first_observed_day == pytest.approx(26.0)


def test_median_first_observed_day_ignores_degenerate_subjects():
    """It must describe the subjects behind mu/Sigma, not the ones discarded."""
    dataset = _flat_dataset(n_flat=1, n_decaying=3)
    # The flat (degenerate) subject is the first participant; start it far later
    # than the rest, so including it would drag the median upward.
    for measurement in dataset["participants"][0]["measurements"]:
        measurement["time"] += 40

    with pytest.warns(UserWarning, match="collapsed onto the bounds"):
        fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert fit.n_degenerate_subjects == 1
    # The three retained subjects all start on day 1; the discarded one on 41.
    assert fit.median_first_observed_day == pytest.approx(1.0)


def _gamma_dataset_sampled_before_detectability(a0=0.5, b0=2.0, c0=12.0):
    """A study whose first sampling day sits below the limit of quantification.

    Sampled from an exact gamma curve so the censored point is consistent data
    rather than a contradiction: at day 0.5 the curve is 4.50 log10, under a
    limit of 4.8, and it clears the limit from day 1 onward. Constructing it any
    other way does not work -- a pure decay is highest at its first reading, so
    an early `negative` would contradict the very curve being fitted.
    """
    limit_log10 = 4.8
    times = [0.5, 1, 2, 3, 4, 6, 8, 10]

    def log10_at(time):
        return (c0 + b0 * np.log(time) - a0 * time) / np.log(10.0)

    def subject(offset):
        measurements = []
        for time in times:
            value = 10 ** (log10_at(time) + offset)
            measurements.append(
                {
                    "analyte": "stool",
                    "time": time,
                    "value": (
                        "negative" if value < 10 ** limit_log10 else float(value)
                    ),
                }
            )
        return {"measurements": measurements}

    return {
        "dataset_id": "late_start_study",
        "analytes": {
            "stool": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "limit_of_quantification": 10**limit_log10,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [subject(o) for o in (-0.1, -0.05, 0.0, 0.05, 0.1)],
    }


def test_median_first_observed_day_counts_censored_observations():
    """A `negative` still means the study looked, and still constrains the curve."""
    dataset = _gamma_dataset_sampled_before_detectability()
    # Guard the fixture: the earliest measurement must actually be censored, or
    # this test would pass without testing anything.
    earliest = min(
        (m for p in dataset["participants"] for m in p["measurements"]),
        key=lambda m: m["time"],
    )
    assert earliest["time"] == 0.5 and earliest["value"] == "negative"

    fit = fit_shedding_model(dataset, analyte="stool", model="gamma")
    assert fit.n_degenerate_subjects == 0
    assert fit.median_first_observed_day == pytest.approx(0.5)


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


def test_to_dict_from_dict_round_trip():
    """SheddingFit.to_dict()/from_dict() is the fit-level persistence story:

    an ABM should be able to fit once, save to YAML, and simulate across many
    runs without refitting. subject_params is deliberately not part of the
    round trip -- everything needed to simulate is.
    """
    mean = np.array([np.log(0.6), np.log(18.0)])
    cov = np.diag([0.04, 0.04])
    fit = _minimal_fit(mean, cov)
    fit.n_degenerate_subjects = 2
    fit.pct_subjects_with_rise = 62.5
    fit.median_first_observed_day = 3.0

    payload = fit.to_dict()
    restored = SheddingFit.from_dict(payload)

    np.testing.assert_allclose(restored.population_mean, fit.population_mean)
    np.testing.assert_allclose(restored.population_cov, fit.population_cov)
    assert restored.sigma == pytest.approx(fit.sigma)
    assert restored.censoring_limit == pytest.approx(fit.censoring_limit)
    assert restored.model == fit.model
    assert restored.method == fit.method
    assert restored.dataset_id == fit.dataset_id
    assert restored.analyte == fit.analyte
    assert restored.biomarker == fit.biomarker
    assert restored.specimen == fit.specimen
    assert restored.reference_event == fit.reference_event
    assert restored.unit == fit.unit
    assert restored.gene_target == fit.gene_target
    assert restored.dose == fit.dose
    assert restored.vaccine_type == fit.vaccine_type
    assert restored.n_subjects == fit.n_subjects
    assert restored.n_degenerate_subjects == fit.n_degenerate_subjects
    assert restored.pct_subjects_with_rise == pytest.approx(fit.pct_subjects_with_rise)
    assert restored.median_first_observed_day == pytest.approx(
        fit.median_first_observed_day
    )
    assert restored.converged == fit.converged
    assert restored.subject_params is None


def test_deserialized_fit_can_still_simulate():
    """The property that lets an ABM fit once and simulate across many runs."""
    from shedding_hub.shedding_simulate import simulate_shedding

    mean = np.array([np.log(0.6), np.log(18.0)])
    fit = _minimal_fit(mean, np.diag([0.04, 0.04]))
    restored = SheddingFit.from_dict(fit.to_dict())
    traj = simulate_shedding(restored, n_individuals=5, times=[1.0, 2.0], seed=0)
    assert len(traj) == 10


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
from shedding_hub.shedding_models import log10_concentration_rowwise, peak_day

DATA = pathlib.Path(__file__).parent.parent / "data"

# Mirrors the historical parameter floor described in shedding_fit.py's
# _DEGENERATE_PARAM comment. Not imported from there because nothing in that
# module reads it directly any more -- it survives only as documentation of
# where _DEGENERATE_PARAM (1e-2) sits relative to it.
_PARAM_FLOOR = 1e-6


@pytest.fixture(scope="module")
def woelfel():
    """The real SARS-CoV-2 study the collapse was first observed on."""
    return sh.load_dataset("woelfel2020virological", local=str(DATA))


# Woelfel sputum under the gamma model is deliberately absent: only 44% of its
# subjects observe a rise, so the rise gate refuses it. That refusal is asserted
# directly by test_woelfel_sputum_gamma_is_refused_for_lack_of_a_rise.
WOELFEL_FITTABLE = [
    ("stool", "exponential"),
    ("stool", "gamma"),
    ("sputum", "exponential"),
]


@pytest.fixture(scope="module")
def woelfel_fits(woelfel):
    """Every analyte/model combination this bug affected that is still fitted."""
    return {
        (analyte, model): fit_shedding_model(woelfel, analyte=analyte, model=model)
        for analyte, model in WOELFEL_FITTABLE
    }


def test_woelfel_sputum_gamma_is_refused_for_lack_of_a_rise(woelfel):
    """
    Real-data guard on the rise gate, on the analyte that motivated it.

    Woelfel sputum was sampled from symptom onset onward and 5 of its 9
    adequately-sampled subjects peak at their first reading. Fitting a rise
    parameter to that leaves b0 unidentifiable, so the gamma model is refused
    outright rather than reported with an arbitrary b0. The exponential fit of
    the same analyte is unaffected and is covered by the tests below.
    """
    with pytest.raises(SheddingDataError) as excinfo:
        fit_shedding_model(woelfel, analyte="sputum", model="gamma")
    assert excinfo.value.reason == "no_rise_observed"

    exponential = fit_shedding_model(woelfel, analyte="sputum", model="exponential")
    assert exponential.pct_subjects_with_rise == pytest.approx(100 * 4 / 9)


def test_woelfel_stool_passes_the_rise_gate(woelfel_fits):
    """Stool clears the gate at exactly the threshold, so pin it deliberately."""
    fit = woelfel_fits[("stool", "gamma")]
    assert fit.pct_subjects_with_rise == pytest.approx(50.0)


@pytest.mark.parametrize("analyte,model", WOELFEL_FITTABLE)
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


@pytest.mark.parametrize("analyte,model", WOELFEL_FITTABLE)
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
    assert smallest > 10 * _PARAM_FLOOR, (
        f"{analyte}/{model} has a subject parameter at {smallest:.3g}, within a "
        f"factor of ten of the floor {_PARAM_FLOOR:.0e}:\n{fit.subject_params}"
    )


def _subject_peak_log10(fit):
    """Each subject's own peak, on the same scale as ``fit.peak_log10``."""
    params = fit.subject_params[list(fit.param_names)].to_numpy(dtype=float)
    peaks = peak_day(fit.model, params)
    return log10_concentration_rowwise(fit.model, params, peaks[:, None])[:, 0]


def test_woelfel_stool_gamma_peak_lies_within_its_subjects_range(woelfel_fits):
    """
    Real-data guard: the population peak must resemble the population.

    Evaluating peak_log10 at exp(population_mean) reported 2.218 for this fit,
    below the peak of every subject but one. That is a coordinate-wise-average
    artefact: at the peak the quantity is
    (c0 + b0*(ln b0 - ln a0) - b0) / ln(10), a nonlinear function of all three
    parameters, so its population median is not the function evaluated at the
    median parameters -- and b0 and c0 trade off along a ridge, so the averaged
    parameter vector belongs to no subject. Taking the median over draws from
    MVN(mu, Sigma) instead puts it among the subjects it is meant to summarize
    (about 3.2, against a subject median near 4.2, with negatives censored at the
    declared LOQ of log10 2.0).

    Pinned as an interval rather than a value: the point of the change is that
    the summary sits among the subjects, not that it equals any number.
    """
    fit = woelfel_fits[("stool", "gamma")]
    subject_peaks = _subject_peak_log10(fit)
    assert subject_peaks.min() < fit.peak_log10 < subject_peaks.max()
    # Also comfortably nearer the subjects' own median than the coordinate-average
    # artefact managed: that put the peak ~2.1 below the median, so a 1.5 bound
    # cleanly separates a sound summary from that failure while tolerating the
    # small shifts a methodology change (e.g. the censoring limit) legitimately
    # produces.
    assert abs(fit.peak_log10 - np.median(subject_peaks)) < 1.5


@pytest.mark.parametrize("analyte,model", WOELFEL_FITTABLE)
def test_woelfel_peak_lies_within_its_subjects_range(woelfel_fits, analyte, model):
    """The same containment property, required of every real fit."""
    fit = woelfel_fits[(analyte, model)]
    subject_peaks = _subject_peak_log10(fit)
    assert subject_peaks.min() < fit.peak_log10 < subject_peaks.max()


def test_exponential_peak_log10_matches_its_closed_form(make_synthetic_dataset):
    """
    Sampling must not disturb the case that has an exact answer.

    For the exponential model the peak is at t = 0, so peak_log10 is
    c0 / ln(10) -- a monotone transform of a single lognormal, whose median is
    exactly the value at exp(mu). Sampling is used anyway to keep one code path
    for both models, so the two must agree; the tolerance is Monte Carlo error
    at _PEAK_LOG10_DRAWS draws, not slack for a real discrepancy.
    """
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=30
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    closed_form = fit.median_params[1] / np.log(10.0)
    assert fit.peak_log10 == pytest.approx(closed_form, abs=0.05)


def test_peak_log10_is_deterministic(woelfel_fits):
    """A catalog rebuild must reproduce byte-for-byte, so the seed is fixed."""
    fit = woelfel_fits[("stool", "gamma")]
    assert fit.peak_log10 == fit.peak_log10


def test_peak_log10_with_zero_covariance_is_the_point_value():
    """A single-subject fit has no spread, so sampling must be a no-op."""
    mean = np.array([np.log(0.5), np.log(2.0), np.log(12.0)])
    fit = _minimal_fit(mean, np.zeros((3, 3)), model="gamma")
    params = np.exp(mean)[None, :]
    expected = log10_concentration_rowwise(
        "gamma", params, peak_day("gamma", params)[:, None]
    )[0, 0]
    assert fit.peak_log10 == pytest.approx(expected)
