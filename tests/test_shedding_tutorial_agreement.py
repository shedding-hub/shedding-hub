"""
Validate the Python port against the published Rstan tutorial.

Reference: https://shedding-hub.github.io/tutorials/Bayesian-workflow-Rstan.html
Subject 3 of woelfel2020virological, stool, censored exponential model, reported
posterior means a0 = 0.74, c0 = 20.37, sig_obs = 0.92 under flat priors.

These tolerances are the only external check that the Python port matches the
published R workflow. If this test fails, the port is wrong: investigate,
don't loosen the numbers. Specific tells worth checking first: c0 off by
roughly a factor of ln(10) means the natural-log scale is being mishandled,
and a0 near 0.55 rather than 0.74 means censored points are being dropped
instead of entering the likelihood (0.55 is the tutorial's own uncensored
result).
"""

import matplotlib

matplotlib.use("Agg")

import pathlib

import pytest

import shedding_hub as sh
from shedding_hub.shedding_fit import fit_shedding_model

DATA = pathlib.Path(__file__).parent.parent / "data"


@pytest.fixture
def woelfel_subject_3():
    dataset = sh.load_dataset("woelfel2020virological", local=str(DATA))
    dataset["participants"] = [dataset["participants"][2]]
    return dataset


def test_subject_3_data_matches_the_tutorial(woelfel_subject_3):
    """Guard the fixture: if the dataset changes, the comparison is void."""
    measurements = [
        m
        for m in woelfel_subject_3["participants"][0]["measurements"]
        if m["analyte"] == "stool"
    ]
    positives = [m for m in measurements if m["value"] != "negative"]
    negatives = [m for m in measurements if m["value"] == "negative"]
    assert len(positives) == 14
    assert sorted(m["time"] for m in negatives) == [20, 22, 23]
    assert woelfel_subject_3["analytes"]["stool"]["limit_of_quantification"] == 100


def test_exponential_fit_agrees_with_published_posterior(woelfel_subject_3):
    fit = fit_shedding_model(woelfel_subject_3, analyte="stool", model="exponential")
    a0, c0 = fit.median_params
    assert a0 == pytest.approx(0.74, abs=0.15)
    assert c0 == pytest.approx(20.37, abs=2.0)
    # sigma's tolerance is proportionally the widest, and the fitted value is
    # expected to land below 0.92: we're comparing an MLE to a posterior mean,
    # and two effects both push the posterior mean of a scale parameter above
    # the MLE — ML variance divides by n rather than n - k, and the posterior
    # for a positive scale parameter is right-skewed so its mean exceeds its
    # mode. A lower sigma here is the expected direction, not a defect.
    assert fit.sigma == pytest.approx(0.92, abs=0.3)
    assert fit.censoring_limit == pytest.approx(2.0)
