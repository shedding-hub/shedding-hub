"""Recovery and invariance checks for cycle-threshold fitting."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from shedding_hub import fit_shedding_model
from shedding_hub.shedding_fit import CT_REFERENCE
from shedding_hub.shedding_models import log10_concentration


def _synthetic(value_type, a0=0.25, b0=1.5, c0=12.0, n_subjects=25, seed=0):
    """
    One cohort, emitted either as concentrations or as the Ct values that the
    SAME curve would have produced through an assay with slope 3.5 and
    intercept 38. Peak time is b0 / a0 = 6.0 days for both.
    """
    rng = np.random.default_rng(seed)
    times = np.array([1.0, 2.0, 4.0, 6.0, 8.0, 11.0, 15.0, 20.0])
    participants = []
    for _ in range(n_subjects):
        params = np.array([[a0, b0, c0 + rng.normal(0.0, 0.8)]])
        y = log10_concentration("gamma", params, times)[0]
        y = y + rng.normal(0.0, 0.2, y.size)
        if value_type == "ct":
            values = 38.0 - 3.5 * y
        else:
            values = 10.0**y
        participants.append(
            {
                "measurements": [
                    {"analyte": "a", "time": float(t), "value": float(v)}
                    for t, v in zip(times, values)
                ]
            }
        )
    unit = "cycle threshold" if value_type == "ct" else "gc/mL"
    return {
        "dataset_id": f"synthetic_{value_type}",
        "analytes": {
            "a": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": unit,
                "limit_of_detection": 40 if value_type == "ct" else 10,
                "limit_of_quantification": "unknown",
            }
        },
        "participants": participants,
    }


def test_ct_fit_recovers_the_true_peak_time():
    fit = fit_shedding_model(_synthetic("ct"), analyte="a", model="gamma")
    assert fit.peak_day == pytest.approx(6.0, rel=0.25)


def test_ct_and_concentration_agree_on_peak_time():
    # The central claim: the standard-curve slope cancels in b0/a0, so both
    # scales must land on the same peak day even though the assay applied a
    # slope of 3.5 and an intercept of 38 to one of them.
    ct = fit_shedding_model(_synthetic("ct"), analyte="a", model="gamma")
    conc = fit_shedding_model(_synthetic("concentration"), analyte="a", model="gamma")
    assert ct.peak_day == pytest.approx(conc.peak_day, rel=0.15)


def test_ct_decay_rate_carries_the_assay_slope():
    # The other half of the claim, stated as a fact rather than a hope: a0 is
    # NOT invariant. It comes back scaled by the assay slope of 3.5, which is
    # exactly why half_life_days is excluded from comparable_with.
    ct = fit_shedding_model(_synthetic("ct"), analyte="a", model="gamma")
    conc = fit_shedding_model(_synthetic("concentration"), analyte="a", model="gamma")
    a0_ct = float(ct.subject_params["a0"].median())
    a0_conc = float(conc.subject_params["a0"].median())
    assert a0_ct / a0_conc == pytest.approx(3.5, rel=0.35)


def test_a_ct_fit_is_a_peak_not_a_trough():
    # If the sign flip were dropped anywhere in the chain, the optimizer would
    # still converge -- on an inverted curve. Assert the fitted median rises
    # from day 1 to the peak and falls after it.
    fit = fit_shedding_model(_synthetic("ct"), analyte="a", model="gamma")
    params = fit.median_params
    heights = log10_concentration("gamma", params[None, :], np.array([1.0, 6.0, 20.0]))[
        0
    ]
    assert heights[1] > heights[0]
    assert heights[1] > heights[2]
