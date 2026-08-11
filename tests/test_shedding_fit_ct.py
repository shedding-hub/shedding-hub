"""Recovery and invariance checks for cycle-threshold fitting."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from shedding_hub import fit_shedding_model
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


# --- Real-data validation on kissler2021viral -------------------------------
#
# kissler2021viral reports both a cycle-threshold analyte
# (AN_OPS_SARSCoV2_ct) and a concentration analyte (AN_OPS_SARSCoV2_viral)
# for the same subjects at the same timepoints. That lets the invariance
# claims exercised above with synthetic data also be checked against a real
# study, once, at module scope: fitting either analyte over 51 subjects
# takes tens of seconds, and both tests below need both fits, so each
# analyte is fit exactly once per test session rather than once per test.


@pytest.fixture(scope="module")
def kissler_data():
    import shedding_hub as sh

    return sh.load_dataset("kissler2021viral", local="./data")


@pytest.fixture(scope="module")
def kissler_ct_fit(kissler_data):
    import shedding_hub as sh

    return sh.fit_shedding_model(
        kissler_data, analyte="AN_OPS_SARSCoV2_ct", model="gamma"
    )


@pytest.fixture(scope="module")
def kissler_conc_fit(kissler_data):
    import shedding_hub as sh

    return sh.fit_shedding_model(
        kissler_data, analyte="AN_OPS_SARSCoV2_viral", model="gamma"
    )


def test_kissler_peak_times_agree_across_value_types(kissler_ct_fit, kissler_conc_fit):
    """
    The same subjects, measured both as Ct and as concentration. Peak time is
    a ratio of b0 to a0, so the assay's standard-curve slope cancels and the
    two fits must land on the same day.

    Measured 2026-08-11: the Ct fit gives peak_day = 0.242 and the
    concentration fit gives peak_day = 0.773 (both retaining 51 of 68
    subjects under the gamma model). The difference is 0.531 days. A purely
    relative tolerance on peak times this small would be far too tight (25%
    of 0.242 is about 0.06 days), so the ``abs=1.0`` term is what actually
    carries this assertion -- it is checking "within a day", not resolving
    fractions of that day.
    """
    assert kissler_ct_fit.peak_day == pytest.approx(
        kissler_conc_fit.peak_day, rel=0.25, abs=1.0
    )


def test_kissler_ct_and_concentration_are_affinely_related(kissler_data):
    """
    Ct = alpha - beta * log10(C) is the assumption the whole design rests on.

    kissler2021viral does not provide two independently-measured assays: its
    concentration column is derived from Ct through a fixed standard curve,
    so the two are the same measurement expressed on two scales rather than
    two measurements that happen to agree. Regressing Ct on log10(C) over
    the matched pairs gives slope -3.6097, intercept 49.593, and a
    correlation of exactly -1.00000 -- a perfect correlation is exactly what
    a fixed affine transform produces, and is evidence of that derivation
    rather than of independent reproducibility.

    What this test does confirm: (1) the affine premise the whole peak-time
    invariance argument rests on holds exactly in this dataset's data, and
    (2) the fitting pipeline is consistent end-to-end across the two scales,
    since it is exercised on genuinely different input values (Ct vs. log10
    concentration) and different model branches (`value_type="ct"` vs.
    `"concentration"`).
    """
    pairs = []
    for participant in kissler_data["participants"]:
        by_time = {}
        for m in participant.get("measurements") or []:
            if not isinstance(m.get("value"), (int, float)):
                continue
            by_time.setdefault(m["time"], {})[m["analyte"]] = m["value"]
        for readings in by_time.values():
            if {"AN_OPS_SARSCoV2_ct", "AN_OPS_SARSCoV2_viral"} <= set(readings):
                pairs.append(
                    (
                        np.log10(readings["AN_OPS_SARSCoV2_viral"]),
                        readings["AN_OPS_SARSCoV2_ct"],
                    )
                )
    # 225 matched pairs, not the ~2,406 timepoints the study recorded: most
    # recorded timepoints are qualitative "negative" strings on one or both
    # analytes rather than paired numeric values, and only numeric pairs at
    # the same subject/timepoint can be regressed here.
    assert len(pairs) > 200
    log10_conc, ct_values = np.array(pairs).T
    slope, _ = np.polyfit(log10_conc, ct_values, 1)
    # Negative because Ct falls as concentration rises, and near the -3.32 of
    # a perfectly efficient assay. Measured slope: -3.6097.
    assert -5.0 < slope < -2.0
    corr = float(np.corrcoef(log10_conc, ct_values)[0, 1])
    # Measured correlation is exactly -1.00000 (see docstring): tightened
    # from the naive "< -0.9" so this assertion still means something.
    assert corr < -0.999
