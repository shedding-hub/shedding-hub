import pytest

from shedding_hub.shedding_select import (
    REFERENCE_EVENT_CLASSES,
    classify_reference_event,
)


@pytest.mark.parametrize(
    "event, expected",
    [
        ("inoculation", "exposure"),
        ("vaccination", "exposure"),
        ("symptom onset", "landmark"),
        ("enrollment", "administrative"),
        ("confirmation date", "administrative"),
        ("hospital admission", "administrative"),
        ("treatment", "administrative"),
    ],
)
def test_each_known_event_has_its_class(event, expected):
    assert classify_reference_event(event) == expected


def test_unknown_event_is_administrative():
    """The conservative default: a new event must not silently acquire a clock."""
    assert classify_reference_event("full moon") == "administrative"
    assert classify_reference_event(None) == "administrative"


def test_every_shipped_reference_event_is_classified_deliberately(shipped_catalog):
    """
    Staleness check, in the shape of test_shipped_catalog_covers_every_dataset.

    A dataset introducing a new reference event must be a decision, not a silent
    fall through to 'administrative'.
    """
    events = set(shipped_catalog.table["reference_event"].dropna())
    missing = events - set(REFERENCE_EVENT_CLASSES)
    assert not missing, (
        f"Reference event(s) {sorted(missing)} appear in the shipped catalog but "
        "are not classified in REFERENCE_EVENT_CLASSES. Add them deliberately: "
        "'exposure' if the event is the exposure itself, 'landmark' if it has a "
        "defined offset from infection, 'administrative' otherwise."
    )


def test_options_lists_groups_best_first(shipped_catalog):
    from shedding_hub.shedding_select import shedding_options

    options = shedding_options(
        catalog=shipped_catalog, biomarker="SARS-CoV-2", specimen="stool"
    )
    assert len(options) > 1
    assert list(options["rank"]) == list(range(1, len(options) + 1))
    # Rule 1: a defensible clock outranks a better curve.
    assert options.iloc[0]["event_class"] in ("exposure", "landmark")


def test_options_prefers_the_rise_capable_model_within_a_clock(shipped_catalog):
    """Rule 3, held apart from rules 1 and 2 by fixing the reference event."""
    from shedding_hub.shedding_select import shedding_options

    options = shedding_options(
        catalog=shipped_catalog,
        biomarker="SARS-CoV-2",
        specimen="stool",
        reference_event="symptom onset",
    )
    order = list(options["model"])
    assert order.index("gamma") < order.index("exponential")


def test_options_counts_one_study_once(shipped_catalog):
    """
    natarajan contributes 14 exponential fits to one group, one per assay.

    Counting fits rather than studies would report a 14-study group that
    make_ensemble refuses to build.
    """
    from shedding_hub.shedding_select import shedding_options

    options = shedding_options(
        catalog=shipped_catalog,
        biomarker="SARS-CoV-2",
        specimen="stool",
        reference_event="enrollment",
        unit="gc/mL",
        model="exponential",
    )
    assert len(options) == 1
    assert options.iloc[0]["n_studies"] == 1


def test_sars_cov_2_stool_ranking_is_pinned(shipped_catalog):
    """Pin the shipped outcome so a rule change shows in a diff, not silently."""
    from shedding_hub.shedding_select import shedding_options

    options = shedding_options(
        catalog=shipped_catalog, biomarker="SARS-CoV-2", specimen="stool"
    )
    top = options.head(2)
    assert list(top["reference_event"]) == ["symptom onset", "symptom onset"]
    assert list(top["unit"]) == ["gc/mL", "gc/mL"]
    assert list(top["model"]) == ["gamma", "exponential"]
    # Rule 4 losing to rule 3 within a settled unit, which is intended: the
    # 2-study gamma is preferred to the 3-study exponential.
    assert list(top["n_studies"]) == [2, 3]
    # Within one event class, every gc/dry gram group ranks below every gc/mL
    # one. Deliberately not asserted globally: rule 1 puts a landmark
    # gc/dry gram group (rank 3) above an administrative gc/mL one (rank 6),
    # which is the clock outranking the unit, exactly as intended.
    landmark = options[options["event_class"] == "landmark"]
    assert landmark[landmark["unit"] == "gc/mL"]["rank"].max() < (
        landmark[landmark["unit"] == "gc/dry gram"]["rank"].min()
    )


def test_a_better_represented_unit_outranks_a_richer_model(shipped_catalog):
    """
    Rule 2 in isolation: without it, unit is settled as a side effect.

    gc/mL is reported by 4 SARS-CoV-2 stool studies and gc/dry gram by 2, so
    gc/mL wins the unit outright -- even though gamma_shifted is identifiable on
    gc/dry gram and not on gc/mL, and would otherwise have taken rank 1.
    """
    from shedding_hub.shedding_select import shedding_options

    options = shedding_options(
        catalog=shipped_catalog, biomarker="SARS-CoV-2", specimen="stool"
    )
    assert options.iloc[0]["unit"] == "gc/mL"
    assert options.iloc[0]["n_unit_studies"] == 4
    shifted = options[options["model"] == "gamma_shifted"].iloc[0]
    assert shifted["unit"] == "gc/dry gram"
    assert shifted["rank"] > options.iloc[0]["rank"]


def test_options_raises_when_nothing_matches(shipped_catalog):
    from shedding_hub.shedding_select import shedding_options

    with pytest.raises(ValueError, match="No fits match"):
        shedding_options(catalog=shipped_catalog, biomarker="not a biomarker")


def test_options_accepts_an_explicit_catalog(make_synthetic_dataset):
    """Not everything is the shipped catalog: private fits must work too."""
    import numpy as np

    from shedding_hub.shedding_catalog import fit_shedding_models
    from shedding_hub.shedding_select import shedding_options

    dataset = make_synthetic_dataset(
        "gamma",
        [0.0, np.log(2.0), np.log(12.0)],
        np.diag([0.04, 0.04, 0.09]),
        n_subjects=12,
        seed=3,
    )
    catalog = fit_shedding_models([dataset], models=("exponential", "gamma"))
    options = shedding_options(catalog=catalog)
    assert set(options["model"]) == {"exponential", "gamma"}
    assert options.iloc[0]["model"] == "gamma"


def test_for_returns_a_simulable_ensemble(shipped_catalog):
    import numpy as np

    from shedding_hub import simulate_shedding
    from shedding_hub.shedding_select import shedding_for

    source = shedding_for("SARS-CoV-2", "stool", catalog=shipped_catalog)
    traj = simulate_shedding(source, n_individuals=20, times=np.arange(1, 11), seed=1)
    assert len(traj) == 200
    # times start at 1, not 0: the gamma models are undefined at their own origin.


def test_for_agrees_with_options_rank_one(shipped_catalog):
    """The property that keeps the two surfaces from drifting apart."""
    from shedding_hub.shedding_select import shedding_for, shedding_options

    best = shedding_options(
        catalog=shipped_catalog, biomarker="SARS-CoV-2", specimen="stool"
    ).iloc[0]
    source = shedding_for("SARS-CoV-2", "stool", catalog=shipped_catalog)
    assert source.selection.picked["model"] == best["model"]
    assert source.selection.picked["reference_event"] == best["reference_event"]
    assert source.selection.picked["unit"] == best["unit"]


def test_for_explains_what_it_passed_over(shipped_catalog):
    from shedding_hub.shedding_select import shedding_for

    selection = shedding_for("SARS-CoV-2", "stool", catalog=shipped_catalog).selection
    assert selection.reason
    assert len(selection.passed_over) > 0
    assert "rank" in selection.passed_over.columns


def test_for_pins_an_overridden_key(shipped_catalog):
    from shedding_hub.shedding_select import shedding_for

    source = shedding_for(
        "SARS-CoV-2", "stool", catalog=shipped_catalog, model="exponential"
    )
    assert source.selection.picked["model"] == "exponential"
    assert source.model == "exponential"


def test_for_records_the_analyte_taken_from_each_study(shipped_catalog):
    """natarajan offers 14 analytes to this group; exactly one must be taken."""
    from shedding_hub.shedding_select import shedding_for

    source = shedding_for(
        "SARS-CoV-2",
        "stool",
        catalog=shipped_catalog,
        reference_event="enrollment",
        model="exponential",
    )
    assert len(source.fits) == 1
    assert set(source.selection.analytes) == {"natarajan2022gastrointestinal"}


def test_every_advertised_group_can_actually_be_built(shipped_catalog):
    """
    The regression guard for the one-analyte-per-study reduction.

    Without it this fails on 13 of the shipped catalog's groups, where one study
    contributes several analytes and make_ensemble refuses the duplicate.
    """
    from shedding_hub.shedding_select import shedding_for, shedding_options

    options = shedding_options(catalog=shipped_catalog)
    for _, row in options.iterrows():
        source = shedding_for(
            row["biomarker"],
            row["specimen"],
            catalog=shipped_catalog,
            reference_event=row["reference_event"],
            unit=row["unit"],
            model=row["model"],
        )
        assert len(source.fits) == row["n_studies"]


def test_reason_does_not_credit_a_rule_that_tied():
    """
    'exposure' and 'landmark' share rank 0, so differing on the label is not
    winning on rule 1. Reading the column instead of the ranked value would
    credit the clock for a decision the alphabetical tie-break actually made.
    """
    import pandas as pd

    from shedding_hub.shedding_select import _reason

    tied = {
        "n_unit_studies": 1,
        "model": "exponential",
        "n_studies": 1,
        "n_subjects": 10,
        "n_measurements": 50,
    }
    best = pd.Series({"event_class": "landmark", **tied})
    runner_up = pd.Series({"event_class": "exposure", **tied})
    assert (
        _reason(best, runner_up) == "it sorted first among otherwise equal candidates"
    )

    # A genuine rule-1 win is still reported as one.
    administrative = pd.Series({"event_class": "administrative", **tied})
    assert _reason(best, administrative) == (
        "its reference event supports an infection time origin"
    )


def test_for_is_deterministic(shipped_catalog):
    from shedding_hub.shedding_select import shedding_for

    first = shedding_for(
        "SARS-CoV-2", "stool", catalog=shipped_catalog
    ).selection.picked
    second = shedding_for(
        "SARS-CoV-2", "stool", catalog=shipped_catalog
    ).selection.picked
    assert first == second


def test_single_component_matches_the_bare_fit(make_synthetic_dataset):
    """A one-study answer must simulate exactly as the fit would."""
    import numpy as np
    import pandas as pd

    from shedding_hub import simulate_shedding
    from shedding_hub.shedding_catalog import fit_shedding_models
    from shedding_hub.shedding_select import shedding_for

    dataset = make_synthetic_dataset(
        "gamma",
        [0.0, np.log(2.0), np.log(12.0)],
        np.diag([0.04, 0.04, 0.09]),
        n_subjects=12,
        seed=3,
    )
    catalog = fit_shedding_models([dataset], models=("gamma",))
    source = shedding_for(catalog=catalog)
    times = np.arange(0, 8)
    pd.testing.assert_frame_equal(
        simulate_shedding(source, n_individuals=10, times=times, seed=7),
        simulate_shedding(catalog.fits[0], n_individuals=10, times=times, seed=7),
        check_like=True,
    )
