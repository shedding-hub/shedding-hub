# Choosing What To Simulate From — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a user ask for "SARS-CoV-2 in stool" and get a simulable, defensibly-chosen shedding source, instead of having to name five catalog keys correctly and getting `No fits match` when they don't.

**Architecture:** One new module, `shedding_hub/shedding_select.py`, holding a reference-event taxonomy and a single ranking implementation exposed through two surfaces: `shedding_options()` (see the choice) and `shedding_for()` (make it). `shedding_for` calls `shedding_options` and takes rank 1, so the two can never disagree. A separate, independent change makes `simulate_shedding` stop claiming an infection time origin it has not earned.

**Tech Stack:** Python 3.10+, numpy, pandas, pytest, black. No new dependencies.

## Global Constraints

- No new install dependencies. numpy, scipy, pandas and pyyaml are already required.
- No change to any fitted number, to `shedding_hub/data/shedding_catalog.yaml`, to `docs/shedding_parameters.{json,csv}`, or to the compatibility keys `make_ensemble` enforces. `tests/test_parameter_export.py` and `test_shipped_catalog_covers_every_dataset` must stay green without regenerating anything.
- `black --check .` must pass. Run `black .` before every commit.
- Reference-event classes are exactly `"exposure"`, `"landmark"`, `"administrative"`. An unrecognised event classifies as `"administrative"`.
- Ranking order is fixed: event class, then unit (by how many studies report it among the candidates sharing its biomarker and specimen), then model, then studies, then subjects, then measurements, then the sorted key tuple. Model order is `gamma_shifted`, `gamma`, `exponential`.
- Every group reported must be buildable: reduce to one fit per study before counting or building.
- Tests must not fit models where a fixture will do, and must not reload the shipped catalog per test — it costs ~2.2s and the suite already takes ~5 minutes. Use the session-scoped `shipped_catalog` fixture (added in Task 1) and the existing `make_synthetic_dataset` fixture, both in `tests/conftest.py`. Never call `load_shedding_catalog()` directly in a test.

---

### Task 1: Reference-event taxonomy

**Files:**
- Create: `shedding_hub/shedding_select.py`
- Modify: `tests/conftest.py`
- Test: `tests/test_shedding_select.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `REFERENCE_EVENT_CLASSES: dict[str, str]`, `classify_reference_event(event: str | None) -> str`. Later tasks rank on the return value, which is one of `"exposure"`, `"landmark"`, `"administrative"`. Also the `shipped_catalog` pytest fixture, which Tasks 2 and 3 use throughout.

- [ ] **Step 1: Add the session-scoped catalog fixture**

Loading the shipped catalog costs about 2.2 seconds, and the tests across Tasks 1-3 need it roughly 97 times — 82 of those inside a single loop. Loading per test would add over three minutes to a suite that currently runs in five. Load it once per session instead.

Append to `tests/conftest.py`:

```python
@pytest.fixture(scope="session")
def shipped_catalog():
    """
    The shipped catalog, loaded once for the whole session.

    Parsing it costs ~2.2s, and the selection tests need it dozens of times.
    Session scope is safe because nothing in these tests mutates a catalog.
    """
    from shedding_hub import load_shedding_catalog

    return load_shedding_catalog()
```

- [ ] **Step 2: Write the failing tests**

Create `tests/test_shedding_select.py`:

```python
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
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `python -m pytest tests/test_shedding_select.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'shedding_hub.shedding_select'`

- [ ] **Step 4: Implement the taxonomy**

Create `shedding_hub/shedding_select.py`:

```python
"""
Choose which fitted estimate to simulate from.

The catalog's five compatibility keys -- biomarker, specimen, reference event,
unit and model -- cut 126 fits into 82 groups, 71 of them a single study. So the
user's problem is not combining many estimates but choosing among incommensurable
ones: gc/mL and gc/dry gram are different quantities, symptom onset and
enrollment are different clocks, and the three models are fitted to different
observation sets and so are not AIC-comparable.

This module makes that choice visible (``shedding_options``) and makes a
documented, overridable default (``shedding_for``). It adds no statistical
assumption: nothing here pools across units or reference events.
"""

# Reference events are not all the same kind of thing, and the difference decides
# whether ``simulate_shedding``'s incubation shift means anything. An 'exposure'
# event *is* the exposure, so the offset from it to infection is zero. A
# 'landmark' has a defined offset -- the incubation period, which the literature
# reports. An 'administrative' date reflects testing behaviour and health-system
# access, and has no fixed relation to infection at all.
REFERENCE_EVENT_CLASSES = {
    "inoculation": "exposure",
    "vaccination": "exposure",
    "symptom onset": "landmark",
    "enrollment": "administrative",
    "confirmation date": "administrative",
    "hospital admission": "administrative",
    "treatment": "administrative",
}


def classify_reference_event(event) -> str:
    """
    Classify a reference event as ``exposure``, ``landmark`` or ``administrative``.

    Unrecognised events -- including ``None`` -- are ``administrative``. That is
    the conservative direction: a reference event the package has never seen has
    not demonstrated a defined offset from infection, so it should not inherit
    one.
    """
    return REFERENCE_EVENT_CLASSES.get(event, "administrative")
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/test_shedding_select.py -q`
Expected: PASS, 9 tests.

- [ ] **Step 6: Format and commit**

```bash
black shedding_hub/shedding_select.py tests/test_shedding_select.py tests/conftest.py
git add shedding_hub/shedding_select.py tests/test_shedding_select.py tests/conftest.py
git commit -m "feat: classify reference events by their relation to infection"
```

---

### Task 2: `shedding_options` — see the choice

**Files:**
- Modify: `shedding_hub/shedding_select.py`
- Test: `tests/test_shedding_select.py`

**Interfaces:**
- Consumes: `classify_reference_event` from Task 1; `load_shedding_catalog`, `SheddingCatalog` from `shedding_hub.shedding_catalog`.
- Produces: `shedding_options(catalog=None, **keys) -> pd.DataFrame` with columns `biomarker, specimen, reference_event, event_class, unit, n_unit_studies, model, n_studies, n_subjects, n_measurements, rank`, sorted best first, `rank` starting at 1. Also `_reduce_to_one_fit_per_study(fits: list) -> list` and `_rank_key(row: dict) -> tuple`, used by Task 3.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_shedding_select.py`:

```python
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
        "gamma", [0.0, np.log(2.0), np.log(12.0)], np.diag([0.04, 0.04, 0.09]),
        n_subjects=12, seed=3,
    )
    catalog = fit_shedding_models([dataset], models=("exponential", "gamma"))
    options = shedding_options(catalog=catalog)
    assert set(options["model"]) == {"exponential", "gamma"}
    assert options.iloc[0]["model"] == "gamma"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_shedding_select.py -q -k options`
Expected: FAIL — `ImportError: cannot import name 'shedding_options'`

- [ ] **Step 3: Implement grouping, reduction and ranking**

Append to `shedding_hub/shedding_select.py`. Add these imports at the top of the file, below the docstring:

```python
from dataclasses import dataclass, field

import pandas as pd

from .shedding_catalog import SheddingCatalog, load_shedding_catalog
```

Then append:

```python
_GROUP_KEYS = ("biomarker", "specimen", "reference_event", "unit", "model")

# Rule 1. 'exposure' and 'landmark' tie: an inoculation is infection, and a
# symptom onset is a fixed, documented offset from it. Both give an agent a
# timeline; an administrative date does not.
_EVENT_CLASS_RANK = {"exposure": 0, "landmark": 0, "administrative": 1}

# Rule 3. Rise-capable first: for a wastewater model the pre-symptomatic rise is
# the epidemiologically interesting part, and an exponential asserts by
# construction that an agent sheds maximally on the day of the reference event. A
# model only reaches the catalog if its gates passed, so presence is already the
# identifiability signal and no separate check is needed here.
#
# Rule 2 -- the unit, by how many studies report it within its own biomarker and
# specimen -- sits between these two and is computed from the candidate set rather
# than being a constant, so it lives on the row as n_unit_studies. Units are
# incommensurable across biomarkers and specimens, and without an explicit rule
# the unit would be settled by whichever one happened to carry the richest model:
# SARS-CoV-2 stool would return a 1-study gc/dry gram gamma_shifted ahead of a
# 3-study gc/mL exponential.
_MODEL_RANK = {"gamma_shifted": 0, "gamma": 1, "exponential": 2}


def _sortable(value) -> str:
    """None and str do not compare, and both appear in these keys."""
    return "" if value is None else str(value)


def _reduce_to_one_fit_per_study(fits: list) -> list:
    """
    Keep one fit per study, so ``make_ensemble`` will accept the group.

    An ensemble may not take two analytes from one study -- that study's subjects
    would enter the mixture twice -- and 13 of the shipped catalog's 82 groups
    contain such a study, one of them contributing 14 fits. This is the narrowing
    ``make_ensemble``'s own error advises, applied by a stated rule: most
    subjects, then most measurements, then the analyte name, which breaks
    remaining ties deterministically.
    """
    by_study = {}
    for fit in fits:
        by_study.setdefault(fit.dataset_id, []).append(fit)
    return [
        sorted(
            group,
            key=lambda f: (-f.n_subjects, -f.n_measurements, _sortable(f.analyte)),
        )[0]
        for _, group in sorted(by_study.items())
    ]


def _rank_key(row) -> tuple:
    """The ranking, in one place. Lower sorts first."""
    return (
        _EVENT_CLASS_RANK[row["event_class"]],
        -row["n_unit_studies"],
        _MODEL_RANK[row["model"]],
        -row["n_studies"],
        -row["n_subjects"],
        -row["n_measurements"],
        tuple(_sortable(row[key]) for key in _GROUP_KEYS),
    )


def _matching_fits(catalog: SheddingCatalog, keys: dict) -> list:
    return [
        fit
        for fit in catalog.fits
        if all(getattr(fit, key, None) == value for key, value in keys.items())
    ]


def _grouped(fits: list) -> dict:
    """Group fits by the five compatibility keys, reduced to one per study."""
    groups = {}
    for fit in fits:
        signature = tuple(getattr(fit, key, None) for key in _GROUP_KEYS)
        groups.setdefault(signature, []).append(fit)
    return {
        signature: _reduce_to_one_fit_per_study(group)
        for signature, group in groups.items()
    }


def shedding_options(catalog: SheddingCatalog | None = None, **keys) -> pd.DataFrame:
    """
    List every group that could be simulated from, best first.

    Each row is one combination of the five keys ``make_ensemble`` requires
    agreement on, so each row can actually be built -- counts are reported after
    reducing each study to a single analyte. ``rank`` 1 is what ``shedding_for``
    returns for the same arguments.

    Args:
        catalog: Catalog to search. Defaults to the shipped one.
        **keys: Attribute filters, e.g. ``biomarker="SARS-CoV-2"``.

    Returns:
        A ``DataFrame`` with one row per group.

    Raises:
        ValueError: If nothing matches ``keys``.
    """
    catalog = load_shedding_catalog() if catalog is None else catalog
    matches = _matching_fits(catalog, keys)
    if not matches:
        raise ValueError(
            f"No fits match {keys}. Browse `catalog.table` for available "
            "combinations, or call shedding_options() with fewer keys."
        )

    # Rule 2's input: how many distinct studies report each unit, counted within
    # one biomarker and specimen. Scoped that way because units are only
    # commensurable there -- counting a unit across the whole candidate set would
    # let gc/mL's SARS-CoV-2 studies decide a rotavirus vaccine row, and would make
    # a group's rank depend on what else the caller happened to leave unfiltered.
    studies_by_unit = {}
    for fit in matches:
        signature = (
            _sortable(fit.biomarker),
            _sortable(fit.specimen),
            _sortable(fit.unit),
        )
        studies_by_unit.setdefault(signature, set()).add(fit.dataset_id)

    rows = []
    for signature, group in _grouped(matches).items():
        row = dict(zip(_GROUP_KEYS, signature))
        row["event_class"] = classify_reference_event(row["reference_event"])
        row["n_unit_studies"] = len(
            studies_by_unit[
                (
                    _sortable(row["biomarker"]),
                    _sortable(row["specimen"]),
                    _sortable(row["unit"]),
                )
            ]
        )
        row["n_studies"] = len(group)
        row["n_subjects"] = sum(fit.n_subjects for fit in group)
        row["n_measurements"] = sum(fit.n_measurements for fit in group)
        rows.append(row)

    rows.sort(key=_rank_key)
    frame = pd.DataFrame(rows)
    frame["rank"] = range(1, len(frame) + 1)
    return frame[
        [
            "biomarker",
            "specimen",
            "reference_event",
            "event_class",
            "unit",
            "n_unit_studies",
            "model",
            "n_studies",
            "n_subjects",
            "n_measurements",
            "rank",
        ]
    ]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_shedding_select.py -q`
Expected: PASS, 16 tests.

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/shedding_select.py tests/test_shedding_select.py
git add shedding_hub/shedding_select.py tests/test_shedding_select.py
git commit -m "feat: list the simulable groups, ranked"
```

---

### Task 3: `shedding_for` — make the choice

**Files:**
- Modify: `shedding_hub/shedding_select.py`
- Modify: `shedding_hub/shedding_ensemble.py` (add the `selection` field to `SheddingEnsemble`)
- Test: `tests/test_shedding_select.py`

**Interfaces:**
- Consumes: `shedding_options`, `_grouped`, `_rank_key`, `_matching_fits` from Task 2; `make_ensemble` from `shedding_hub.shedding_ensemble`.
- Produces: `shedding_for(biomarker=None, specimen=None, *, catalog=None, weights="n_subjects", method="mixture", **keys) -> SheddingEnsemble` whose `.selection` is a `Selection` dataclass with fields `picked: dict`, `passed_over: pd.DataFrame`, `reason: str`, `analytes: dict[str, str]`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_shedding_select.py`:

```python
def test_for_returns_a_simulable_ensemble(shipped_catalog):
    import numpy as np

    from shedding_hub import simulate_shedding
    from shedding_hub.shedding_select import shedding_for

    source = shedding_for("SARS-CoV-2", "stool", catalog=shipped_catalog)
    traj = simulate_shedding(
        source, n_individuals=20, times=np.arange(1, 11), seed=1
    )
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
    assert _reason(best, runner_up) == "it sorted first among otherwise equal candidates"

    # A genuine rule-1 win is still reported as one.
    administrative = pd.Series({"event_class": "administrative", **tied})
    assert _reason(best, administrative) == (
        "its reference event can be placed on an infection timeline"
    )


def test_for_is_deterministic(shipped_catalog):
    from shedding_hub.shedding_select import shedding_for

    first = shedding_for("SARS-CoV-2", "stool", catalog=shipped_catalog).selection.picked
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
        "gamma", [0.0, np.log(2.0), np.log(12.0)], np.diag([0.04, 0.04, 0.09]),
        n_subjects=12, seed=3,
    )
    catalog = fit_shedding_models([dataset], models=("gamma",))
    source = shedding_for(catalog=catalog)
    times = np.arange(0, 8)
    pd.testing.assert_frame_equal(
        simulate_shedding(source, n_individuals=10, times=times, seed=7),
        simulate_shedding(catalog.fits[0], n_individuals=10, times=times, seed=7),
        check_like=True,
    )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_shedding_select.py -q -k "for_ or advertised or single_component"`
Expected: FAIL — `ImportError: cannot import name 'shedding_for'`

- [ ] **Step 3: Add the `selection` field to `SheddingEnsemble`**

In `shedding_hub/shedding_ensemble.py`, change the dataclass fields (currently `fits`, `weights`, `method`) to add a fourth. Find:

```python
@dataclass
class SheddingEnsemble:
    """An ensemble of per-study fits sharing a biomarker, specimen, and unit."""

    fits: list[SheddingFit]
    weights: np.ndarray
    method: str
```

Replace with:

```python
@dataclass
class SheddingEnsemble:
    """An ensemble of per-study fits sharing a biomarker, specimen, and unit."""

    fits: list[SheddingFit]
    weights: np.ndarray
    method: str
    # Provenance, not statistics: set by shedding_for to record why this
    # combination was chosen over the alternatives. Typed loosely to keep
    # shedding_ensemble free of an import from shedding_select, which imports it.
    # Deliberately absent from to_dict: it describes a choice made against one
    # catalog, and would be misleading if restored beside fits from another.
    selection: object = None
```

- [ ] **Step 4: Implement `shedding_for`**

Append to `shedding_hub/shedding_select.py`:

```python
@dataclass
class Selection:
    """Why ``shedding_for`` returned what it did."""

    picked: dict
    passed_over: pd.DataFrame
    reason: str
    analytes: dict = field(default_factory=dict)

    def __str__(self) -> str:
        keys = self.picked
        return (
            f"picked {keys['reference_event']} / {keys['unit']} / {keys['model']} "
            f"({keys['n_studies']} study/studies, {keys['n_subjects']} subjects); "
            f"{self.reason}"
        )


# The ranking rules in order, paired with how to describe the one that decided.
# Each reads the *ranked* quantity, not the raw column, because those differ:
# _EVENT_CLASS_RANK maps 'exposure' and 'landmark' both to 0, so a winner and
# runner-up that differ in event_class may still tie on rule 1. Comparing the
# labels would then credit rule 1 for a decision some later rule actually made.
_RULE_REASONS = (
    (
        lambda row: _EVENT_CLASS_RANK[row["event_class"]],
        "its reference event can be placed on an infection timeline",
    ),
    (lambda row: row["n_unit_studies"], "its unit is reported by more studies"),
    (lambda row: _MODEL_RANK[row["model"]], "its model resolves the rise"),
    (lambda row: row["n_studies"], "it rests on more studies"),
    (lambda row: row["n_subjects"], "it rests on more subjects"),
    (lambda row: row["n_measurements"], "it rests on more measurements"),
)


def _reason(best, runner_up) -> str:
    """Name the first rule on which the winner beat the runner-up."""
    if runner_up is None:
        return "it was the only candidate"
    for ranked_value, description in _RULE_REASONS:
        if ranked_value(best) != ranked_value(runner_up):
            return description
    return "it sorted first among otherwise equal candidates"


def shedding_for(
    biomarker=None,
    specimen=None,
    *,
    catalog: SheddingCatalog | None = None,
    weights="n_subjects",
    method: str = "mixture",
    **keys,
):
    """
    Return the best-ranked simulable source for a biomarker and specimen.

    Equivalent to taking ``rank`` 1 from ``shedding_options`` and building it,
    and implemented that way so the two cannot disagree. Any of
    ``reference_event``, ``unit`` or ``model`` may be passed to pin that key and
    rank within the remainder.

    Always returns a ``SheddingEnsemble``, single-component when one study
    matched -- ``make_ensemble`` guarantees a one-component ensemble consumes the
    generator exactly as the underlying fit does, so callers keep one code path
    however many studies backed the answer.

    Args:
        biomarker: e.g. ``"SARS-CoV-2"``. Optional, but omitting it will usually
            leave candidates from different biomarkers to rank against one
            another.
        specimen: e.g. ``"stool"``.
        catalog: Catalog to choose from. Defaults to the shipped one.
        weights: Passed to ``make_ensemble``.
        method: Passed to ``make_ensemble``.
        **keys: Further filters, e.g. ``model="exponential"``.

    Returns:
        A ``SheddingEnsemble`` whose ``selection`` records the choice.

    Raises:
        ValueError: If nothing matches.
    """
    from .shedding_ensemble import make_ensemble

    if biomarker is not None:
        keys["biomarker"] = biomarker
    if specimen is not None:
        keys["specimen"] = specimen

    catalog = load_shedding_catalog() if catalog is None else catalog
    options = shedding_options(catalog=catalog, **keys)
    best = options.iloc[0]
    runner_up = options.iloc[1] if len(options) > 1 else None

    groups = _grouped(_matching_fits(catalog, keys))
    # Matched on the sortable form rather than by indexing with a tuple built
    # from the row: pandas may hand back NaN where the fit held None, and NaN
    # does not equal itself, so a direct dict lookup would miss the group for
    # any fit with no unit or no reference event.
    target = tuple(_sortable(best[key]) for key in _GROUP_KEYS)
    components = next(
        group
        for signature, group in groups.items()
        if tuple(_sortable(value) for value in signature) == target
    )

    ensemble = make_ensemble(components, weights=weights, method=method)
    ensemble.selection = Selection(
        picked={key: best[key] for key in list(_GROUP_KEYS) + [
            "event_class", "n_unit_studies", "n_studies", "n_subjects",
            "n_measurements",
        ]},
        passed_over=options.iloc[1:].reset_index(drop=True),
        reason=_reason(best, runner_up),
        analytes={fit.dataset_id: fit.analyte for fit in components},
    )
    return ensemble
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/test_shedding_select.py tests/test_shedding_ensemble.py -q`
Expected: PASS. `test_every_advertised_group_can_actually_be_built` is the slow one; it iterates all 82 groups.

- [ ] **Step 6: Format and commit**

```bash
black shedding_hub/shedding_select.py shedding_hub/shedding_ensemble.py tests/test_shedding_select.py
git add shedding_hub/shedding_select.py shedding_hub/shedding_ensemble.py tests/test_shedding_select.py
git commit -m "feat: pick a defensible shedding source, and say why"
```

---

### Task 4: Stop claiming an infection time origin that was not earned

**Files:**
- Modify: `shedding_hub/shedding_simulate.py:107-150`
- Test: `tests/test_shedding_simulate.py`

**Interfaces:**
- Consumes: `classify_reference_event` from Task 1.
- Produces: no new public names. `simulate_shedding` gains no arguments; `result.attrs["time_origin"]` and `result.attrs["reference_event_class"]` change meaning.

This task is independent of Tasks 2 and 3 and may be done before them.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_shedding_simulate.py`:

```python
def _fit_with_reference_event(make_synthetic_dataset, event):
    import numpy as np

    from shedding_hub.shedding_fit import fit_shedding_model

    dataset = make_synthetic_dataset(
        "exponential", np.array([np.log(0.6), np.log(18.0)]), np.diag([0.04, 0.04]),
        n_subjects=12, seed=5,
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_shedding_simulate.py -q -k "origin or administrative or exposure_event or symptom_onset_shifts"`
Expected: FAIL — `assert 'infection' == 'enrollment_shifted'` on the administrative test, and `KeyError: 'reference_event_class'`.

- [ ] **Step 3: Implement the honest time origin**

In `shedding_hub/shedding_simulate.py`, add to the imports at the top:

```python
from .shedding_select import classify_reference_event
```

Then find the `frame.attrs` assignment near the end of `simulate_shedding`:

```python
    frame.attrs = {
        "time_origin": "infection" if incubation_applied else source.reference_event,
        "incubation_applied": incubation_applied,
        "model": source.model,
        "unit": source.unit,
        "biomarker": getattr(source, "biomarker", None),
        "specimen": getattr(source, "specimen", None),
    }
    return frame
```

Replace with:

```python
    event = source.reference_event
    event_class = classify_reference_event(event)
    time_origin = event
    if incubation_applied:
        if event_class == "landmark":
            time_origin = "infection"
        else:
            time_origin = f"{event}_shifted"
            if event_class == "exposure":
                warnings.warn(
                    f"{event!r} is already the exposure, so there is no "
                    "incubation period to bridge: shifting moves the origin to "
                    f"before the exposure itself. time_origin is recorded as "
                    f"{time_origin!r}, not 'infection'.",
                    UserWarning,
                    stacklevel=2,
                )
            else:
                warnings.warn(
                    f"{event!r} is an administrative reference event, which has "
                    "no fixed offset from infection -- it reflects testing "
                    "behaviour and health-system access. time_origin is recorded "
                    f"as {time_origin!r}, not 'infection'.",
                    UserWarning,
                    stacklevel=2,
                )

    frame.attrs = {
        "time_origin": time_origin,
        "reference_event_class": event_class,
        "incubation_applied": incubation_applied,
        "model": source.model,
        "unit": source.unit,
        "biomarker": getattr(source, "biomarker", None),
        "specimen": getattr(source, "specimen", None),
    }
    return frame
```

If `warnings` is not already imported in this module, add `import warnings` at the top.

- [ ] **Step 4: Update the `incubation_period` docstring**

In the same function's docstring, find the `incubation_period` argument description beginning `Days from infection to the reference event.` and append to it:

```
            Only meaningful for a reference event that is a natural-history
            landmark (symptom onset). Applying it to an administrative event
            (enrollment, confirmation date, hospital admission, treatment) or to
            the exposure itself (inoculation, vaccination) warns, and
            ``attrs["time_origin"]`` records ``"<event>_shifted"`` rather than
            ``"infection"``.
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/test_shedding_simulate.py -q`
Expected: PASS. If a pre-existing test asserts `attrs["time_origin"] == "infection"` on a non-landmark fixture, update that fixture to use `"symptom onset"` rather than loosening the assertion — the whole point is that the old behaviour was wrong.

- [ ] **Step 6: Format and commit**

```bash
black shedding_hub/shedding_simulate.py tests/test_shedding_simulate.py
git add shedding_hub/shedding_simulate.py tests/test_shedding_simulate.py
git commit -m "fix: only claim an infection time origin when the event earns it"
```

---

### Task 5: Exports, README, and the methods note

**Files:**
- Modify: `shedding_hub/__init__.py`
- Modify: `README.md`
- Modify: `docs/modeling-methods.md`
- Test: `tests/test_shedding_select.py`

**Interfaces:**
- Consumes: everything from Tasks 1-4.
- Produces: `sh.shedding_options`, `sh.shedding_for`, `sh.classify_reference_event` at package level.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_shedding_select.py`:

```python
def test_public_exports():
    import shedding_hub as sh

    assert callable(sh.shedding_options)
    assert callable(sh.shedding_for)
    assert callable(sh.classify_reference_event)
    for name in ("shedding_options", "shedding_for", "classify_reference_event"):
        assert name in sh.__all__
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_shedding_select.py -q -k exports`
Expected: FAIL — `AttributeError: module 'shedding_hub' has no attribute 'shedding_options'`

- [ ] **Step 3: Add the exports**

In `shedding_hub/__init__.py`, after the line `from .shedding_ensemble import SheddingEnsemble, make_ensemble`, add:

```python
from .shedding_select import (
    REFERENCE_EVENT_CLASSES,
    Selection,
    classify_reference_event,
    shedding_for,
    shedding_options,
)
```

And add these five names to the `__all__` list:

```python
    "REFERENCE_EVENT_CLASSES",
    "Selection",
    "classify_reference_event",
    "shedding_for",
    "shedding_options",
```

- [ ] **Step 4: Run the test**

Run: `python -m pytest tests/test_shedding_select.py -q -k exports`
Expected: PASS

- [ ] **Step 5: Add the README section**

In `README.md`, inside the `### Simulating Shedding` section, immediately before the sentence beginning `Three models are available.`, insert:

````markdown
Picking a fit by hand means naming five keys — biomarker, specimen, reference
event, unit and model — and those keys cut the catalog into 82 groups, 71 of
which hold a single study. To see the choice, and to have it made for you:

```python
>>> import shedding_hub as sh
>>> options = sh.shedding_options(biomarker='SARS-CoV-2', specimen='stool')
>>> list(options.columns)
['biomarker', 'specimen', 'reference_event', 'event_class', 'unit', 'n_unit_studies', 'model', 'n_studies', 'n_subjects', 'n_measurements', 'rank']
>>> source = sh.shedding_for('SARS-CoV-2', 'stool')
>>> source.selection.picked['event_class']
'landmark'

```

`shedding_for` takes rank 1 from `shedding_options`, preferring a reference event
that can be placed on an infection timeline, then the unit most studies report
for that biomarker and specimen, then a model that resolves the rise, then the
weight of evidence. Pass `model=`, `unit=` or `reference_event=` to pin any of
them, and read `source.selection` for what was chosen and what it beat.
````

- [ ] **Step 6: Verify the README doctests pass**

Run: `python -m doctest -o ELLIPSIS -o NORMALIZE_WHITESPACE README.md`
Expected: no output (success). These run against the shipped catalog and are deliberately not skipped, so if a future rebuild changes the top-ranked group, this fails loudly — fix it by updating the example, not by adding `+SKIP`.

- [ ] **Step 7: Note the taxonomy in the methods doc**

In `docs/modeling-methods.md`, in section `## 5. Simulation`, append:

```markdown
### Reference events are not interchangeable

The catalog spans seven reference events, in three classes. `inoculation` and
`vaccination` *are* the exposure, so nothing separates them from time zero.
`symptom onset` is a natural-history landmark, a defined and documented offset
from infection. `enrollment`, `confirmation date`, `hospital admission` and
`treatment` are administrative: they record when a subject entered a study, was
tested, or was admitted, which depends on testing behaviour and health-system
access rather than on their infection.

Only the landmark class earns an infection time origin. `simulate_shedding`
warns when an incubation period is applied to either of the others, and records
`time_origin` as `"<event>_shifted"` rather than `"infection"`. `shedding_for`
prefers the classes that can be anchored, which is why it will pass over a
better-supported fit measured from an administrative date.
```

- [ ] **Step 8: Run the whole suite**

Run: `python -m pytest -q` and `black --check .`
Expected: all pass. Confirm `tests/test_parameter_export.py` and
`test_shipped_catalog_covers_every_dataset` are still green — nothing in this
plan regenerates the catalog or the parameter table.

- [ ] **Step 9: Commit**

```bash
black .
git add shedding_hub/__init__.py README.md docs/modeling-methods.md tests/test_shedding_select.py
git commit -m "feat: export the selection API, and document what it prefers"
```

---

## Self-Review Notes

**Spec coverage.** Taxonomy → Task 1. `shedding_options` → Task 2. `shedding_for`,
`Selection`, single-component equivalence → Task 3. One-analyte-per-study
amendment → Task 2 (`_reduce_to_one_fit_per_study`) with its regression guard in
Task 3. Time-origin fix → Task 4. Exports and docs → Task 5. Every testing bullet
in the spec maps to a named test.

**Ordering.** Task 4 depends only on Task 1, so it can be done at any point after
it. Tasks 2 and 3 are strictly sequential.

**Naming consistency.** `classify_reference_event`, `shedding_options`,
`shedding_for`, `Selection`, `_reduce_to_one_fit_per_study`, `_rank_key`,
`_grouped`, `_matching_fits`, `_GROUP_KEYS` are used identically wherever they
appear across tasks.

**Known risk.** Task 3 sets `ensemble.selection` after construction. The field is
declared on the dataclass in Task 3 Step 3 so this is an ordinary assignment, not
a dynamic attribute — but `SheddingEnsemble.to_dict` deliberately does not
serialize it, and `from_dict` therefore returns `selection=None`. That asymmetry
is intended: a selection describes a choice made against one catalog and would
mislead if restored beside components from another.
