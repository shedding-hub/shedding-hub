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

from dataclasses import dataclass, field

import pandas as pd

from .shedding_catalog import SheddingCatalog, load_shedding_catalog

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
# Rule 2 -- the unit, by how many studies report it -- sits between these two and
# is computed per candidate set rather than being a constant, so it lives on the
# row as n_unit_studies. Units are incommensurable, and without an explicit rule
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

    # Rule 2's input: how many distinct studies report each unit anywhere in the
    # candidate set. A property of the unit, not of any one group, so it is
    # counted once here and carried on every row that shares the unit.
    studies_by_unit = {}
    for fit in matches:
        studies_by_unit.setdefault(_sortable(fit.unit), set()).add(fit.dataset_id)

    rows = []
    for signature, group in _grouped(matches).items():
        row = dict(zip(_GROUP_KEYS, signature))
        row["event_class"] = classify_reference_event(row["reference_event"])
        row["n_unit_studies"] = len(studies_by_unit[_sortable(row["unit"])])
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
