"""
Choose which fitted estimate to simulate from.

The catalog's six compatibility keys -- biomarker, specimen, reference event,
unit, value type and model -- cut 126 fits into 82 groups, 71 of them a single
study. So the user's problem is not combining many estimates but choosing among
incommensurable ones: gc/mL and gc/dry gram are different quantities, symptom
onset and enrollment are different clocks, a Ct fit's height and a
concentration fit's are different scales, and the three models are fitted to
different observation sets and so are not AIC-comparable.

This module makes that choice visible (``shedding_options``) and makes a
documented, overridable default (``shedding_for``). It adds no statistical
assumption: nothing here pools across units or reference events.

Examples:
    >>> import shedding_hub as sh
    >>> sorted(sh.REFERENCE_EVENT_CLASSES)[:3]
    ['confirmation date', 'enrollment', 'hospital admission']
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


def classify_reference_event(event: str | None) -> str:
    """
    Classify a reference event as ``exposure``, ``landmark`` or ``administrative``.

    Unrecognised events -- including ``None`` -- are ``administrative``. That is
    the conservative direction: a reference event the package has never seen has
    not demonstrated a defined offset from infection, so it should not inherit
    one.

    Examples:
        >>> import shedding_hub as sh
        >>> sh.classify_reference_event('enrollment')
        'administrative'
    """
    return REFERENCE_EVENT_CLASSES.get(event, "administrative")


_GROUP_KEYS = (
    "biomarker",
    "specimen",
    "reference_event",
    "unit",
    "value_type",
    "model",
)

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
    """Group fits by the six compatibility keys, reduced to one per study."""
    groups = {}
    for fit in fits:
        signature = tuple(getattr(fit, key, None) for key in _GROUP_KEYS)
        groups.setdefault(signature, []).append(fit)
    return {
        signature: _reduce_to_one_fit_per_study(group)
        for signature, group in groups.items()
    }


def shedding_options(
    biomarker=None,
    specimen=None,
    *,
    catalog: SheddingCatalog | None = None,
    **keys,
) -> pd.DataFrame:
    """
    List every group that could be simulated from, best first.

    Each row is one combination of the six keys ``make_ensemble`` requires
    agreement on, so each row can actually be built -- counts are reported after
    reducing each study to a single analyte. ``rank`` 1 is what ``shedding_for``
    returns for the same arguments.

    Args:
        biomarker (str | None): e.g. ``"SARS-CoV-2"``. Positional, matching
            ``shedding_for``.
        specimen (str | None): e.g. ``"stool"``.
        catalog: Catalog to search. Defaults to the shipped one.
        **keys (Any): Further attribute filters, e.g. ``model="gamma"``.

    Returns:
        A ``DataFrame`` with one row per group.

    Raises:
        ValueError: If nothing matches ``keys``.

    Examples:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> options = sh.shedding_options('SARS-CoV-2', 'stool', catalog=catalog)
        >>> options.shape
        (10, 12)
        >>> options.iloc[0]['model']
        'gamma'
    """
    if biomarker is not None:
        keys["biomarker"] = biomarker
    if specimen is not None:
        keys["specimen"] = specimen

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
            "value_type",
            "n_unit_studies",
            "model",
            "n_studies",
            "n_subjects",
            "n_measurements",
            "rank",
        ]
    ]


@dataclass
class Selection:
    """
    Why ``shedding_for`` returned what it did.

    Examples:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> selection = sh.shedding_for(
        ...     'SARS-CoV-2', 'stool', catalog=catalog
        ... ).selection
        >>> selection.reason
        'its model resolves the rise'
        >>> str(selection)
        'picked symptom onset / gc/mL / gamma (2 study/studies, 16 subjects); its model resolves the rise'
    """

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
        biomarker (str | None): e.g. ``"SARS-CoV-2"``. Optional, but omitting
            it will usually leave candidates from different biomarkers to
            rank against one another.
        specimen (str | None): e.g. ``"stool"``.
        catalog: Catalog to choose from. Defaults to the shipped one.
        weights (str | Sequence[float]): Passed to ``make_ensemble``. An
            explicit array is applied in component order, which is by
            ``dataset_id`` -- see ``ensemble.components``.
        method: Passed to ``make_ensemble``.
        **keys (Any): Further filters, e.g. ``model="exponential"``.

    Returns:
        source (SheddingEnsemble): Whose ``selection`` records the choice.

    Raises:
        ValueError: If nothing matches.

    Examples:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> source = sh.shedding_for('SARS-CoV-2', 'stool', catalog=catalog)
        >>> source.model
        'gamma'
        >>> source.selection.reason
        'its model resolves the rise'
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
    # Matched on the sortable form for the same reason _sortable exists at all:
    # None and str do not compare, and both appear in these keys, so a group
    # signature holding None must be compared in the normalised form rather than
    # raw.
    target = tuple(_sortable(best[key]) for key in _GROUP_KEYS)
    components = next(
        group
        for signature, group in groups.items()
        if tuple(_sortable(value) for value in signature) == target
    )

    ensemble = make_ensemble(components, weights=weights, method=method)
    ensemble.selection = Selection(
        picked={
            key: best[key]
            for key in list(_GROUP_KEYS)
            + [
                "event_class",
                "n_unit_studies",
                "n_studies",
                "n_subjects",
                "n_measurements",
            ]
        },
        passed_over=options.iloc[1:].reset_index(drop=True),
        reason=_reason(best, runner_up),
        analytes={fit.dataset_id: fit.analyte for fit in components},
    )
    return ensemble
