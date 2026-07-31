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
