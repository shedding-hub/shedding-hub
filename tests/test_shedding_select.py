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
