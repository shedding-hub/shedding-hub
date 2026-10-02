import csv
from pathlib import Path

INPUTS = Path(__file__).resolve().parent.parent / "benchmark" / "manual_inputs"


def _notes():
    with (INPUTS / "NOTES.csv").open(encoding="utf-8") as fh:
        return {row["study_id"]: row for row in csv.DictReader(fh)}


def test_every_input_csv_has_a_note_and_the_standard_columns():
    notes = _notes()
    files = sorted(p for p in INPUTS.glob("*.csv") if p.name != "NOTES.csv")
    assert {p.stem for p in files} == set(notes)
    for path in files:
        with path.open(encoding="utf-8") as fh:
            rows = list(csv.reader(fh))
        assert rows[0][:3] == ["PatientID", "time", "value"], path.name
        note = notes[path.stem]
        # One row per measurement: the note records how many were written.
        assert len(rows) - 1 == int(note["rows"]), path.name
        assert ";".join(rows[0]) == note["columns"], path.name
        assert note["note"].strip(), path.name


def test_notes_have_no_em_dash():
    assert all("—" not in row["note"] for row in _notes().values())
