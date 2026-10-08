import csv
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

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


# The curation step: non-detect coding and log10, by each extraction script's rule.

REPO_ROOT = INPUTS.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import manual_curation as mc  # noqa: E402

# Studies whose script goes on to do something the curation step leaves alone,
# so their values are not expected to equal the reference's. NOTES.csv says what.
NOT_REPRODUCED = {
    "hakki2022onset",  # the script drops some rows
    "kimse2020viral",  # digitized cycles rounded
    "lavezzo2020suppression",  # dates worked out from the daily swab columns
}


def _manual_rows():
    with (REPO_ROOT / "benchmark" / "manifest.csv").open(encoding="utf-8") as fh:
        return [row for row in csv.DictReader(fh) if row["era"] == "manual"]


def _values(study):
    with (INPUTS / f"{study}.csv").open(encoding="utf-8", newline="") as fh:
        return [row["value"] for row in csv.DictReader(fh)]


def _number(text):
    try:
        return float(text)
    except ValueError:
        return None


def test_every_manual_study_has_a_curation_rule():
    assert set(mc.RULES) == {row["study_id"] for row in _manual_rows()}


def test_curated_values_are_numbers_or_schema_words():
    for row in _manual_rows():
        words = {v for v in _values(row["study_id"]) if _number(v) is None}
        assert words <= {"negative", "positive", "inconclusive"}, row["study_id"]


def test_curated_values_reproduce_the_reference_where_the_script_stops_there():
    """Counts of non-detects, and every number, as in the reference dataset.

    This reads released data to check the builders. It never feeds a CSV.
    """
    for row in _manual_rows():
        study = row["study_id"]
        if study in NOT_REPRODUCED:
            continue
        reference = yaml.safe_load(
            (REPO_ROOT / row["reference_yaml"]).read_text(encoding="utf-8")
        )
        expected = [
            m["value"] for p in reference["participants"] for m in p["measurements"]
        ]
        values = _values(study)
        assert len(values) == len(expected), study
        for word in ("negative", "positive", "inconclusive"):
            assert values.count(word) == expected.count(word), (study, word)
        numbers = sorted(v for v in map(_number, values) if v is not None)
        wanted = sorted(float(v) for v in expected if not isinstance(v, str))
        assert numbers == pytest.approx(wanted, rel=1e-9), study


def test_combined_dataset_limit_of_detection_names_its_unit():
    for row in _manual_rows():
        with (INPUTS / f"{row['study_id']}.csv").open(encoding="utf-8") as fh:
            header = fh.readline().strip().split(",")
        assert "LOD" not in header, row["study_id"]


def _table(values, **columns):
    return pd.DataFrame({"value": values, **columns})


def test_curation_rules_on_small_tables():
    table, _ = mc.curate("tan2021early", _table(["35.41", "38", "38.0"]))
    assert list(table["value"]) == ["35.41", "negative", "negative"]

    table, _ = mc.curate("peiris2003clinical", _table(["3", "4.5"]))
    assert list(table["value"]) == ["1000.0", repr(10**4.5)]

    table, _ = mc.curate("woelfel2020virological", _table(["0.01", "2"]))
    assert list(table["value"]) == ["negative", "100.0"]

    table, _ = mc.curate(
        "covid2020clinical", _table(["30.2", "44.9", "55", "negative"])
    )
    assert list(table["value"]) == ["30.2", "inconclusive", "negative", "negative"]

    table, _ = mc.curate(
        "gautret2020hydroxychloroquine",
        _table(["35", "22"], CombinedDataset_value=["1", "4718.5"]),
    )
    assert list(table["value"]) == ["negative", "22"]

    table, note = mc.curate(
        "arts2023longitudinal",
        _table(
            ["5", "7", "9", "2"],
            analyte=[
                "N_conc (gc/mg-dw)",
                "N_conc (gc/mg-dw)",
                "crAss_conc",
                "crAss_conc",
            ],
            det=["True", "False", "True", "False"],
            quant=[None, None, "False", "False"],
        ),
    )
    assert list(table["value"]) == ["5", "negative", "positive", "negative"]
    assert "2 cells written as negative" in note
