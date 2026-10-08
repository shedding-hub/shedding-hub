"""The comparison code recovers injected errors and reports nothing else.

Each test builds a synthetic draft by perturbing a copy of a released dataset.
Every draft also has its participants reordered and its analyte keys renamed,
since neither carries meaning.
"""

import copy
import random
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import align  # noqa: E402
import discrepancies as disc  # noqa: E402
import match_findings as mf  # noqa: E402

EXCLUDED = {"jones2021estimating"}  # population-scale; not in the analysis set
STUDIES = sorted(
    p.parent.name
    for p in (REPO_ROOT / "data").glob("*/*.yaml")
    if p.stem not in EXCLUDED
)
PERTURBED = [
    "woelfel2020virological",
    "tan2021early",
    "iwakiri2009quantitative",
    "antinori2022epidemiological",
]
BOTH = [align.REPORTED, align.FIGURE]


def load(study):
    path = REPO_ROOT / "data" / study / f"{study}.yaml"
    return yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader)


def disguise(dataset, seed=0):
    """A copy with participants shuffled and analyte keys renamed.

    Returns the copy and, for each reference participant, its new position.
    """
    draft = copy.deepcopy(dataset)
    rng = random.Random(seed)
    names = {key: f"renamed_{n}" for n, key in enumerate(draft["analytes"])}
    draft["analytes"] = {names[k]: v for k, v in draft["analytes"].items()}
    order = list(range(len(draft["participants"])))
    rng.shuffle(order)
    for participant in draft["participants"]:
        for m in participant["measurements"]:
            m["analyte"] = names[m["analyte"]]
    draft["participants"] = [draft["participants"][i] for i in order]
    return draft, {old: new for new, old in enumerate(order)}


def unique_readings(dataset, numeric=None, minimum=4):
    """(participant, measurement) positions that cannot be confused.

    A reading qualifies when no other measurement of its participant has the
    same analyte and either the same time or the same value.
    """
    out = []
    for i, participant in enumerate(dataset["participants"]):
        ms = participant["measurements"]
        if len(ms) < minimum:
            continue
        for j, m in enumerate(ms):
            if numeric is not None and align.is_number(m["value"]) != numeric:
                continue
            if not align.is_number(m["time"]):
                continue
            clash = any(
                k != j
                and o["analyte"] == m["analyte"]
                and (o["time"] == m["time"] or o["value"] == m["value"])
                for k, o in enumerate(ms)
            )
            if not clash:
                out.append((i, j))
    return out


def found(study, draft, reference, tol):
    result = disc.compare(study, draft, reference, tol)
    return sorted((r["type"], r["s2_path"]) for r in result.rows), result


def spread(items, n=3):
    """A few items from different participants, chosen the same way each time."""
    out, seen = [], set()
    for i, j in items:
        if i not in seen:
            seen.add(i)
            out.append((i, j))
    assert len(out) >= 1
    return out[:n]


@pytest.mark.parametrize("study", STUDIES)
def test_a_disguised_copy_has_no_discrepancies(study):
    reference = load(study)
    draft, _ = disguise(reference)
    for tol in BOTH:
        rows, result = found(study, draft, reference, tol)
        assert rows == []
        counts = result.counts
        assert counts["correct_measurements"] == counts["s2_measurements"]
        assert counts["matched_participants"] == counts["s2_participants"]


@pytest.mark.parametrize("study", PERTURBED)
@pytest.mark.parametrize("tol", BOTH)
def test_dropped_measurements_are_omissions(study, tol):
    reference = load(study)
    draft, where = disguise(reference)
    picks = spread(unique_readings(reference))
    for i, j in sorted(picks, reverse=True):
        del draft["participants"][where[i]]["measurements"][j]
    rows, _ = found(study, draft, reference, tol)
    assert rows == sorted(
        ("omission", f"/participants/{i}/measurements/{j}") for i, j in picks
    )


@pytest.mark.parametrize("study", PERTURBED)
@pytest.mark.parametrize("tol", BOTH)
def test_shifted_values_are_wrong_values(study, tol):
    reference = load(study)
    draft, where = disguise(reference)
    picks = spread(unique_readings(reference, numeric=True))
    for i, j in picks:
        draft["participants"][where[i]]["measurements"][j]["value"] *= 31.7
    rows, _ = found(study, draft, reference, tol)
    assert rows == sorted(
        ("wrong value", f"/participants/{i}/measurements/{j}/value") for i, j in picks
    )


@pytest.mark.parametrize("study", PERTURBED)
@pytest.mark.parametrize("tol", BOTH)
def test_shifted_times_are_wrong_times(study, tol):
    reference = load(study)
    draft, where = disguise(reference)
    picks = spread(unique_readings(reference, numeric=True))
    for i, j in picks:
        draft["participants"][where[i]]["measurements"][j]["time"] += 1000.25
    rows, _ = found(study, draft, reference, tol)
    assert rows == sorted(
        ("wrong time", f"/participants/{i}/measurements/{j}/time") for i, j in picks
    )


@pytest.mark.parametrize("study", PERTURBED)
def test_small_shifts_pass_only_the_figure_tolerances(study):
    reference = load(study)
    draft, where = disguise(reference)
    (i, j), (k, l) = unique_readings(reference, numeric=True)[:2]
    unit = reference["analytes"][
        reference["participants"][i]["measurements"][j]["analyte"]
    ]["unit"]
    nudged = draft["participants"][where[i]]["measurements"][j]
    nudged["value"] = (
        nudged["value"] + 0.4 if unit == align.CT_UNIT else nudged["value"] * 1.2
    )
    draft["participants"][where[k]]["measurements"][l]["time"] += 0.4
    assert found(study, draft, reference, align.FIGURE)[0] == []
    assert found(study, draft, reference, align.REPORTED)[0] == sorted(
        [
            ("wrong value", f"/participants/{i}/measurements/{j}/value"),
            ("wrong time", f"/participants/{k}/measurements/{l}/time"),
        ]
    )


@pytest.mark.parametrize("study", PERTURBED)
@pytest.mark.parametrize("tol", BOTH)
def test_flipped_non_detects_are_reported(study, tol):
    reference = load(study)
    draft, where = disguise(reference)
    i, j = unique_readings(reference, numeric=True)[0]
    draft["participants"][where[i]]["measurements"][j]["value"] = "negative"
    expected = [
        ("wrong non-detect status", f"/participants/{i}/measurements/{j}/value")
    ]
    negatives = [
        (p, q)
        for p, participant in enumerate(reference["participants"])
        for q, m in enumerate(participant["measurements"])
        if m["value"] == "negative" and p != i
    ]
    if negatives:
        p, q = negatives[0]
        draft["participants"][where[p]]["measurements"][q]["value"] = "positive"
        expected.append(
            ("wrong non-detect status", f"/participants/{p}/measurements/{q}/value")
        )
    assert found(study, draft, reference, tol)[0] == sorted(expected)


@pytest.mark.parametrize("study", PERTURBED)
@pytest.mark.parametrize("tol", BOTH)
def test_swapped_analyte_fields_are_field_errors(study, tol):
    reference = load(study)
    draft, _ = disguise(reference)
    keys = list(reference["analytes"])
    first, last = keys[0], keys[-1]
    old_unit = reference["analytes"][first]["unit"]
    old_event = reference["analytes"][last]["reference_event"]
    draft["analytes"]["renamed_0"]["unit"] = (
        "pfu/mL" if old_unit != "pfu/mL" else "gc/mL"
    )
    draft["analytes"][f"renamed_{len(keys) - 1}"]["reference_event"] = (
        "enrollment" if old_event != "enrollment" else "exposure"
    )
    assert found(study, draft, reference, tol)[0] == sorted(
        [
            ("wrong unit", f"/analytes/{first}/unit"),
            ("wrong reference_event", f"/analytes/{last}/reference_event"),
        ]
    )


@pytest.mark.parametrize("study", PERTURBED)
def test_a_dropped_participant_is_missing(study):
    reference = load(study)
    draft, where = disguise(reference)
    del draft["participants"][where[2]]
    rows, result = found(study, draft, reference, align.REPORTED)
    assert rows == [
        ("missing", "/participants/2"),
        ("participant-count mismatch", "/participants"),
    ]
    missing = next(r for r in result.rows if r["type"] == "missing")
    assert missing["n_measurements"] == len(
        reference["participants"][2]["measurements"]
    )
    assert missing["s0_path"] == "" and missing["s0_anchor"] == "/participants"


def test_an_invented_participant_is_extra():
    reference = load("woelfel2020virological")
    draft, _ = disguise(reference)
    key = next(iter(draft["analytes"]))
    draft["participants"].append(
        {
            "measurements": [
                {"analyte": key, "time": 500 + n, "value": 3.3e4 * (n + 1)}
                for n in range(3)
            ]
        }
    )
    rows, result = found("x", draft, reference, align.FIGURE)
    assert [r[0] for r in rows] == ["extra", "participant-count mismatch"]
    extra = next(r for r in result.rows if r["type"] == "extra")
    assert extra["level"] == "participant" and extra["n_measurements"] == 3
    assert extra["s0_path"] == f"/participants/{len(reference['participants'])}"


def test_missing_and_extra_analytes():
    reference = load("woelfel2020virological")
    draft, _ = disguise(reference)
    gone = "renamed_1"
    del draft["analytes"][gone]
    for participant in draft["participants"]:
        participant["measurements"] = [
            m for m in participant["measurements"] if m["analyte"] != gone
        ]
    draft["analytes"]["urine"] = {
        **draft["analytes"]["renamed_0"],
        "specimen": "urine",
        "biomarker": "mpox",
    }
    rows, result = found("x", draft, reference, align.REPORTED)
    reference_key = list(reference["analytes"])[1]
    assert rows == [
        ("extra analyte", ""),
        ("missing analyte", f"/analytes/{reference_key}"),
    ]
    missing = next(r for r in result.rows if r["type"] == "missing analyte")
    assert missing["n_measurements"] == sum(
        m["analyte"] == reference_key
        for p in reference["participants"]
        for m in p["measurements"]
    )


def test_a_wrong_specimen_is_one_field_error():
    reference = load("tan2021early")
    draft, _ = disguise(reference)
    draft["analytes"]["renamed_0"]["specimen"] = "stool"
    key = next(iter(reference["analytes"]))
    assert found("x", draft, reference, align.REPORTED)[0] == [
        ("wrong specimen", f"/analytes/{key}/specimen")
    ]


def test_attribute_errors():
    reference = load("tan2021early")
    draft, where = disguise(reference)
    draft["participants"][where[0]]["attributes"]["age"] += 1
    moved = draft["participants"][where[1]]["attributes"]
    moved["age_group"] = moved.pop("sex")
    del draft["participants"][where[2]]["attributes"]["sex"]
    draft["participants"][where[3]]["attributes"]["hospitalized"] = True
    rows, result = found("x", draft, reference, align.REPORTED)
    assert rows == sorted(
        [
            ("wrong attribute value", "/participants/0/attributes/age"),
            ("attribute in the wrong field", "/participants/1/attributes/sex"),
            ("missing attribute", "/participants/2/attributes/sex"),
            ("extra attribute", ""),
        ]
    )
    extra = next(r for r in result.rows if r["type"] == "extra attribute")
    assert extra["s0_path"] == f"/participants/{where[3]}/attributes/hospitalized"


def test_a_draft_with_every_value_wrong_keeps_its_participants():
    reference = load("woelfel2020virological")
    draft, _ = disguise(reference)
    for participant in draft["participants"]:
        for m in participant["measurements"]:
            if align.is_number(m["value"]):
                m["value"] *= 1000
    rows, result = found("x", draft, reference, align.FIGURE)
    assert {kind for kind, _ in rows} == {"wrong value"}
    assert result.counts["matched_participants"] == len(reference["participants"])


@pytest.mark.parametrize(
    "draft", [{}, {"participants": "none", "analytes": []}, {"participants": [1, None]}]
)
def test_malformed_drafts_are_compared_not_crashed_on(draft):
    reference = load("tan2021early")
    result = disc.compare("x", draft, reference, align.REPORTED)
    kinds = {r["type"] for r in result.rows}
    assert "missing analyte" in kinds and result.counts["correct_measurements"] == 0


def test_no_draft_is_one_row_counting_every_reference_measurement():
    reference = load("tan2021early")
    result = disc.no_draft("tan2021early", reference, "no draft")
    assert [r["type"] for r in result.rows] == ["no draft"]
    assert result.rows[0]["n_measurements"] == result.counts["s2_measurements"] == 82


# --------------------------------------------------------------------- findings


def _rows_with_errors():
    reference = load("tan2021early")
    draft, where = disguise(reference)
    i, j = unique_readings(reference, numeric=True)[0]
    draft["participants"][where[i]]["measurements"][j]["value"] += 9
    del draft["participants"][where[i]]["measurements"][j + 1]
    draft["participants"][where[i]]["attributes"]["age"] += 1
    draft["analytes"]["renamed_0"]["unit"] = "gc/mL"
    del draft["participants"][where[5]]
    result = disc.compare("tan2021early", draft, reference, align.REPORTED)
    rows = {r["type"]: r for r in result.rows}
    return rows, int(rows["wrong value"]["s0_path"].split("/")[2])


def _finding(check, path, severity="warning"):
    return {
        "id": "f",
        "check": check,
        "severity": severity,
        "path": path,
        "message": "m",
    }


def test_findings_hit_by_location_and_type():
    rows, person = _rows_with_errors()
    assert set(rows) == {
        "wrong value",
        "omission",
        "wrong attribute value",
        "wrong unit",
        "missing",
        "participant-count mismatch",
    }
    cases = [
        ("analyte.value_vs_paper", "/analytes/renamed_0/unit", ["wrong unit"]),
        ("analyte.value_vs_paper", "/analytes/renamed_0", ["wrong unit"]),
        ("analyte.value_vs_paper", "/analytes/renamed_0/specimen", []),
        (
            "attr.value_vs_paper",
            f"/participants/{person}/attributes/age",
            ["wrong attribute value"],
        ),
        ("attr.value_vs_paper", f"/participants/{person}", ["wrong attribute value"]),
        ("sample.count_vs_paper", f"/participants/{person}/measurements", ["omission"]),
        ("sample.count_vs_paper", "/", ["omission"]),
        (
            "participant.count_vs_paper",
            "/participants",
            ["missing", "participant-count mismatch"],
        ),
        ("participant.count_vs_paper", "", ["missing", "participant-count mismatch"]),
        ("value.censored_string", rows["wrong value"]["s0_path"], ["wrong value"]),
        # The right place but a check that cannot speak about the error there.
        ("analyte.value_vs_paper", rows["wrong value"]["s0_path"], ["wrong value"]),
        ("attr.value_vs_paper", rows["wrong value"]["s0_path"], []),
        # A path that merely starts with the same characters is not an ancestor.
        ("participant.count_vs_paper", "/participant", []),
    ]
    out = mf.match(
        "tan2021early", [_finding(c, p) for c, p, _ in cases], list(rows.values())
    )
    for row, (check, path, expected) in zip(out, cases):
        ids = sorted(rows[kind]["id"] for kind in expected)
        assert sorted(filter(None, row["matched_discrepancies"].split(";"))) == ids, (
            check,
            path,
        )
        assert row["hit"] == ("yes" if expected else "no")
        assert row["category"] == "data"
    wrong_check = out[11]
    assert wrong_check["same_location_other_type"] == rows["wrong value"]["id"]
    assert [r["adjudicate"] for r in out] == [
        "no" if expected else "yes" for _, _, expected in cases
    ]


def test_gene_target_wording_and_extra_attributes_are_listed_but_not_scored():
    reference = load("tan2021early")
    draft, where = disguise(reference)
    draft["analytes"]["renamed_0"]["gene_target"] = "N gene (CDC N1 assay)"
    draft["participants"][where[0]]["attributes"]["hospitalized"] = True
    draft["participants"][where[1]]["attributes"]["age"] += 1
    result = disc.compare("x", draft, reference, align.REPORTED)
    assert {r["type"]: r["scored"] for r in result.rows} == {
        "wrong gene_target": "no",
        "extra attribute": "no",
        "wrong attribute value": "yes",
    }
    # A finding about a difference that is not scored has not found an error.
    extra = next(r for r in result.rows if r["type"] == "extra attribute")
    out = mf.match(
        "x", [_finding("attr.value_vs_paper", extra["s0_path"])], result.rows
    )
    assert out[0]["hit"] == "no" and out[0]["adjudicate"] == "yes"
    assert out[0]["same_location_other_type"] == extra["id"]


def test_schema_gaps_notes_and_unknown_checks_are_not_scored():
    rows, _ = _rows_with_errors()
    findings = [
        _finding("schema.enum_gap_vs_paper", "/analytes/renamed_0/unit"),
        _finding("model.observation", "/", "info"),
        _finding("analyte.unused", "/analytes/renamed_0"),
    ]
    out = mf.match("tan2021early", findings, list(rows.values()))
    assert [r["category"] for r in out] == ["schema_gap", "observation", "unmapped"]
    assert {r["hit"] for r in out} == {"no"} == {r["adjudicate"] for r in out}
    assert out[0]["same_location_other_type"] == rows["wrong unit"]["id"]


def test_review_files_are_checked_for_form(tmp_path):
    good = tmp_path / "good.yaml"
    good.write_text(
        yaml.safe_dump(
            {"dataset": "x", "findings": [_finding("model.observation", "/", "info")]}
        ),
        encoding="utf-8",
    )
    assert mf.read_findings(good) == ([_finding("model.observation", "/", "info")], [])
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        yaml.safe_dump(
            {"findings": [{"check": "c", "severity": "fatal", "path": "/"}, 3]}
        ),
        encoding="utf-8",
    )
    findings, problems = mf.read_findings(bad)
    assert len(findings) == 1 and len(problems) == 3
    assert mf.resolves({"a": [{"b": 1}]}, "/a/0/b") and not mf.resolves(
        {"a": []}, "/a/0"
    )
