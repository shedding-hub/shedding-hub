"""Runs are compared with one another in the reference's numbering."""

import copy
import csv
import io
import math
import random
import statistics
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import adjudication as adj  # noqa: E402
import align  # noqa: E402
import discrepancies as disc  # noqa: E402
import match_findings as mf  # noqa: E402
import reproducibility as rep  # noqa: E402

STUDY = "tan2021early"


def load(study):
    path = REPO_ROOT / "data" / study / f"{study}.yaml"
    return yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader)


def as_csv(rows, columns):
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    buffer.seek(0)
    return list(csv.DictReader(buffer))


def shuffled(dataset, seed):
    """A copy with participants reordered; returns it and the new positions."""
    draft = copy.deepcopy(dataset)
    order = list(range(len(draft["participants"])))
    random.Random(seed).shuffle(order)
    draft["participants"] = [draft["participants"][i] for i in order]
    return draft, {old: new for new, old in enumerate(order)}


def run(draft, reference, findings=()):
    """A run as compare_runs takes it, with findings given in draft paths."""
    result = disc.compare(STUDY, draft, reference, align.REPORTED)
    rows = as_csv(result.rows, disc.COLUMNS)
    return {
        "study": adj.Study(draft, reference, align.REPORTED),
        "summary": as_csv([result.counts], disc.SUMMARY_COLUMNS)[0],
        "pairs": as_csv(result.pairs, disc.PAIR_COLUMNS),
        "findings": as_csv(mf.match(STUDY, list(findings), rows), mf.COLUMNS),
    }


def finding(id, path):
    return {
        "id": id,
        "check": "analyte.value_vs_paper",
        "severity": "warning",
        "path": path,
        "message": "m",
    }


def three_runs(change_unit=False):
    reference = load(STUDY)
    n = sum(len(p["measurements"]) for p in reference["participants"])
    # The reading the second and third runs get wrong: numeric, with a time
    # and value no other reading of its participant has.
    i, j = next(
        (i, j)
        for i, p in enumerate(reference["participants"])
        for j, m in enumerate(p["measurements"])
        if align.is_number(m["value"])
        and not any(
            k != j and (o["time"] == m["time"] or o["value"] == m["value"])
            for k, o in enumerate(p["measurements"])
        )
    )
    first, _ = shuffled(reference, 1)
    second, where2 = shuffled(reference, 2)
    third, where3 = shuffled(reference, 3)
    for draft, where in ((second, where2), (third, where3)):
        draft["participants"][where[i]]["measurements"][j]["value"] += 9
    if change_unit:
        key = next(iter(third["analytes"]))
        third["analytes"][key]["unit"] = "gc/mL"
    runs = [
        run(first, reference),
        run(
            second,
            reference,
            [finding("a", f"/participants/{where2[i]}/measurements/{j}/value")],
        ),
        run(
            third,
            reference,
            [finding("b", f"/participants/{where3[i]}/measurements/{j}/value")],
        ),
    ]
    return runs, n, reference


def test_runs_that_agree_after_renumbering_are_the_same():
    runs, n, reference = three_runs()
    out = rep.compare_runs(runs)
    fields = len(reference["analytes"]) * 6
    assert (out["fields_compared"], out["fields_same"]) == (fields, fields)
    # Runs 2 and 3 miss the same reading; run 1 misses none.
    expected = (2 * (n - 1) / n + 1) / 3
    assert math.isclose(out["reproduced_jaccard"], expected)
    recalls = [1.0, (n - 1) / n, (n - 1) / n]
    assert out["recall_by_run"] == ";".join(f"{x:.4f}" for x in recalls)
    assert math.isclose(out["recall_sd"], statistics.stdev(recalls))
    # The two findings sit at different draft paths but the same reading.
    assert math.isclose(out["finding_jaccard"], 1 / 3)


def test_a_field_one_run_changes_is_not_the_same():
    runs, _, reference = three_runs(change_unit=True)
    out = rep.compare_runs(runs)
    assert out["fields_compared"] - out["fields_same"] == 1


def test_a_run_without_a_draft_shares_nothing():
    runs, _, reference = three_runs()
    runs[2] = {
        "study": None,
        "summary": {
            "s2_measurements": runs[0]["summary"]["s2_measurements"],
            "correct_measurements": "0",
        },
        "pairs": [],
        "findings": [],
    }
    out = rep.compare_runs(runs)
    assert out["fields_compared"] == 0
    assert out["recall_by_run"].endswith(";0.0000")


def test_places_are_translated_to_the_reference():
    reference = load(STUDY)
    draft, where = shuffled(reference, 7)
    r = run(draft, reference)
    key = next(iter(draft["analytes"]))
    assert (
        rep.place(
            f"/participants/{where[4]}/measurements/1/time", r["study"], r["pairs"]
        )
        == "/participants/4/measurements/1/time"
    )
    assert (
        rep.place(f"/participants/{where[4]}/attributes/age", r["study"], r["pairs"])
        == "/participants/4"
    )
    assert (
        rep.place(f"/analytes/{key}/unit", r["study"], r["pairs"])
        == f"/analytes/{key}/unit"
    )
    assert rep.place("/participants/999", r["study"], r["pairs"]) == "unmatched"
    assert rep.place("/analytes/nowhere/unit", r["study"], r["pairs"]) == "unmatched"


def test_jaccard_of_empty_sets_is_one():
    assert rep.jaccard(set(), set()) == 1.0
    assert rep.mean_pairwise([{1}, {1}, set()]) == 1 / 3


def test_the_command_line_compares_runs(tmp_path):
    reference = load(STUDY)
    for name, seed in (("r1", 1), ("r2", 2)):
        draft, _ = shuffled(reference, seed)
        (tmp_path / name / STUDY).mkdir(parents=True)
        (tmp_path / name / STUDY / "s0.yaml").write_text(
            yaml.safe_dump(draft), encoding="utf-8"
        )
        assert (
            disc.main(
                [
                    str(tmp_path / name),
                    "--out",
                    str(tmp_path / f"{name}.cmp"),
                    "--ref",
                    "HEAD",
                ]
            )
            == 0
        )
        assert (
            mf.main([str(tmp_path / name), "--out", str(tmp_path / f"{name}.cmp")]) == 0
        )
    argv = []
    for name in ("r1", "r2"):
        argv += ["--run", str(tmp_path / name), str(tmp_path / f"{name}.cmp")]
    out = tmp_path / "out"
    assert (
        rep.main(argv + ["--out", str(out), "--ref", "HEAD", "--replicates", "50"]) == 0
    )
    rows = list(
        csv.DictReader((out / "reproducibility_studies.csv").open(encoding="utf-8"))
    )
    assert [
        (r["study_id"], float(r["reproduced_jaccard"]), float(r["recall_sd"]))
        for r in rows
    ] == [(STUDY, 1.0, 0.0)]
    summary = {
        r["metric"]: r
        for r in csv.DictReader((out / "reproducibility.csv").open(encoding="utf-8"))
    }
    assert float(summary["analyte fields the same in every run"]["estimate"]) == 1.0
