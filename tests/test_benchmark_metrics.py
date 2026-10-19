"""Metrics follow the adjudicated truth, and the bootstrap behaves.

The comparison comes from a draft made by perturbing a copy of a released
dataset, with one error of each kind at a different participant.
"""

import copy
import csv
import io
import math
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import adjudication as adj  # noqa: E402
import align  # noqa: E402
import discrepancies as disc  # noqa: E402
import match_findings as mf  # noqa: E402
import metrics  # noqa: E402

STUDY = "tan2021early"


def load(study):
    path = REPO_ROOT / "data" / study / f"{study}.yaml"
    return yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader)


def numeric_reading(participant):
    """The index of a numeric reading whose time and value are unique."""
    ms = participant["measurements"]
    for j, m in enumerate(ms):
        if not align.is_number(m["value"]):
            continue
        if any(
            k != j and (o["time"] == m["time"] or o["value"] == m["value"])
            for k, o in enumerate(ms)
        ):
            continue
        return j
    return None


def perturbed():
    """A draft with a wrong value, an omission, an extra reading and a
    flipped non-detect, each at its own participant, and a wrong age."""
    reference = load(STUDY)
    draft = copy.deepcopy(reference)
    people = [
        (i, j)
        for i, p in enumerate(reference["participants"])
        if len(p["measurements"]) >= 3 and (j := numeric_reading(p)) is not None
    ]
    (a, ja), (b, jb), (c, _), (d, jd) = people[:4]
    draft["participants"][a]["measurements"][ja]["value"] += 9
    draft["participants"][a]["attributes"]["age"] += 1
    del draft["participants"][b]["measurements"][jb]
    extra = dict(draft["participants"][c]["measurements"][0], time=999, value=21.0)
    draft["participants"][c]["measurements"].append(extra)
    draft["participants"][d]["measurements"][jd]["value"] = "negative"
    return draft, reference, (a, ja)


def as_csv(rows, columns):
    """Rows as the CSV files give them back: every value a string."""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    buffer.seek(0)
    return list(csv.DictReader(buffer))


def comparison():
    draft, reference, where = perturbed()
    result = disc.compare(STUDY, draft, reference, align.REPORTED)
    rows = as_csv(result.rows, disc.COLUMNS)
    summary = as_csv([result.counts], disc.SUMMARY_COLUMNS)[0]
    pairs = as_csv(result.pairs, disc.PAIR_COLUMNS)
    return rows, summary, pairs, where


def counts_with(changes=None, shared=()):
    rows, summary, pairs, _ = comparison()
    verdict = metrics.verdicts(rows, [], {})
    by_type = {r["type"]: r["id"] for r in rows}
    for kind, v in (changes or {}).items():
        verdict[by_type[kind]] = v
    return metrics.study_counts(summary, rows, pairs, verdict, list(shared)), summary


def test_without_adjudication_the_reference_is_truth():
    counts, summary = counts_with()
    truth = int(summary["s2_measurements"])
    assert counts["truth"] == truth
    assert counts["draft_measurements"] == truth  # one dropped, one added
    # The wrong value, the omission and the flipped non-detect.
    assert counts["reproduced"] == truth - 3
    assert counts["attribute_errors"] == 1
    assert counts["analyte_field_errors"] == 0
    assert counts["participant_count_agrees"] and counts["analyte_count_agrees"]


def test_each_verdict_moves_the_counts_as_documented():
    base, _ = counts_with()

    def delta(changes):
        after, _ = counts_with(changes)
        return {
            k: after[k] - base[k]
            for k in ("truth", "draft_measurements", "reproduced", "excluded")
            if after[k] != base[k]
        }

    assert delta({"wrong value": "S0"}) == {"reproduced": 1}
    assert delta({"wrong value": "both"}) == {"reproduced": 1}
    assert delta({"wrong value": "neither"}) == {}
    assert delta({"wrong non-detect status": "curator"}) == {
        "truth": -1,
        "draft_measurements": -1,
        "excluded": 1,
    }
    assert delta({"omission": "S0"}) == {"truth": -1}
    assert delta({"omission": "curator"}) == {"truth": -1, "excluded": 1}
    assert delta({"extra": "S0"}) == {"truth": 1, "reproduced": 1}
    assert delta({"extra": "both"}) == {"draft_measurements": -1}
    after, _ = counts_with({"wrong attribute value": "S0"})
    assert after["attribute_errors"] == 0


def test_answers_map_through_the_key():
    rows, _, _, _ = comparison()
    wrong = next(r for r in rows if r["type"] == "wrong value")
    omitted = next(r for r in rows if r["type"] == "omission")
    key = [
        {
            "packet": "discrepancies",
            "ref_id": wrong["id"],
            "item_id": "D1",
            "A_is": "S2",
            "B_is": "S0",
        },
        {
            "packet": "discrepancies",
            "ref_id": omitted["id"],
            "item_id": "D2",
            "A_is": "S0",
            "B_is": "S2",
        },
    ]
    verdict = metrics.verdicts(rows, key, {"D1": "b", "D2": "a"})
    assert verdict[wrong["id"]] == "S0"
    assert verdict[omitted["id"]] == "S0"
    assert metrics.verdicts(rows, key, {"D1": "a"})[wrong["id"]] == "S2"
    assert metrics.verdicts(rows, key, {})[wrong["id"]] == "reference"
    # Unscored differences get no verdict.
    unscored = [r["id"] for r in rows if r["scored"] == "no"]
    assert not set(unscored) & set(verdict)


def test_outcomes_are_read_loosely_and_checked(tmp_path):
    path = tmp_path / "packet.csv"
    path.write_text(
        "item_id,outcome\nD1, A \nD2,Both\nD3,\nD4,maybe\n", encoding="utf-8"
    )
    answers, problems = metrics.read_answers(path, metrics.DISCREPANCY_OUTCOMES)
    assert answers == {"D1": "a", "D2": "both"}
    assert problems == ["packet.csv D4: 'maybe'"]


def test_a_shared_error_at_a_measurement_is_not_reproduced():
    rows, summary, pairs, _ = comparison()
    agreeing = next(p["s0_path"] for p in pairs if p["status"] == "agree")
    base, _ = counts_with()
    after, _ = counts_with(shared=[agreeing + "/value", "/participants/0/attributes"])
    assert after["reproduced"] == base["reproduced"] - 1
    assert after["shared_errors"] == 2


def test_non_detect_agreement_follows_the_verdict():
    rows, _, pairs, _ = comparison()
    flipped = next(r for r in rows if r["type"] == "wrong non-detect status")
    place = flipped["s2_path"].removesuffix("/value")
    agree, total = metrics.nondetect_agreement(pairs, {})
    assert total - agree == 1
    agree2, total2 = metrics.nondetect_agreement(pairs, {(STUDY, place): "S0"})
    assert (agree2, total2) == (agree + 1, total)
    agree3, total3 = metrics.nondetect_agreement(pairs, {(STUDY, place): "curator"})
    assert (agree3, total3) == (agree, total - 1)


def test_value_errors_are_in_cycles_for_ct():
    _, _, pairs, _ = comparison()
    errors = metrics.value_errors(pairs)
    assert errors["log10"] == []
    assert max(errors["cycles"]) == 9
    assert sorted(errors["cycles"])[:-1] == [0] * (len(errors["cycles"]) - 1)


def test_the_bootstrap_resamples_studies():
    # Identical studies: no spread.
    assert metrics.bootstrap([9, 9, 9], [10, 10, 10]) == (0.9, 0.9)
    # One study: the interval is the point.
    low, high = metrics.bootstrap([5], [10])
    assert low == high == 0.5
    # Unequal studies: the interval holds the pooled estimate, and repeats.
    num, den = [90, 5, 300, 40], [100, 10, 400, 50]
    low, high = metrics.bootstrap(num, den)
    assert low < sum(num) / sum(den) < high
    assert metrics.bootstrap(num, den) == (low, high)
    # With few studies there are few distinct resamples, so seeds can agree;
    # with twelve unequal studies they do not.
    many_num = [90, 5, 300, 40, 7, 66, 12, 250, 3, 81, 140, 19]
    many_den = [100, 10, 400, 50, 9, 70, 30, 260, 8, 90, 200, 20]
    assert metrics.bootstrap(many_num, many_den, seed=1) != metrics.bootstrap(
        many_num, many_den
    )


def test_kappa():
    assert metrics.kappa(["a", "b", "a"], ["a", "b", "a"]) == 1.0
    assert metrics.kappa(["a", "a", "b", "b"], ["a", "b", "b", "b"]) == 0.5


def test_the_command_line_reads_answers_and_reports_kappa(tmp_path):
    draft, _, _ = perturbed()
    run = tmp_path / "run"
    (run / STUDY).mkdir(parents=True)
    (run / STUDY / "s0.yaml").write_text(yaml.safe_dump(draft), encoding="utf-8")
    cmp_dir, adj_dir, out = tmp_path / "cmp", tmp_path / "adj", tmp_path / "out"
    assert disc.main([str(run), "--out", str(cmp_dir), "--ref", "HEAD"]) == 0
    assert mf.main([str(run), "--out", str(cmp_dir)]) == 0
    argv = [str(run), "--comparison", str(cmp_dir), "--out", str(adj_dir)]
    assert adj.main(argv + ["--ref", "HEAD"]) == 0

    def recall(folder):
        rows = list(csv.DictReader((folder / "metrics.csv").open(encoding="utf-8")))
        return (
            next(
                float(r["estimate"])
                for r in rows
                if r["group"] == "A: manual" and r["metric"] == "recall"
            ),
            rows,
        )

    common = ["--comparison", str(cmp_dir), "--out", str(out), "--replicates", "50"]
    assert metrics.main(common) == 0
    before, _ = recall(out)

    # The adjudicator finds the draft right about the wrong value, and the
    # second adjudicator agrees on every rechecked item.
    key = list(csv.DictReader((adj_dir / "key.csv").open(encoding="utf-8")))
    for name in ("discrepancy_packet.csv", "recheck_discrepancy_packet.csv"):
        path = adj_dir / name
        items = list(csv.DictReader(path.open(encoding="utf-8")))
        for item in items:
            k = next(k for k in key if k["item_id"] == item["item_id"])
            draft_label = "A" if k["A_is"] == "S0" else "B"
            reference_label = "B" if draft_label == "A" else "A"
            item["outcome"] = (
                draft_label if k["type"] == "wrong value" else reference_label
            )
        disc.write_csv(path, adj.DISCREPANCY_COLUMNS, items)
    assert metrics.main(common + ["--adjudication", str(adj_dir)]) == 0
    after, rows = recall(out)
    summary = list(
        csv.DictReader((cmp_dir / "comparison_summary.csv").open(encoding="utf-8"))
    )[0]
    truth = int(summary["s2_measurements"])
    assert math.isclose(before, (truth - 3) / truth)
    assert math.isclose(after, (truth - 2) / truth)
    kappas = [r for r in rows if r["group"] == "adjudication"]
    assert [(r["metric"], float(r["estimate"])) for r in kappas] == [
        ("Cohen's kappa, discrepancy", 1.0)
    ]

    # An outcome that is not one of the codes stops the run.
    path = adj_dir / "discrepancy_packet.csv"
    items = list(csv.DictReader(path.open(encoding="utf-8")))
    items[0]["outcome"] = "draft"
    disc.write_csv(path, adj.DISCREPANCY_COLUMNS, items)
    assert metrics.main(common + ["--adjudication", str(adj_dir)]) == 1


def review(verdict_changes=None, answers=None):
    """Review metrics for the perturbed study with five findings."""
    rows, summary, pairs, _ = comparison()
    verdict = metrics.verdicts(rows, [], {})
    by_type = {r["type"]: r for r in rows}
    for kind, v in (verdict_changes or {}).items():
        verdict[by_type[kind]["id"]] = v
    agreeing = [p["s0_path"] for p in pairs if p["status"] == "agree"]
    omitted_at = by_type["omission"]["s0_anchor"].removesuffix("/measurements")

    def finding(id, check, path, severity="warning"):
        return {
            "id": id,
            "check": check,
            "severity": severity,
            "path": path,
            "message": "m",
        }

    findings = [
        finding("f1", "analyte.value_vs_paper", by_type["wrong value"]["s0_path"]),
        finding("f2", "sample.count_vs_paper", omitted_at),
        finding("f3", "analyte.value_vs_paper", agreeing[0] + "/value"),
        finding("f4", "schema.missing_field", "/analytes"),
        finding("f5", "analyte.value_vs_paper", agreeing[1] + "/value", "info"),
    ]
    matched = as_csv(mf.match(STUDY, findings, rows), mf.COLUMNS)
    key = [
        {"packet": "findings", "ref_id": f"{STUDY}:{f}", "item_id": f}
        for f in ("f3", "f5")
    ]
    status = metrics.finding_status(matched, verdict, key, answers or {})
    counts = metrics.study_counts(summary, rows, pairs, verdict, [])
    study = {"study_id": STUDY, **counts}
    out = metrics.review_metrics(
        "A", [study], rows, verdict, status, matched, replicates=50
    )
    return {r["metric"]: r for r in out}, {f["finding_id"]: f["status"] for f in status}


def test_review_sensitivity_counts_errors_hit_by_a_finding():
    out, status = review()
    # Five errors: wrong value, omission, extra, flipped non-detect, wrong age.
    assert (
        out["review sensitivity"]["numerator"],
        out["review sensitivity"]["denominator"],
    ) == (2, 5)
    assert out["review sensitivity | wrong value"]["estimate"] == 1.0
    assert out["review sensitivity | extra"]["estimate"] == 0.0
    # When the draft was right about the value, there is one error fewer, and
    # the finding about it was a false alarm.
    out, status = review({"wrong value": "S0"})
    assert (
        out["review sensitivity"]["numerator"],
        out["review sensitivity"]["denominator"],
    ) == (1, 4)
    assert status["f1"] == "wrong"


def test_review_precision_uses_adjudicated_findings_only():
    out, status = review()
    assert status == {
        "f1": "right",
        "f2": "right",
        "f3": "not adjudicated",
        "f5": "not adjudicated",
    }
    assert (
        out["review precision"]["numerator"],
        out["review precision"]["denominator"],
    ) == (2, 2)
    assert (
        out["scored findings left out of precision | not adjudicated"]["numerator"] == 2
    )

    out, status = review(answers={"f3": "correct", "f5": "incorrect"})
    assert (
        out["review precision"]["numerator"],
        out["review precision"]["denominator"],
    ) == (3, 4)
    without = out["review precision, without info"]
    assert (without["numerator"], without["denominator"]) == (3, 3)
    assert (
        out["errors shared by draft and reference, found by findings"]["numerator"] == 1
    )

    out, status = review({"wrong value": "curator"})
    assert status["f1"] == "curator"


def test_review_workload_counts_every_finding():
    out, _ = review()
    assert out["scored findings per study"]["numerator"] == 4
    assert out["all findings per study"]["numerator"] == 5
    assert out["findings | schema_gap"]["numerator"] == 1
    assert out["findings | data"]["numerator"] == 4
