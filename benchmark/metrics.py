"""Accuracy and agreement metrics for Analyses A and B (benchmark/PLAN.md).

Truth is the reference (S2) corrected by adjudication. Each scored
discrepancy gets a verdict:

    S0          the draft is right; the reference was wrong
    both        both versions are acceptable
    S2          the reference is right; the draft made an error
    neither     both are wrong; the draft made an error
    curator     a curator's decision the sources cannot settle; left out of
                the counts and reported apart
    reference   not adjudicated (not sampled, or not answered yet); the
                reference stands, so the difference counts as an error

`S0` and `both` mean the difference is not an error of the draft.

Measurement counts per study start from the comparison (truth = S2
measurements, draft = S0 measurements, reproduced = pairs that agree) and
change with each verdict:

    wrong value, time or non-detect status   S0/both: reproduced + 1
                                             curator: truth - 1, draft - 1
    omission (reference only)                S0/both/curator: truth - 1
    extra (draft only)                       S0: truth + 1, reproduced + 1
                                             both/curator: draft - 1
    missing participant or analyte           as omission, for n measurements
    extra participant or analyte             as extra, for n measurements

A finding judged correct where the two versions agree is an error they
share. When it points at one measurement, that measurement no longer counts
as reproduced. Other shared errors are counted apart.

Recall is reproduced / truth, precision is reproduced / draft, both pooled
over the studies of a group, with a 95% interval from a cluster bootstrap
over studies (2,000 replicates, seed 20261008). Field, attribute and count
agreement are pooled the same way. Per-study recall is also summarized
(median and interquartile range); that summary is descriptive and was not
pre-specified. Value errors are over pairs matched on time and use the
reference values as recorded.

When the recheck packets are filled in, Cohen's kappa between the two
adjudicators is reported for each packet.

Usage (from the repository root):
    python benchmark/metrics.py --comparison CMP_DIR [--adjudication ADJ_DIR]
        --out OUT_DIR
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import align  # noqa: E402
import discrepancies as disc  # noqa: E402

SEED = 20261008
REPLICATES = 2000
NOT_ERRORS = {"S0", "both"}
EXCLUDED = {"curator"}
ERRORS = {"S2", "neither", "reference"}
NOT_ADJUDICATED = {"no draft", "unreadable draft"}
SCORED_FIELDS = tuple(f for f in align.ANALYTE_FIELDS if f != "gene_target")
ATTRIBUTE_TYPES = {
    "wrong attribute value",
    "attribute in the wrong field",
    "missing attribute",
}
TIMED_TYPES = {"wrong value", "wrong time", "wrong non-detect status"}
DISCREPANCY_OUTCOMES = {"a", "b", "both", "neither", "curator"}
FINDING_OUTCOMES = {"correct", "incorrect", "cannot tell"}
MEASUREMENT = re.compile(
    r"^/participants/(\d+)/measurements/(\d+)(?:/(?:value|time))?$"
)

METRIC_COLUMNS = [
    "group",
    "metric",
    "numerator",
    "denominator",
    "estimate",
    "ci_low",
    "ci_high",
    "studies",
]
STUDY_COLUMNS = [
    "study_id",
    "era",
    "input_csv_origin",
    "pipeline_version",
    "extraction_model",
    "draft",
    "truth",
    "draft_measurements",
    "reproduced",
    "recall",
    "excluded",
    "shared_errors",
    "analyte_pairs",
    "analyte_field_errors",
    "attributes",
    "attribute_errors",
    "participant_count_agrees",
    "analyte_count_agrees",
]
VERDICT_COLUMNS = ["id", "study_id", "level", "type", "n_measurements", "verdict"]


# ------------------------------------------------------------------ verdicts


def read_answers(path: Path, allowed: set) -> tuple[dict, list[str]]:
    """{item id: outcome} from a filled packet, and any outcome not allowed."""
    if not path.is_file():
        return {}, []
    answers, problems = {}, []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            outcome = (row.get("outcome") or "").strip().lower()
            if not outcome:
                continue
            if outcome not in allowed:
                problems.append(f"{path.name} {row['item_id']}: {outcome!r}")
                continue
            answers[row["item_id"]] = outcome
    return answers, problems


def verdicts(rows: list[dict], key: list[dict], answers: dict) -> dict:
    """{discrepancy id: verdict} for every scored, comparable discrepancy."""
    by_ref = {k["ref_id"]: k for k in key if k["packet"] == "discrepancies"}
    out = {}
    for row in rows:
        if row["scored"] != "yes" or row["type"] in NOT_ADJUDICATED:
            continue
        item = by_ref.get(row["id"])
        outcome = answers.get(item["item_id"]) if item else None
        if outcome in ("a", "b"):
            out[row["id"]] = item["A_is" if outcome == "a" else "B_is"]
        elif outcome in ("both", "neither", "curator"):
            out[row["id"]] = outcome
        else:
            out[row["id"]] = "reference"
    return out


def shared_errors(finding_rows: list[dict], key: list[dict], answers: dict) -> dict:
    """{study: [S0 paths of findings judged correct that hit nothing]}."""
    by_ref = {k["ref_id"]: k for k in key if k["packet"] == "findings"}
    out = defaultdict(list)
    for row in finding_rows:
        item = by_ref.get(f"{row['study_id']}:{row['finding_id']}")
        if item and answers.get(item["item_id"]) == "correct":
            out[row["study_id"]].append(row["path"])
    return out


# ------------------------------------------------------------------- studies


def study_counts(
    summary: dict,
    rows: list[dict],
    pairs: list[dict],
    verdict: dict,
    shared: list[str],
) -> dict:
    """The truth-adjusted counts of one study."""
    if summary["draft"] != "compared":
        truth = int(summary["s2_measurements"])
        return {
            "draft": summary["draft"],
            "truth": truth,
            "draft_measurements": 0,
            "reproduced": 0,
            "excluded": 0,
            "shared_errors": 0,
            "analyte_pairs": 0,
            "analyte_field_errors": 0,
            "attributes": 0,
            "attribute_errors": 0,
            "participant_count_agrees": False,
            "analyte_count_agrees": False,
        }
    truth = int(summary["s2_measurements"])
    draft = int(summary["s0_measurements"])
    reproduced = int(summary["correct_measurements"])
    excluded = field_errors = attribute_errors = 0
    analyte_pairs = int(summary["matched_analytes"])
    attributes = int(summary["s2_attributes"])
    fields_excluded = 0
    participants_agree = summary["s0_participants"] == summary["s2_participants"]
    analytes_agree = summary["s0_analytes"] == summary["s2_analytes"]
    analyte_set_errors = 0

    for row in rows:
        v = verdict.get(row["id"])
        if v is None:
            continue
        kind, n = row["type"], int(row["n_measurements"] or 0)
        fine, out = v in NOT_ERRORS, v in EXCLUDED
        if kind in TIMED_TYPES:
            if fine:
                reproduced += 1
            elif out:
                truth, draft, excluded = truth - 1, draft - 1, excluded + 1
        elif kind in ("omission", "missing", "missing analyte"):
            if fine or out:
                truth -= n
                excluded += n if out else 0
        elif kind in ("extra", "extra analyte"):
            if v == "S0":
                truth, reproduced = truth + n, reproduced + n
            elif v == "both" or out:
                draft -= n
                excluded += n if out else 0
        elif kind == "participant-count mismatch":
            participants_agree = fine or out
        if kind.startswith("wrong ") and row["level"] == "analyte":
            if out:
                fields_excluded += 1
            elif not fine:
                field_errors += 1
        if kind in ATTRIBUTE_TYPES:
            if out:
                attributes -= 1
            elif not fine:
                attribute_errors += 1
        if kind in ("missing analyte", "extra analyte") and not (fine or out):
            analyte_set_errors += 1

    if not analytes_agree and analyte_set_errors == 0:
        analytes_agree = True  # every analyte difference was judged not an error

    agreeing = {p["s0_path"] for p in pairs if p["status"] == "agree"}
    located = set()
    for path in shared:
        found = MEASUREMENT.match(path)
        if found:
            place = f"/participants/{found.group(1)}/measurements/{found.group(2)}"
            if place in agreeing and place not in located:
                located.add(place)
                reproduced -= 1
    return {
        "draft": "compared",
        "truth": truth,
        "draft_measurements": draft,
        "reproduced": reproduced,
        "excluded": excluded,
        "shared_errors": len(shared),
        "analyte_pairs": analyte_pairs * len(SCORED_FIELDS) - fields_excluded,
        "analyte_field_errors": field_errors,
        "attributes": attributes,
        "attribute_errors": attribute_errors,
        "participant_count_agrees": participants_agree,
        "analyte_count_agrees": analytes_agree,
    }


# ------------------------------------------------------------------- metrics


def bootstrap(
    numerators,
    denominators,
    seed: int = SEED,
    replicates: int = REPLICATES,
) -> tuple[float, float]:
    """95% percentile interval of sum(num) / sum(den), resampling studies."""
    num = np.asarray(numerators, dtype=float)
    den = np.asarray(denominators, dtype=float)
    if len(num) == 0 or den.sum() == 0:
        return math.nan, math.nan
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(num), size=(replicates, len(num)))
    totals = den[draws].sum(axis=1)
    ratios = np.divide(
        num[draws].sum(axis=1),
        totals,
        out=np.full(replicates, np.nan),
        where=totals > 0,
    )
    low, high = np.nanpercentile(ratios, [2.5, 97.5])
    return float(low), float(high)


def ratio(group, metric, num, den, seed, replicates) -> dict:
    total_num, total_den = sum(num), sum(den)
    low, high = bootstrap(num, den, seed, replicates)
    return {
        "group": group,
        "metric": metric,
        "numerator": total_num,
        "denominator": total_den,
        "estimate": total_num / total_den if total_den else math.nan,
        "ci_low": low,
        "ci_high": high,
        "studies": len(num),
    }


def describe(group, metric, values, studies: int) -> dict:
    """Median and interquartile range, in the place of an interval.

    `denominator` is the number of values summarized.
    """
    values = [v for v in values if not math.isnan(v)]
    out = {"group": group, "metric": metric, "denominator": len(values)}
    out["studies"] = studies
    if values:
        q1, median, q3 = np.percentile(values, [25, 50, 75])
        out.update(estimate=float(median), ci_low=float(q1), ci_high=float(q3))
    return out


def value_errors(pairs: list[dict]) -> dict:
    """Absolute differences of matched numeric readings, by kind of unit.

    Concentrations in log10 units, cycle thresholds in cycles.
    """
    out = {"log10": [], "cycles": []}
    for p in pairs:
        if p["status"] not in ("agree", "differs"):
            continue
        try:
            v0, v2 = float(p["value0"]), float(p["value2"])
        except ValueError:
            continue
        if p["unit"] == align.CT_UNIT:
            out["cycles"].append(abs(v0 - v2))
        elif v0 > 0 and v2 > 0:
            out["log10"].append(abs(math.log10(v0) - math.log10(v2)))
    return out


def detected(value: str) -> bool | None:
    """Whether a reading as written in the pairs file reports a detection."""
    if value == "negative":
        return False
    if value in align.POSITIVE:
        return True
    try:
        return float(value) > 0
    except ValueError:
        return None


def nondetect_agreement(pairs: list[dict], verdict_by_s2: dict) -> tuple[int, int]:
    """(agreeing, compared) non-detect status over pairs matched on time."""
    agree = compared = 0
    for p in pairs:
        if p["status"] not in ("agree", "differs"):
            continue
        d0, d2 = detected(p["value0"]), detected(p["value2"])
        if d0 is None or d2 is None:
            continue
        v = verdict_by_s2.get((p["study_id"], p["s2_path"]))
        if d0 != d2 and v in EXCLUDED:
            continue
        compared += 1
        agree += d0 == d2 or v in NOT_ERRORS
    return agree, compared


def group_metrics(
    group: str,
    studies: list[dict],
    pairs: list[dict],
    verdict_by_s2: dict,
    seed: int = SEED,
    replicates: int = REPLICATES,
) -> list[dict]:
    """The metrics of section 4 for one group of studies."""
    out = []
    if not studies:
        return out
    r = [s["reproduced"] for s in studies]
    out.append(
        ratio(group, "recall", r, [s["truth"] for s in studies], seed, replicates)
    )
    compared = [s for s in studies if s["draft"] == "compared"]
    out.append(
        ratio(
            group,
            "precision",
            [s["reproduced"] for s in compared],
            [s["draft_measurements"] for s in compared],
            seed,
            replicates,
        )
    )
    out.append(
        describe(
            group,
            "recall per study (median, IQR)",
            [s["reproduced"] / s["truth"] if s["truth"] else math.nan for s in studies],
            len(studies),
        )
    )
    out.append(
        ratio(
            group,
            "analyte field agreement",
            [s["analyte_pairs"] - s["analyte_field_errors"] for s in compared],
            [s["analyte_pairs"] for s in compared],
            seed,
            replicates,
        )
    )
    out.append(
        ratio(
            group,
            "attribute agreement",
            [s["attributes"] - s["attribute_errors"] for s in compared],
            [s["attributes"] for s in compared],
            seed,
            replicates,
        )
    )
    for name in ("participant_count_agrees", "analyte_count_agrees"):
        agree = [int(s[name]) for s in studies]
        out.append(
            ratio(
                group, name.replace("_", " "), agree, [1] * len(agree), seed, replicates
            )
        )
    ids = {s["study_id"] for s in studies}
    mine = [p for p in pairs if p["study_id"] in ids]
    agree, total = nondetect_agreement(mine, verdict_by_s2)
    out.append(
        {
            "group": group,
            "metric": "non-detect status agreement",
            "numerator": agree,
            "denominator": total,
            "estimate": agree / total if total else math.nan,
            "studies": len(ids),
        }
    )
    errors = value_errors(mine)
    for unit in ("log10", "cycles"):
        out.append(
            describe(
                group, f"value error, {unit} (median, IQR)", errors[unit], len(ids)
            )
        )
    return out


def kappa(first: list[str], second: list[str]) -> float:
    """Cohen's kappa for two raters over the same items."""
    n = len(first)
    if n == 0:
        return math.nan
    observed = sum(a == b for a, b in zip(first, second)) / n
    c1, c2 = Counter(first), Counter(second)
    expected = sum(c1[k] * c2[k] for k in c1.keys() | c2.keys()) / n**2
    return 1.0 if expected == 1 else (observed - expected) / (1 - expected)


def groups(studies: list[dict]) -> dict:
    """The groups to report: each era, and its strata."""
    out = {}
    eras = {"A: manual": "manual", "B: ai_assisted": "ai_assisted"}
    for label, era in eras.items():
        members = [s for s in studies if s["era"] == era]
        out[label] = members
        strata = ["input_csv_origin"]
        if era == "ai_assisted":
            strata += ["pipeline_version", "extraction_model"]
        for column in strata:
            for value in sorted({s[column] for s in members}):
                out[f"{label} | {column}={value}"] = [
                    s for s in members if s[column] == value
                ]
    return out


# ------------------------------------------------------------------------ CLI


def read_csv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--comparison", type=Path, required=True)
    ap.add_argument("--adjudication", type=Path, help="Folder with filled packets")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--replicates", type=int, default=REPLICATES)
    args = ap.parse_args(argv)

    rows = read_csv(args.comparison / "discrepancies.csv")
    summary = read_csv(args.comparison / "comparison_summary.csv")
    pairs = read_csv(args.comparison / "measurement_pairs.csv")
    finding_rows = read_csv(args.comparison / "finding_matches.csv")
    key, answers, problems, rechecks = [], {}, [], []
    if args.adjudication:
        key = read_csv(args.adjudication / "key.csv")
        d_answers, p1 = read_answers(
            args.adjudication / "discrepancy_packet.csv", DISCREPANCY_OUTCOMES
        )
        f_answers, p2 = read_answers(
            args.adjudication / "findings_packet.csv", FINDING_OUTCOMES
        )
        answers, problems = {**d_answers, **f_answers}, p1 + p2
        for packet, allowed in (
            ("discrepancy_packet.csv", DISCREPANCY_OUTCOMES),
            ("findings_packet.csv", FINDING_OUTCOMES),
        ):
            again, p3 = read_answers(args.adjudication / f"recheck_{packet}", allowed)
            problems += p3
            both = sorted(i for i in again if i in answers)
            if both:
                rechecks.append(
                    {
                        "group": "adjudication",
                        "metric": f"Cohen's kappa, {packet.split('_')[0]}",
                        "estimate": kappa(
                            [answers[i] for i in both], [again[i] for i in both]
                        ),
                        "denominator": len(both),
                    }
                )
    if problems:
        for problem in problems:
            print(f"UNKNOWN OUTCOME {problem}", file=sys.stderr)
        return 1

    verdict = verdicts(rows, key, answers)
    shared = shared_errors(finding_rows, key, answers)
    manifest = disc.read_manifest()
    by_study = defaultdict(list)
    for row in rows:
        by_study[row["study_id"]].append(row)
    pairs_by_study = defaultdict(list)
    for p in pairs:
        pairs_by_study[p["study_id"]].append(p)

    studies = []
    for s in summary:
        sid = s["study_id"]
        counts = study_counts(
            s, by_study[sid], pairs_by_study[sid], verdict, shared.get(sid, [])
        )
        info = manifest[sid]
        studies.append(
            {
                "study_id": sid,
                **{
                    c: info[c]
                    for c in (
                        "era",
                        "input_csv_origin",
                        "pipeline_version",
                        "extraction_model",
                    )
                },
                **counts,
                "recall": (
                    counts["reproduced"] / counts["truth"] if counts["truth"] else ""
                ),
            }
        )

    # Pairs name a whole measurement; a discrepancy names its value or time.
    verdict_by_s2 = {
        (r["study_id"], re.sub(r"/(value|time)$", "", r["s2_path"])): verdict[r["id"]]
        for r in rows
        if r["id"] in verdict and r["s2_path"]
    }
    metrics = []
    for label, members in groups(studies).items():
        metrics += group_metrics(
            label, members, pairs, verdict_by_s2, args.seed, args.replicates
        )

    metrics += rechecks
    args.out.mkdir(parents=True, exist_ok=True)
    disc.write_csv(args.out / "metrics.csv", METRIC_COLUMNS, metrics)
    disc.write_csv(args.out / "studies.csv", STUDY_COLUMNS, studies)
    disc.write_csv(
        args.out / "verdicts.csv",
        VERDICT_COLUMNS,
        [
            {
                **{c: r[c] for c in VERDICT_COLUMNS if c != "verdict"},
                "verdict": verdict[r["id"]],
            }
            for r in rows
            if r["id"] in verdict
        ],
    )
    tally = Counter(verdict.values())
    print(
        f"{len(studies)} studies, {len(verdict)} scored discrepancies "
        f"({', '.join(f'{k} {v}' for k, v in sorted(tally.items()))}) -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
