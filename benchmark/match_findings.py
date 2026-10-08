"""Relate the review agent's findings to the discrepancies of a study.

A finding hits a discrepancy when two things hold:

- location: the finding's path is the discrepancy's place in S0, or an
  ancestor of it. For something absent from S0 the place is where it would
  be (`s0_anchor` in `discrepancies.csv`);
- type: the finding's check is one that can speak about that kind of
  discrepancy (`CHECK_TO_TYPES`).

Findings from `schema.*` checks are schema gaps and `model.observation`
findings are notes. Both are listed but not scored against data errors, and
neither are findings from a check this file has no mapping for.

A scored finding that hits no discrepancy is marked for adjudication. The
comparison only sees where the draft and the reference differ, so it cannot
tell a false alarm from a finding about an error they share.

Usage (from the repository root):
    python benchmark/match_findings.py RUN_DIR --out OUT_DIR

`OUT_DIR` must already hold `discrepancies.csv` for the same run.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import discrepancies as disc  # noqa: E402

ANALYTE_FIELD_ERRORS = {f"wrong {name}" for name in disc.align.ANALYTE_FIELDS}
ATTRIBUTE_ERRORS = {
    "wrong attribute value",
    "attribute in the wrong field",
    "missing attribute",
    "extra attribute",
}
VALUE_ERRORS = {"wrong value", "wrong non-detect status"}

# The first seven rows are the mapping of the analysis plan, with one addition:
# the review agent also uses `analyte.value_vs_paper` for a single reading
# that differs from the paper, so that check can speak about a measurement's
# value or time as well as an analyte field. The `value.*` rows are rule-based
# checks that the plan does not list; they speak about a measurement's value.
CHECK_TO_TYPES = {
    "analyte.value_vs_paper": ANALYTE_FIELD_ERRORS | VALUE_ERRORS | {"wrong time"},
    "attr.value_vs_paper": ATTRIBUTE_ERRORS,
    "attr.misassigned_vs_paper": ATTRIBUTE_ERRORS,
    "analyte.missing_vs_paper": {"missing analyte"},
    "participant.count_vs_paper": {"participant-count mismatch", "missing", "extra"},
    "sample.count_vs_paper": {"omission", "extra"},
    "sample.types_vs_paper": {"omission", "extra"},
    "value.censored_string": VALUE_ERRORS,
    "value.nonpositive": VALUE_ERRORS,
    "value.nan": VALUE_ERRORS,
}
# `missing` and `extra` name both a participant and a measurement type, so the
# level is part of the test.
CHECK_LEVELS = {
    "participant.count_vs_paper": {"study", "participant"},
    "sample.count_vs_paper": {"measurement"},
    "sample.types_vs_paper": {"measurement"},
}
FINDING_FIELDS = ("check", "severity", "path", "message")
SEVERITIES = {"error", "warning", "info"}

COLUMNS = [
    "study_id",
    "finding_id",
    "check",
    "severity",
    "path",
    "category",
    "hit",
    "matched_discrepancies",
    "same_location_other_type",
    "adjudicate",
]


def category(check: str) -> str:
    """data (scored), schema_gap, observation, or unmapped."""
    if check.startswith("schema."):
        return "schema_gap"
    if check == "model.observation":
        return "observation"
    return "data" if check in CHECK_TO_TYPES else "unmapped"


def covers(finding_path: str, place: str) -> bool:
    """True if the finding's path is `place` or one of its ancestors."""
    finding_path = finding_path.rstrip("/")
    return place == finding_path or place.startswith(finding_path + "/")


def compatible(check: str, row: dict) -> bool:
    return row["type"] in CHECK_TO_TYPES.get(check, ()) and row["level"] in (
        CHECK_LEVELS.get(check, {row["level"]})
    )


def match(study_id: str, findings: list[dict], rows: list[dict]) -> list[dict]:
    """One output row per finding of one study."""
    out = []
    for finding in findings:
        check, path = str(finding.get("check", "")), str(finding.get("path", ""))
        here = [r for r in rows if r["s0_anchor"] and covers(path, r["s0_anchor"])]
        kind = category(check)
        hits = [r["id"] for r in here if compatible(check, r)] if kind == "data" else []
        others = [r["id"] for r in here if r["id"] not in hits]
        out.append(
            {
                "study_id": study_id,
                "finding_id": finding.get("id", ""),
                "check": check,
                "severity": finding.get("severity", ""),
                "path": path,
                "category": kind,
                "hit": "yes" if hits else "no",
                "matched_discrepancies": ";".join(hits),
                "same_location_other_type": ";".join(others),
                # A scored finding that hits nothing may still be right: the
                # draft and the reference can share an error. An adjudicator
                # decides, so that such a finding is not counted as a false
                # alarm by default.
                "adjudicate": "yes" if kind == "data" and not hits else "no",
            }
        )
    return out


def resolves(document, path: str) -> bool:
    """True if a JSON pointer addresses something in the document."""
    node = document
    for part in filter(None, path.split("/")):
        part = part.replace("~1", "/").replace("~0", "~")
        if isinstance(node, dict) and part in node:
            node = node[part]
        elif isinstance(node, list) and part.isdigit() and int(part) < len(node):
            node = node[int(part)]
        else:
            return False
    return True


def read_findings(path: Path) -> tuple[list[dict], list[str]]:
    """The findings of a review file and anything wrong with their form."""
    problems = []
    try:
        review = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        return [], [f"not YAML: {exc}"]
    findings = review.get("findings") if isinstance(review, dict) else None
    if not isinstance(findings, list):
        return [], ["no `findings` list"]
    for number, finding in enumerate(findings):
        if not isinstance(finding, dict):
            problems.append(f"finding {number} is not a mapping")
            continue
        for name in FINDING_FIELDS:
            if not isinstance(finding.get(name), str):
                problems.append(f"finding {number} has no text `{name}`")
        if finding.get("severity") not in SEVERITIES:
            problems.append(
                f"finding {number} has severity {finding.get('severity')!r}"
            )
    return [f for f in findings if isinstance(f, dict)], problems


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir", type=Path, help="benchmark/runs/<run_id>")
    ap.add_argument("--out", type=Path, required=True, help="Folder for the tables")
    args = ap.parse_args(argv)

    with (args.out / "discrepancies.csv").open(encoding="utf-8", newline="") as fh:
        by_study: dict = {}
        for row in csv.DictReader(fh):
            by_study.setdefault(row["study_id"], []).append(row)

    manifest = disc.read_manifest()
    out, problems, reviewed = [], [], 0
    for study_dir in sorted(args.run_dir.iterdir()):
        review = study_dir / "findings.yaml"
        if study_dir.name not in manifest or not review.is_file():
            continue
        reviewed += 1
        findings, wrong = read_findings(review)
        problems += [f"{study_dir.name}: {p}" for p in wrong]
        out += match(study_dir.name, findings, by_study.get(study_dir.name, []))

    disc.write_csv(args.out / "finding_matches.csv", COLUMNS, out)
    hits = sum(1 for row in out if row["hit"] == "yes")
    print(f"{reviewed} review files, {len(out)} findings, {hits} hits -> {args.out}")
    for problem in problems:
        print(f"MALFORMED {problem}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
