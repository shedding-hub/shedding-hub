"""Between-run reproducibility of the agents (benchmark/PLAN.md, Reproducibility).

The studies in `replicates.csv` are each run several times from the same
frozen evidence. Every run is compared with the reference by
`discrepancies.py` and `match_findings.py`; this file compares the runs with
one another, study by study:

- analyte fields: for each scored field of each reference analyte that every
  run matched, whether all runs give the same value; pooled over studies;
- reproduced measurements: the Jaccard similarity of the sets of reference
  measurements two runs reproduce, averaged over the pairs of runs;
- recall against the reference: its standard deviation across runs;
- findings: the Jaccard similarity of the locations of two runs' scored
  findings (check and place), averaged over the pairs of runs.

Drafts number their participants and name their analytes freely, so a
finding's place is translated to the reference's numbering first, through the
same pairing the comparison uses. A place with no partner in the reference
is kept as "unmatched" under its check.

Replicate runs are not adjudicated, so recall here is against the reference.
Jaccard similarity of two empty sets is taken as 1.

Usage (from the repository root):
    python benchmark/reproducibility.py --run RUN_DIR CMP_DIR --run ... --out OUT_DIR
        [--studies replicates.csv] [--ref v1.1.0]
"""

from __future__ import annotations

import argparse
import csv
import itertools
import re
import statistics
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adjudication as adj  # noqa: E402
import align  # noqa: E402
import discrepancies as disc  # noqa: E402
import metrics  # noqa: E402

PARTICIPANT = re.compile(r"^/participants/(\d+)(/.*)?$")
MEASUREMENT = re.compile(r"^(/participants/\d+/measurements/\d+)(/.*)?$")
ANALYTE = re.compile(r"^/analytes/([^/]+)(/.*)?$")

STUDY_COLUMNS = [
    "study_id",
    "runs",
    "fields_compared",
    "fields_same",
    "reproduced_jaccard",
    "recall_by_run",
    "recall_sd",
    "finding_jaccard",
]


def jaccard(a: set, b: set) -> float:
    return 1.0 if not a and not b else len(a & b) / len(a | b)


def mean_pairwise(sets: list[set]) -> float:
    pairs = list(itertools.combinations(sets, 2))
    return sum(jaccard(a, b) for a, b in pairs) / len(pairs) if pairs else 1.0


def field_values(study: adj.Study) -> dict:
    """{(reference analyte, field): the draft's value} for matched analytes."""
    return {
        (k2, name): study.a0[k0].get(name)
        for k0, k2 in study.pairs.items()
        for name in metrics.SCORED_FIELDS
    }


def place(path: str, study: adj.Study, pairs: list[dict]) -> str:
    """A draft path in the reference's numbering, where it has a partner."""
    to_s2 = {p["s0_path"]: p["s2_path"] for p in pairs}
    found = MEASUREMENT.match(path)
    if found and found.group(1) in to_s2:
        return to_s2[found.group(1)] + (found.group(2) or "")
    found = PARTICIPANT.match(path)
    if found:
        partner = study.people.get(int(found.group(1)))
        return "unmatched" if partner is None else f"/participants/{partner}"
    found = ANALYTE.match(path)
    if found:
        partner = study.pairs.get(found.group(1))
        if partner is None:
            return "unmatched"
        return f"/analytes/{partner}{found.group(2) or ''}"
    return path


def compare_runs(runs: list[dict]) -> dict:
    """One study's reproducibility.

    Each run is a dict with `study` (an adjudication.Study, or None when the
    run left no usable draft), `summary`, `pairs` and `findings` (the rows
    of the comparison files for this study).
    """
    values = [field_values(r["study"]) if r["study"] else {} for r in runs]
    shared = set.intersection(*(set(v) for v in values)) if values else set()
    same = sum(len({repr(v[k]) for v in values}) == 1 for k in shared)

    reproduced = [
        {p["s2_path"] for p in r["pairs"] if p["status"] == "agree"} for r in runs
    ]
    recalls = []
    for r in runs:
        truth = int(r["summary"].get("s2_measurements") or 0)
        correct = int(r["summary"].get("correct_measurements") or 0)
        recalls.append(correct / truth if truth else 0.0)

    locations = []
    for r in runs:
        here = set()
        for f in r["findings"]:
            if f["category"] != "data":
                continue
            where = (
                place(f["path"], r["study"], r["pairs"]) if r["study"] else "unmatched"
            )
            here.add((f["check"], where))
        locations.append(here)

    return {
        "runs": len(runs),
        "fields_compared": len(shared),
        "fields_same": same,
        "reproduced_jaccard": mean_pairwise(reproduced),
        "recall_by_run": ";".join(f"{x:.4f}" for x in recalls),
        "recall_sd": statistics.stdev(recalls) if len(recalls) > 1 else 0.0,
        "finding_jaccard": mean_pairwise(locations),
    }


def summarize(rows: list[dict], seed: int, replicates: int) -> list[dict]:
    group = "reproducibility"
    out = [
        metrics.ratio(
            group,
            "analyte fields the same in every run",
            [r["fields_same"] for r in rows],
            [r["fields_compared"] for r in rows],
            seed,
            replicates,
        )
    ]
    for name, column in (
        ("Jaccard of reproduced measurements (median, IQR)", "reproduced_jaccard"),
        ("SD of recall across runs (median, IQR)", "recall_sd"),
        ("Jaccard of finding locations (median, IQR)", "finding_jaccard"),
    ):
        out.append(metrics.describe(group, name, [r[column] for r in rows], len(rows)))
    out.append(
        {
            "group": group,
            "metric": "largest SD of recall across runs",
            "estimate": max((r["recall_sd"] for r in rows), default=float("nan")),
            "studies": len(rows),
        }
    )
    return out


# ------------------------------------------------------------------------ CLI


def load_run(run_dir: Path, cmp_dir: Path, study_id: str, s2: dict, tol) -> dict:
    def mine(name):
        return [
            r for r in metrics.read_csv(cmp_dir / name) if r["study_id"] == study_id
        ]

    summary = mine("comparison_summary.csv")
    draft = run_dir / study_id / "s0.yaml"
    s0 = None
    if draft.is_file():
        try:
            s0 = yaml.load(draft.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader)
        except yaml.YAMLError:
            s0 = None
    usable = isinstance(s0, dict) and summary and summary[0]["draft"] == "compared"
    return {
        "study": adj.Study(s0, s2, tol) if usable else None,
        "summary": summary[0] if summary else {},
        "pairs": mine("measurement_pairs.csv"),
        "findings": mine("finding_matches.csv"),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--run",
        nargs=2,
        action="append",
        type=Path,
        required=True,
        metavar=("RUN_DIR", "CMP_DIR"),
        help="A replicate run and its comparison folder; give two or more",
    )
    ap.add_argument("--studies", type=Path, help="CSV with a study_id column")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ref", default=disc.RELEASE)
    ap.add_argument("--seed", type=int, default=metrics.SEED)
    ap.add_argument("--replicates", type=int, default=metrics.REPLICATES)
    args = ap.parse_args(argv)
    if len(args.run) < 2:
        ap.error("give at least two runs")

    manifest = disc.read_manifest()
    present = [
        {d.name for d in run.iterdir() if d.is_dir() and d.name in manifest}
        for run, _ in args.run
    ]
    studies = sorted(set.intersection(*present))
    if args.studies:
        wanted = {r["study_id"] for r in metrics.read_csv(args.studies)}
        studies = [s for s in studies if s in wanted]
        missing = sorted(wanted - set(studies))
        if missing:
            print(f"NOT IN EVERY RUN {', '.join(missing)}", file=sys.stderr)

    rows = []
    for study_id in studies:
        row = manifest[study_id]
        s2 = disc.load_reference(row["reference_yaml"], args.ref)
        tol = align.tolerances_for(row["input_csv_origin"])
        runs = [load_run(r, c, study_id, s2, tol) for r, c in args.run]
        rows.append({"study_id": study_id, **compare_runs(runs)})

    args.out.mkdir(parents=True, exist_ok=True)
    disc.write_csv(args.out / "reproducibility_studies.csv", STUDY_COLUMNS, rows)
    disc.write_csv(
        args.out / "reproducibility.csv",
        metrics.METRIC_COLUMNS,
        summarize(rows, args.seed, args.replicates),
    )
    print(f"{len(rows)} studies over {len(args.run)} runs -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
