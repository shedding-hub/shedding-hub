"""List every difference between an agent draft (S0) and its reference (S2).

One row per discrepancy, with the taxonomy of the analysis plan:

    study        missing analyte, extra analyte, participant-count mismatch
    analyte      wrong biomarker / specimen / unit / reference_event /
                 gene_target / limit_of_detection / limit_of_quantification
    participant  missing, extra, wrong attribute value, attribute in the wrong
                 field, missing attribute, extra attribute
    measurement  omission, extra, wrong time, wrong value, wrong non-detect
                 status

A study with no usable draft gets a single study-level row (`no draft` or
`unreadable draft`), since a failed extraction is a result.

A row describes a difference. It does not say which side is right: that is
decided at adjudication.

Each row carries the JSON pointer of the element in S0 and in S2 where it
exists. For something absent from S0 (an omission, a missing participant),
`s0_anchor` is the place in S0 where it would be, so that a review finding
about that place can be related to it. `n_measurements` is the number of
measurements the row stands for.

Usage (from the repository root), on the outputs of a benchmark run:
    python benchmark/discrepancies.py RUN_DIR --out OUT_DIR [--ref v1.1.0]
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import align  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "benchmark" / "manifest.csv"
RELEASE = "v1.1.0"

COLUMNS = [
    "id",
    "study_id",
    "level",
    "type",
    "s0_path",
    "s0_anchor",
    "s2_path",
    "s0_value",
    "s2_value",
    "n_measurements",
]
SUMMARY_COLUMNS = [
    "study_id",
    "draft",
    "s0_analytes",
    "s2_analytes",
    "matched_analytes",
    "s0_participants",
    "s2_participants",
    "matched_participants",
    "s0_measurements",
    "s2_measurements",
    "matched_measurements",
    "correct_measurements",
    "discrepancies",
]


def pointer(*parts) -> str:
    """A JSON pointer (RFC 6901) from path segments."""
    return "".join("/" + str(p).replace("~", "~0").replace("/", "~1") for p in parts)


@dataclass
class Comparison:
    """The discrepancies of one study and the counts behind the metrics."""

    study_id: str
    rows: list[dict] = field(default_factory=list)
    counts: dict = field(default_factory=dict)

    def add(self, level, kind, *, s0=None, anchor=None, s2=None, v0="", v2="", n=0):
        self.rows.append(
            {
                "id": f"{self.study_id}:{len(self.rows) + 1:04d}",
                "study_id": self.study_id,
                "level": level,
                "type": kind,
                "s0_path": s0 or "",
                "s0_anchor": anchor if anchor is not None else (s0 or ""),
                "s2_path": s2 or "",
                "s0_value": v0,
                "s2_value": v2,
                "n_measurements": n,
            }
        )


def _count_measurements(participants, analyte=None) -> int:
    return sum(
        1
        for p in participants
        for m in align.measurements_of(p)
        if analyte is None or m.get("analyte") == analyte
    )


def _compare_attributes(out: Comparison, i0, i2, p0: dict, p2: dict) -> None:
    attrs0 = p0.get("attributes") if isinstance(p0.get("attributes"), dict) else {}
    attrs2 = p2.get("attributes") if isinstance(p2.get("attributes"), dict) else {}
    only0 = {f for f in attrs0 if f not in attrs2}
    only2 = {f for f in attrs2 if f not in attrs0}

    for name in sorted(attrs0.keys() & attrs2.keys()):
        if attrs0[name] != attrs2[name]:
            out.add(
                "participant",
                "wrong attribute value",
                s0=pointer("participants", i0, "attributes", name),
                s2=pointer("participants", i2, "attributes", name),
                v0=attrs0[name],
                v2=attrs2[name],
            )
    # The same value under another name is one misplaced attribute, not a
    # missing one plus an extra one.
    for name2 in sorted(only2):
        twin = next(
            (n for n in sorted(only0) if attrs0[n] == attrs2[name2]),
            None,
        )
        if twin is not None:
            only0.discard(twin)
            out.add(
                "participant",
                "attribute in the wrong field",
                s0=pointer("participants", i0, "attributes", twin),
                s2=pointer("participants", i2, "attributes", name2),
                v0=f"{twin}: {attrs0[twin]}",
                v2=f"{name2}: {attrs2[name2]}",
            )
        else:
            out.add(
                "participant",
                "missing attribute",
                anchor=pointer("participants", i0, "attributes"),
                s2=pointer("participants", i2, "attributes", name2),
                v2=attrs2[name2],
            )
    for name in sorted(only0):
        out.add(
            "participant",
            "extra attribute",
            s0=pointer("participants", i0, "attributes", name),
            v0=attrs0[name],
        )


def compare(study_id: str, s0, s2: dict, tol: align.Tolerances) -> Comparison:
    """Every discrepancy between a draft and its reference."""
    out = Comparison(study_id)
    a0, a2 = align.analytes_of(s0), align.analytes_of(s2)
    p0, p2 = align.participants_of(s0), align.participants_of(s2)

    analyte_pairs, extra_analytes, missing_analytes = align.match_analytes(a0, a2)
    participant_pairs, extra_people, missing_people = align.match_participants(
        p0, p2, analyte_pairs, a2, tol
    )
    matched0 = [p0[i] for i in participant_pairs]
    matched2 = [p2[j] for j in participant_pairs.values()]

    # Study level.
    for key in missing_analytes:
        out.add(
            "study",
            "missing analyte",
            anchor=pointer("analytes"),
            s2=pointer("analytes", key),
            v2=_describe(a2[key]),
            n=_count_measurements(matched2, key),
        )
    for key in extra_analytes:
        out.add(
            "study",
            "extra analyte",
            s0=pointer("analytes", key),
            v0=_describe(a0[key]),
            n=_count_measurements(matched0, key),
        )
    if len(p0) != len(p2):
        out.add(
            "study",
            "participant-count mismatch",
            s0=pointer("participants"),
            s2=pointer("participants"),
            v0=len(p0),
            v2=len(p2),
        )

    # Analyte fields.
    for k0, k2 in analyte_pairs.items():
        for name in align.ANALYTE_FIELDS:
            if a0[k0].get(name) != a2[k2].get(name):
                out.add(
                    "analyte",
                    f"wrong {name}",
                    s0=pointer("analytes", k0, name) if name in a0[k0] else None,
                    anchor=pointer("analytes", k0, name),
                    s2=pointer("analytes", k2, name) if name in a2[k2] else None,
                    v0=a0[k0].get(name, ""),
                    v2=a2[k2].get(name, ""),
                )

    # Participants.
    for j in missing_people:
        out.add(
            "participant",
            "missing",
            anchor=pointer("participants"),
            s2=pointer("participants", j),
            n=_count_measurements([p2[j]]),
        )
    for i in extra_people:
        out.add(
            "participant",
            "extra",
            s0=pointer("participants", i),
            n=_count_measurements([p0[i]]),
        )

    # Attributes and measurements of matched participants.
    matched = correct = 0
    for i0, i2 in participant_pairs.items():
        _compare_attributes(out, i0, i2, p0[i0], p2[i2])
        groups0, groups2 = align.by_analyte(p0[i0]), align.by_analyte(p2[i2])
        list0, list2 = align.measurements_of(p0[i0]), align.measurements_of(p2[i2])

        def at0(index, leaf=None):
            parts = ["participants", i0, "measurements", index]
            return pointer(*parts, *([leaf] if leaf else []))

        def at2(index, leaf=None):
            parts = ["participants", i2, "measurements", index]
            return pointer(*parts, *([leaf] if leaf else []))

        for k0, k2 in analyte_pairs.items():
            unit = a2[k2].get("unit")
            result = align.pair_measurements(
                groups0.get(k0, []), groups2.get(k2, []), unit, tol
            )
            for m0, m2, agree in result.matched:
                matched += 1
                correct += agree
                if agree:
                    continue
                v0, v2 = list0[m0].get("value"), list2[m2].get("value")
                same_kind = align.value_kind(v0) == align.value_kind(v2)
                out.add(
                    "measurement",
                    "wrong value" if same_kind else "wrong non-detect status",
                    s0=at0(m0, "value"),
                    s2=at2(m2, "value"),
                    v0=v0,
                    v2=v2,
                    n=1,
                )
            for m0, m2 in result.shifted:
                out.add(
                    "measurement",
                    "wrong time",
                    s0=at0(m0, "time"),
                    s2=at2(m2, "time"),
                    v0=list0[m0].get("time"),
                    v2=list2[m2].get("time"),
                    n=1,
                )
            for m0 in result.extra:
                out.add("measurement", "extra", s0=at0(m0), v0=_reading(list0[m0]), n=1)
            for m2 in result.missing:
                out.add(
                    "measurement",
                    "omission",
                    anchor=pointer("participants", i0, "measurements"),
                    s2=at2(m2),
                    v2=_reading(list2[m2]),
                    n=1,
                )
        # A draft measurement naming an analyte the draft does not define.
        for key, items in groups0.items():
            if key not in a0:
                for m0, _, _ in items:
                    out.add(
                        "measurement", "extra", s0=at0(m0), v0=_reading(list0[m0]), n=1
                    )

    out.counts = {
        "study_id": study_id,
        "draft": "compared",
        "s0_analytes": len(a0),
        "s2_analytes": len(a2),
        "matched_analytes": len(analyte_pairs),
        "s0_participants": len(p0),
        "s2_participants": len(p2),
        "matched_participants": len(participant_pairs),
        "s0_measurements": _count_measurements(p0),
        "s2_measurements": _count_measurements(p2),
        "matched_measurements": matched,
        "correct_measurements": correct,
        "discrepancies": len(out.rows),
    }
    return out


def no_draft(study_id: str, s2: dict, kind: str) -> Comparison:
    """The single row of a study whose extraction gave nothing usable."""
    out = Comparison(study_id)
    n = _count_measurements(align.participants_of(s2))
    out.add("study", kind, n=n)
    out.counts = {
        "study_id": study_id,
        "draft": kind,
        "s2_analytes": len(align.analytes_of(s2)),
        "s2_participants": len(align.participants_of(s2)),
        "s2_measurements": n,
        "discrepancies": 1,
    }
    return out


def _describe(analyte: dict) -> str:
    return "; ".join(
        f"{name}={analyte[name]}" for name in align.ANALYTE_FIELDS if name in analyte
    )


def _reading(measurement: dict) -> str:
    return f"time={measurement.get('time')}; value={measurement.get('value')}"


# -------------------------------------------------------------------------- CLI


def read_manifest() -> dict:
    with MANIFEST.open(encoding="utf-8", newline="") as fh:
        return {row["study_id"]: row for row in csv.DictReader(fh)}


def load_reference(path: str, ref: str | None) -> dict:
    """A reference dataset: from the release when `ref` is given, else on disk."""
    if ref is None:
        text = (REPO_ROOT / path).read_text(encoding="utf-8")
    else:
        text = subprocess.run(
            ["git", "show", f"{ref}:{path}"],
            cwd=REPO_ROOT,
            capture_output=True,
            check=True,
        ).stdout.decode("utf-8")
    return yaml.load(text, Loader=yaml.CSafeLoader)


def compare_study(study_dir: Path, row: dict, ref: str | None) -> Comparison:
    study_id = row["study_id"]
    s2 = load_reference(row["reference_yaml"], ref)
    draft = study_dir / "s0.yaml"
    if not draft.is_file():
        return no_draft(study_id, s2, "no draft")
    try:
        s0 = yaml.load(draft.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader)
    except yaml.YAMLError:
        s0 = None
    if not isinstance(s0, dict):
        return no_draft(study_id, s2, "unreadable draft")
    return compare(study_id, s0, s2, align.tolerances_for(row["input_csv_origin"]))


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir", type=Path, help="benchmark/runs/<run_id>")
    ap.add_argument("--out", type=Path, required=True, help="Folder for the tables")
    ap.add_argument("--ref", default=RELEASE, help="Release to read S2 from")
    args = ap.parse_args(argv)

    manifest = read_manifest()
    studies = sorted(
        d for d in args.run_dir.iterdir() if d.is_dir() and d.name in manifest
    )
    if not studies:
        raise SystemExit(f"no study folders in {args.run_dir}")

    rows, summary = [], []
    for study_dir in studies:
        result = compare_study(study_dir, manifest[study_dir.name], args.ref)
        rows.extend(result.rows)
        summary.append(result.counts)

    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "discrepancies.csv", COLUMNS, rows)
    write_csv(args.out / "comparison_summary.csv", SUMMARY_COLUMNS, summary)
    print(f"{len(studies)} studies, {len(rows)} discrepancies -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
