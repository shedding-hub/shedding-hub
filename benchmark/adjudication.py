"""Build the blinded packets for adjudication (section 8 of benchmark/PLAN.md).

Two packets come out of one benchmark run:

- discrepancies: scored differences between the draft (S0) and the reference
  (S2). Each shows the study, where to look and the two values labelled A and
  B in random order. Nothing in it says which version came from the agent,
  and no finding is shown. The type of difference is worded so that it does
  not give the direction away: an omission and an extra measurement are both
  "a measurement in one version only".
- findings: scored findings of the review agent that hit no discrepancy, to
  be judged against the source. These quote the draft, so they are done after
  the discrepancy packet.

Which items go in, as the plan says: every scored discrepancy of the
manual-era studies, or a random sample of 300 if there are more, and a random
sample of 200 from the AI era; every unmatched finding of the manual era and
a random sample of 200 from the AI era. A random 10% of each packet, drawn
within each era, goes to the second adjudicator in a copy with no answers.

One kind of item cannot be fully blinded. A missing attribute is always
empty in the draft, because an attribute that only the draft records is not
scored and so never reaches the packet.

`no draft` and `unreadable draft` rows are not adjudicated: a failed
extraction is a result, with nothing to compare.

The draws use one seed. The plan fixes a seed only for the bootstrap, and
the same value is used here so that the packets can be rebuilt exactly.

Usage (from the repository root):
    python benchmark/adjudication.py RUN_DIR --comparison CMP_DIR --out OUT_DIR

`CMP_DIR` holds `discrepancies.csv` and `finding_matches.csv` for RUN_DIR.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import align  # noqa: E402
import discrepancies as disc  # noqa: E402
import match_findings as mf  # noqa: E402

SEED = 20261008
MANUAL_CAP = 300  # every manual-era discrepancy, unless there are more
AI_SAMPLE = 200
RECHECK_SHARE = 0.10
NOT_ADJUDICATED = {"no draft", "unreadable draft"}
PRIMARY = "Yuke Wang"
SECOND = "Youwei Hu"

ABSENT = "(not in this version)"
NOT_GIVEN = "(not given)"
LOCATING_FIELDS = ("biomarker", "specimen", "unit", "reference_event")

DISCREPANCY_COLUMNS = [
    "item_id",
    "study_id",
    "source",
    "difference",
    "where",
    "A",
    "B",
    "outcome",
    "correct_value",
    "note",
]
FINDING_COLUMNS = [
    "item_id",
    "study_id",
    "source",
    "check",
    "severity",
    "confidence",
    "path",
    "where",
    "draft_content",
    "message",
    "suggested_patch",
    "outcome",
    "note",
]
KEY_COLUMNS = [
    "item_id",
    "packet",
    "ref_id",
    "study_id",
    "era",
    "type",
    "A_is",
    "B_is",
    "recheck",
]

INSTRUCTIONS = """\
Adjudication packets for benchmark run {run}.
Analysis plan: benchmark/PLAN.md, tag benchmark-plan-v1, section 8.

Primary adjudicator: {primary}. Second adjudicator (recheck): {second}.

Do discrepancy_packet.csv first and findings_packet.csv after it. The findings
quote the draft, so reading them first would show which version is which.

discrepancy_packet.csv
  Each row is a place where two versions of a dataset differ. Judge against
  the article and its supplementary files only. Do not open the released
  dataset or the run outputs. "{absent}" means the item is missing from that
  version; "{not_given}" means the field or attribute is left empty.
  In `outcome` write one of:
    A        A is right
    B        B is right
    both     both are acceptable
    neither  both are wrong; write the right value in `correct_value`
    curator  a curator's decision that the sources cannot settle
  Use `note` for anything that explains the answer.

findings_packet.csv
  Each row is a review finding that the comparison could not confirm, because
  the two versions agree at that place. They may share an error. In `outcome`
  write correct, incorrect or cannot tell.

recheck_discrepancy_packet.csv, recheck_findings_packet.csv
  A random 10% of the items, for the second adjudicator, without the first
  adjudicator's answers. Fill them in independently.

key.csv and packet_info.json are not for adjudicators: the key says which
version each label stands for.
"""


# ------------------------------------------------------------------- sampling


def sample(items: list, n: int | None, seed: int, label: str) -> list:
    """All items if there are at most `n` (or `n` is None), else a random `n`."""
    if n is None or len(items) <= n:
        return list(items)
    return random.Random(f"{seed}:{label}").sample(items, n)


def recheck_ids(items: list[dict], seed: int, packet: str) -> set:
    """A random share of the items of each era, rounded up."""
    chosen = set()
    for era in sorted({item["era"] for item in items}):
        ids = sorted(item["item_id"] for item in items if item["era"] == era)
        k = math.ceil(RECHECK_SHARE * len(ids))
        chosen.update(random.Random(f"{seed}:recheck:{packet}:{era}").sample(ids, k))
    return chosen


def select(rows: list[dict], era_of: dict, seed: int, packet: str, eligible) -> list:
    """The rows that go into a packet, in random order."""
    picked = []
    for era, n in (
        ("manual", MANUAL_CAP if packet == "D" else None),
        ("ai_assisted", AI_SAMPLE),
    ):
        pool = sorted(
            (r for r in rows if era_of.get(r["study_id"]) == era and eligible(r)),
            key=lambda r: r["_ref_id"],
        )
        picked += sample(pool, n, seed, f"{packet}:{era}")
    random.Random(f"{seed}:order:{packet}").shuffle(picked)
    return picked


# ------------------------------------------------------------------ rendering


def text(value) -> str:
    """A value as the adjudicator reads it."""
    if value is None or value == "":
        return NOT_GIVEN
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, (list, tuple)):
        return ", ".join(text(v) for v in value)
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True, default=str)
    return str(value)


def parts(path: str) -> list[str]:
    """The segments of a JSON pointer."""
    return [p.replace("~1", "/").replace("~0", "~") for p in path.split("/")[1:]]


def node_at(document, path: str):
    node = document
    for part in parts(path):
        if isinstance(node, dict):
            node = node.get(part)
        elif isinstance(node, list) and part.isdigit() and int(part) < len(node):
            node = node[int(part)]
        else:
            return None
    return node


def index_after(path: str, name: str) -> int | None:
    """The list index that follows `name` in a pointer, if any."""
    segments = parts(path)
    for here, after in zip(segments, segments[1:]):
        if here == name and after.isdigit():
            return int(after)
    return None


def reading(measurement: dict) -> str:
    return (
        f"time={text(measurement.get('time'))}; value={text(measurement.get('value'))}"
    )


def attributes_of(participant: dict) -> dict:
    attributes = participant.get("attributes")
    return attributes if isinstance(attributes, dict) else {}


def shared_attributes(p0: dict, p2: dict, leave_out=()) -> str:
    """The attributes both versions give the same value, to find a participant."""
    a0, a2 = attributes_of(p0), attributes_of(p2)
    same = sorted(
        k for k in a0.keys() & a2.keys() if a0[k] == a2[k] and k not in leave_out
    )
    return "; ".join(f"{k}={text(a0[k])}" for k in same)


def shared_analyte(a0: dict, a2: dict) -> str:
    """The locating fields both versions give the same value."""
    return "; ".join(
        f"{f}={text(a0.get(f))}"
        for f in LOCATING_FIELDS
        if a0.get(f) == a2.get(f) and a0.get(f) is not None
    )


def whole_participant(participant: dict) -> str:
    """A participant present in one version only, with a few readings."""
    measurements = align.measurements_of(participant)
    out = [f"{k}={text(v)}" for k, v in sorted(attributes_of(participant).items())]
    out.append(f"{len(measurements)} measurements")
    out += [f"[{reading(m)}]" for m in measurements[:3]]
    return "; ".join(out)


def join(*pieces) -> str:
    return " | ".join(p for p in pieces if p)


class Study:
    """The draft and the reference of one study, paired as the comparison does."""

    def __init__(self, s0: dict, s2: dict, tol: align.Tolerances):
        self.s0, self.s2 = s0, s2
        self.a0, self.a2 = align.analytes_of(s0), align.analytes_of(s2)
        self.p0, self.p2 = align.participants_of(s0), align.participants_of(s2)
        self.pairs = align.match_analytes(self.a0, self.a2)[0]
        self.back = {k2: k0 for k0, k2 in self.pairs.items()}
        # An extra measurement carries only its place in S0, so the partner of
        # its participant comes from here.
        self.people = align.match_participants(
            self.p0, self.p2, self.pairs, self.a2, tol
        )[0]

    def analyte_where(self, k0=None, k2=None) -> str:
        k0 = k0 if k0 is not None else self.back.get(k2)
        k2 = k2 if k2 is not None else self.pairs.get(k0)
        if k0 in self.a0 and k2 in self.a2:
            return "analyte: " + shared_analyte(self.a0[k0], self.a2[k2])
        return "analyte: (not defined in one version)"

    def participant_where(self, i0, i2, leave_out=()) -> str:
        if i0 is None or i2 is None or i0 >= len(self.p0) or i2 >= len(self.p2):
            return ""
        return join(
            f"participant at position {i2 + 1}",
            shared_attributes(self.p0[i0], self.p2[i2], leave_out),
        )


def describe(row: dict, study: Study) -> tuple[str, str, str, str]:
    """(difference, where, value in S0, value in S2) for one discrepancy row."""
    kind, level = row["type"], row["level"]
    s0_path, s2_path = row["s0_path"], row["s2_path"]
    anchor = row["s0_anchor"] or s0_path
    i0, i2 = index_after(anchor, "participants"), index_after(s2_path, "participants")
    if i2 is None and i0 is not None:
        i2 = study.people.get(i0)

    if kind in ("missing analyte", "extra analyte"):
        v0 = text(row["s0_value"]) if kind == "extra analyte" else ABSENT
        v2 = text(row["s2_value"]) if kind == "missing analyte" else ABSENT
        return "an analyte in one version only", "the study's analytes", v0, v2
    if kind == "participant-count mismatch":
        return (
            "number of participants",
            "the study",
            text(row["s0_value"]),
            text(row["s2_value"]),
        )
    if level == "analyte":
        field = kind.removeprefix("wrong ")
        k0 = parts(anchor)[1]
        a0, a2 = study.a0.get(k0, {}), study.a2.get(study.pairs.get(k0), {})
        return (
            f"analyte field {field}",
            study.analyte_where(k0=k0),
            text(a0.get(field)),
            text(a2.get(field)),
        )
    if level == "participant" and kind in ("missing", "extra"):
        v0 = whole_participant(study.p0[i0]) if kind == "extra" else ABSENT
        v2 = whole_participant(study.p2[i2]) if kind == "missing" else ABSENT
        return "a participant in one version only", "the study's participants", v0, v2
    if kind == "wrong attribute value":
        name = parts(s2_path)[-1]
        return (
            f"participant attribute {name}",
            study.participant_where(i0, i2, {name}),
            text(row["s0_value"]),
            text(row["s2_value"]),
        )
    if kind == "attribute in the wrong field":
        names = {parts(s0_path)[-1], parts(s2_path)[-1]}
        return (
            "an attribute recorded under different names",
            study.participant_where(i0, i2, names),
            text(row["s0_value"]),
            text(row["s2_value"]),
        )
    if kind == "missing attribute":
        name = parts(s2_path)[-1]
        return (
            f"participant attribute {name}",
            study.participant_where(i0, i2, {name}),
            NOT_GIVEN,
            text(row["s2_value"]),
        )
    if level == "measurement":
        j0 = index_after(s0_path, "measurements")
        j2 = index_after(s2_path, "measurements")
        m0 = align.measurements_of(study.p0[i0])[j0] if j0 is not None else None
        m2 = align.measurements_of(study.p2[i2])[j2] if j2 is not None else None
        where = join(
            study.participant_where(i0, i2),
            study.analyte_where(
                k0=m0.get("analyte") if m0 else None,
                k2=m2.get("analyte") if m2 and not m0 else None,
            ),
        )
        label = {
            "omission": "a measurement in one version only",
            "extra": "a measurement in one version only",
            "wrong time": "the time of a measurement",
            "wrong value": "the value of a measurement",
            "wrong non-detect status": "whether a measurement detected the analyte",
        }[kind]
        return (
            label,
            where,
            reading(m0) if m0 is not None else ABSENT,
            reading(m2) if m2 is not None else ABSENT,
        )
    raise ValueError(f"no wording for {level} / {kind}")


def finding_where(path: str, s0: dict) -> str:
    """What a finding's path points at in the draft, in words."""
    segments = parts(path)
    i = index_after(path, "participants")
    if segments[:1] == ["analytes"] and len(segments) > 1:
        analyte = align.analytes_of(s0).get(segments[1], {})
        return "analyte: " + "; ".join(
            f"{f}={text(analyte.get(f))}" for f in LOCATING_FIELDS
        )
    if i is not None:
        people = align.participants_of(s0)
        if i < len(people):
            attributes = "; ".join(
                f"{k}={text(v)}" for k, v in sorted(attributes_of(people[i]).items())
            )
            return join(f"draft participant at position {i + 1}", attributes)
    return "the study"


def content(node, limit: int = 300) -> str:
    if isinstance(node, list):
        return f"(a list of {len(node)} items)"
    rendered = text(node) if node is not None else "(nothing at this path)"
    return rendered if len(rendered) <= limit else rendered[: limit - 3] + "..."


# ------------------------------------------------------------------- packets


def build(
    discrepancy_rows: list[dict],
    finding_rows: list[dict],
    findings: dict,
    studies: dict,
    manifest: dict,
    seed: int = SEED,
):
    """The two packets and the key.

    `findings` maps (study, finding id) to the finding as the review file has
    it, and `studies` maps a study to its `Study`; both need only cover the
    studies with rows. Returns (discrepancy packet, finding packet, key).
    """
    era_of = {s: row["era"] for s, row in manifest.items()}
    for row in discrepancy_rows:
        row["_ref_id"] = row["id"]
    for row in finding_rows:
        row["_ref_id"] = f"{row['study_id']}:{row['finding_id']}"

    chosen_d = select(
        discrepancy_rows,
        era_of,
        seed,
        "D",
        lambda r: r["scored"] == "yes" and r["type"] not in NOT_ADJUDICATED,
    )
    chosen_f = select(
        finding_rows, era_of, seed, "F", lambda r: r["adjudicate"] == "yes"
    )

    flips = random.Random(f"{seed}:ab")
    d_packet, f_packet, key = [], [], []
    for number, row in enumerate(chosen_d, 1):
        study = studies[row["study_id"]]
        difference, where, v0, v2 = describe(row, study)
        s0_first = flips.random() < 0.5
        item = {
            "item_id": f"D{number:04d}",
            "study_id": row["study_id"],
            "source": source(manifest[row["study_id"]], study.s2),
            "difference": difference,
            "where": where,
            "A": v0 if s0_first else v2,
            "B": v2 if s0_first else v0,
            "outcome": "",
            "correct_value": "",
            "note": "",
        }
        d_packet.append(item)
        key.append(
            {
                "item_id": item["item_id"],
                "packet": "discrepancies",
                "ref_id": row["id"],
                "study_id": row["study_id"],
                "era": era_of[row["study_id"]],
                "type": row["type"],
                "A_is": "S0" if s0_first else "S2",
                "B_is": "S2" if s0_first else "S0",
            }
        )
    for number, row in enumerate(chosen_f, 1):
        study = studies[row["study_id"]]
        finding = findings[(row["study_id"], row["finding_id"])]
        patch = finding.get("suggested_patch")
        item = {
            "item_id": f"F{number:04d}",
            "study_id": row["study_id"],
            "source": source(manifest[row["study_id"]], study.s2),
            "check": row["check"],
            "severity": row["severity"],
            "confidence": text(finding.get("confidence")),
            "path": row["path"],
            "where": finding_where(row["path"], study.s0),
            "draft_content": content(node_at(study.s0, row["path"])),
            "message": str(finding.get("message", "")),
            "suggested_patch": "" if patch is None else text(patch),
            "outcome": "",
            "note": "",
        }
        f_packet.append(item)
        key.append(
            {
                "item_id": item["item_id"],
                "packet": "findings",
                "ref_id": row["_ref_id"],
                "study_id": row["study_id"],
                "era": era_of[row["study_id"]],
                "type": row["check"],
                "A_is": "",
                "B_is": "",
            }
        )

    for packet, prefix in (("discrepancies", "D"), ("findings", "F")):
        rows = [k for k in key if k["packet"] == packet]
        again = recheck_ids(rows, seed, prefix)
        for k in rows:
            k["recheck"] = "yes" if k["item_id"] in again else "no"
    return d_packet, f_packet, key


def source(manifest_row: dict, s2: dict) -> str:
    doi = manifest_row.get("doi") or (s2.get("doi") if isinstance(s2, dict) else "")
    if doi:
        return f"https://doi.org/{doi}"
    return str(s2.get("url", "")) if isinstance(s2, dict) else ""


# ------------------------------------------------------------------------ CLI


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir", type=Path, help="benchmark/runs/<run_id>")
    ap.add_argument("--comparison", type=Path, required=True, help="Comparison folder")
    ap.add_argument("--out", type=Path, required=True, help="Folder for the packets")
    ap.add_argument("--ref", default=disc.RELEASE, help="Release to read S2 from")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args(argv)

    d_file = args.comparison / "discrepancies.csv"
    f_file = args.comparison / "finding_matches.csv"
    discrepancy_rows, finding_rows = read_csv(d_file), read_csv(f_file)
    manifest = disc.read_manifest()

    studies, findings = {}, {}
    for study_id in sorted({r["study_id"] for r in discrepancy_rows + finding_rows}):
        draft = args.run_dir / study_id / "s0.yaml"
        s0 = (
            yaml.load(draft.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader)
            if draft.is_file()
            else None
        )
        s2 = disc.load_reference(manifest[study_id]["reference_yaml"], args.ref)
        tol = align.tolerances_for(manifest[study_id]["input_csv_origin"])
        studies[study_id] = Study(s0 if isinstance(s0, dict) else {}, s2, tol)
        review = args.run_dir / study_id / "findings.yaml"
        if review.is_file():
            for finding in mf.read_findings(review)[0]:
                findings[(study_id, str(finding.get("id", "")))] = finding

    d_packet, f_packet, key = build(
        discrepancy_rows, finding_rows, findings, studies, manifest, args.seed
    )
    again = {k["item_id"] for k in key if k["recheck"] == "yes"}

    args.out.mkdir(parents=True, exist_ok=True)
    disc.write_csv(args.out / "discrepancy_packet.csv", DISCREPANCY_COLUMNS, d_packet)
    disc.write_csv(args.out / "findings_packet.csv", FINDING_COLUMNS, f_packet)
    disc.write_csv(
        args.out / "recheck_discrepancy_packet.csv",
        DISCREPANCY_COLUMNS,
        [i for i in d_packet if i["item_id"] in again],
    )
    disc.write_csv(
        args.out / "recheck_findings_packet.csv",
        FINDING_COLUMNS,
        [i for i in f_packet if i["item_id"] in again],
    )
    disc.write_csv(args.out / "key.csv", KEY_COLUMNS, key)
    (args.out / "README.txt").write_text(
        INSTRUCTIONS.format(
            run=args.run_dir.name,
            primary=PRIMARY,
            second=SECOND,
            absent=ABSENT,
            not_given=NOT_GIVEN,
        ),
        encoding="utf-8",
    )

    def counts(packet):
        rows = [k for k in key if k["packet"] == packet]
        return {
            era: {
                "items": sum(k["era"] == era for k in rows),
                "recheck": sum(k["era"] == era and k["recheck"] == "yes" for k in rows),
            }
            for era in ("manual", "ai_assisted")
        }

    info = {
        "run": args.run_dir.name,
        "reference": args.ref,
        "seed": args.seed,
        "plan": "benchmark-plan-v1",
        "inputs": {
            "discrepancies.csv": sha256(d_file),
            "finding_matches.csv": sha256(f_file),
        },
        "discrepancies": counts("discrepancies"),
        "findings": counts("findings"),
        "primary_adjudicator": PRIMARY,
        "second_adjudicator": SECOND,
    }
    (args.out / "packet_info.json").write_text(
        json.dumps(info, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"{len(d_packet)} discrepancies and {len(f_packet)} findings to adjudicate, "
        f"{len(again)} to recheck -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
