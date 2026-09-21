"""
Write the `curation` provenance block into every dataset.

Two of the four fields can be derived, and two cannot. Being explicit about
which is which matters more here than filling every slot:

`method` and `extracted_on` come from git. The first commit that added a
dataset's YAML gives its date, and its subject gives the era -- the agent
batches all say "agent-extracted" or "AI-extracted". Those two signals are
independent and agree on all 146 datasets, and the date ranges do not overlap
(manual ends 2025-10-28, agent-assisted begins 2026-02-05), so either alone
would classify correctly. `--check` reports any disagreement rather than
silently preferring one.

`data_source` is not in the repository. A digitized CSV and a supplementary
spreadsheet look alike on disk, and an emailed one looks like both. It is read
from a CSV of curator knowledge instead (--sources), with `main_text` added to
every dataset because the article text is always read alongside whichever
object carries the numbers. A dataset absent from that CSV gets `[main_text]`,
which is true but incomplete; `--check` lists those so the gap stays visible
rather than passing for a complete record.

`reviewers` is left empty. The obvious proxy, the commit author, is whoever
pushed the branch, which for every agent batch is one person and not the
reviewer. Recording a confident wrong answer for 107 datasets is worse than
recording nothing.

Run via `make curation`.
"""

import argparse
import csv
import pathlib
import re
import subprocess
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = REPO_ROOT / "data"

# The agent batches name themselves in the commit subject. Anchored on that
# rather than a date cutoff so a manual dataset landing today is still manual.
AI_SUBJECT = re.compile(r"agent-extracted|AI-extracted|agent extracted", re.I)

# Secondary signal, used only to cross-check. The eras do not overlap.
ERA_BOUNDARY = "2026-01-01"

VALID_SOURCES = {
    "main_text",
    "table",
    "figure",
    "supplement",
    "author_shared",
    "data_repository",
}
# Kept in enum order so the written arrays read consistently.
SOURCE_ORDER = [
    "main_text",
    "table",
    "figure",
    "supplement",
    "author_shared",
    "data_repository",
]


def first_commit(path: pathlib.Path) -> tuple[str, str]:
    """Date and subject of the commit that introduced this file."""
    out = subprocess.run(
        [
            "git",
            "log",
            "--diff-filter=A",
            "--format=%ad|%s",
            "--date=short",
            "--reverse",
            "--",
            str(path.relative_to(REPO_ROOT)),
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        encoding="utf-8",
    ).stdout
    line = out.splitlines()[0] if out.strip() else ""
    date, _, subject = line.partition("|")
    return date.strip(), subject.strip()


def read_sources(path: pathlib.Path | None) -> dict:
    """study_id -> list of data_source values, from curator knowledge."""
    if path is None:
        return {}
    sources = {}
    with path.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            study = (row.get("study_id") or "").strip()
            raw = (
                row.get("data_source") or row.get("proposed_data_source") or ""
            ).strip()
            if not study:
                continue
            vals = [v.strip() for v in raw.split(";") if v.strip()]
            bad = set(vals) - VALID_SOURCES
            if bad:
                raise SystemExit(f"{study}: unknown data_source {sorted(bad)}")
            sources[study] = vals
    return sources


def render(method, sources, extracted_on) -> str:
    """The YAML block, written textually so the rest of the file is untouched."""
    lines = ["curation:", f"  method: {method}", "  data_source:"]
    lines += [f"    - {v}" for v in sources]
    lines.append(f"  extracted_on: '{extracted_on}'")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sources",
        type=pathlib.Path,
        default=None,
        help="CSV with study_id and data_source (semicolon-separated).",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Report what would change and what is incomplete; write nothing.",
    )
    args = parser.parse_args()

    sources = read_sources(args.sources)
    disagreed, incomplete, written, already, no_main_text = [], [], [], [], []

    for path in sorted(DATA.glob("*/*.yaml")):
        study = path.parent.name
        text = path.read_text(encoding="utf-8")
        if "\ncuration:\n" in text or text.startswith("curation:"):
            already.append(study)
            continue

        date, subject = first_commit(path)
        by_subject = "ai_assisted" if AI_SUBJECT.search(subject) else "manual"
        by_date = "ai_assisted" if date >= ERA_BOUNDARY else "manual"
        if by_subject != by_date:
            disagreed.append((study, date, by_subject, by_date, subject))

        # main_text is injected only where a curator has said nothing at all.
        # An explicit value is the curator's judgement and is kept as written:
        # cdc2024nhphrn is a surveillance network with no paper, so its data
        # come from a repository and no article text exists to have been read.
        # Forcing main_text onto it would record something false.
        if study in sources:
            vals = list(sources[study])
            if "main_text" not in vals:
                no_main_text.append(study)
        else:
            vals = ["main_text"]
        vals = [v for v in SOURCE_ORDER if v in vals]
        if vals == ["main_text"]:
            incomplete.append(study)

        block = render(by_subject, vals, date)
        anchor = "\nanalytes:\n"
        if anchor not in text:
            raise SystemExit(f"{study}: no top-level `analytes:` to insert before")
        if not args.check:
            path.write_text(
                text.replace(anchor, "\n" + block + "analytes:\n", 1), encoding="utf-8"
            )
        written.append(study)

    print(f"datasets written        : {len(written)}")
    print(f"already had a block     : {len(already)}")
    print(
        f"data_source incomplete  : {len(incomplete)} (main_text only, awaiting curator input)"
    )
    if incomplete:
        for s in incomplete:
            print(f"    {s}")
    if no_main_text:
        print(f"no main_text, as the curator recorded: {len(no_main_text)}")
        for s in no_main_text:
            print(f"    {s}")
    if disagreed:
        print(
            f"\nERA SIGNALS DISAGREE on {len(disagreed)} dataset(s) -- resolve before trusting:"
        )
        for row in disagreed:
            print(f"    {row}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
