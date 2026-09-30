"""
Fill `pipeline_version`, `models` and `reviewers` in every `curation` block.

`backfill_curation.py` wrote the fields git can derive and left these three
empty. They are assignments rather than derivations, confirmed by the
maintainers on 30 Sep 2026, and this script records them by era:

- manual: reviewers only. No pipeline, no models.
- ai_assisted, first added by PR #149 (merged 5 Feb 2026): pipeline `v1`,
  extraction only, no review agent.
- every other ai_assisted dataset: pipeline `v2`, extraction then the review
  agent.

`reviewers` is a phase-level assignment (who reviewed that era's datasets),
not a per-dataset review log. ORCID iDs are written as URLs so they cannot be
mistaken for the GitHub handle used where a reviewer has no ORCID iD.

`models` lists the model names as the code requested them. No dated snapshot
was recorded at the time. The extraction model was not uniform within either
pipeline version: an extraction notebook that set `model_name="gpt-5-mini"`
explicitly was used for two batches, and the evidence is the agents repo's
history:

- v1: the four drafts first committed in c3d9488 (5 Feb 2026), which switched
  the notebook from gpt-4o to gpt-5-mini, plus yen2011detection, re-run in the
  same commit. The other eleven were drafted while the notebook used gpt-4o.
- v2: the 8 Jul 2026 batch (d333572, the same sixteen datasets first added
  here by a8c5bc5), drafted with the notebook at gpt-5-mini, except the two
  re-extracted through the pipeline on 22 Jul. Later batches ran through the
  pipeline, whose default is gpt-5.2.

Existing values are replaced, so rerunning is safe. `--check` writes nothing
and exits 1 if any dataset would change.
"""

import argparse
import pathlib
import re
import subprocess
import sys

import yaml

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = REPO_ROOT / "data"

ORCID_WANG = "https://orcid.org/0000-0002-9615-7859"
ORCID_HOFFMANN = "https://orcid.org/0000-0003-4403-0722"
GITHUB_HU = "hww228"  # no ORCID iD

REVIEWERS = {
    "manual": [ORCID_WANG, ORCID_HOFFMANN],
    "ai_assisted": [ORCID_WANG, GITHUB_HU],
}

# Merge commit of PR #149. The datasets it added, relative to its first
# parent, are pipeline v1.
V1_MERGE = "02eba2a"
V1_COUNT = 16

# Commit that first added the 8 Jul 2026 batch, the first v2 batch.
V2_FIRST_BATCH = "a8c5bc5"
V2_FIRST_BATCH_COUNT = 16

V1_GPT5_MINI = {
    "kirby2016vomiting",
    "teunis2015shedding",
    "tjon2006high",
    "tu2008norovirus",
    "yen2011detection",
}
# First-batch datasets re-extracted through the pipeline (gpt-5.2) on 22 Jul.
V2_REEXTRACTED = {"piralla2013different", "pereira2022standardization"}

REVIEW = "claude-opus-4-8 (review)"


def added_by(commit: str) -> set[str]:
    """Dataset ids whose YAML the commit added, relative to its first parent."""
    out = subprocess.run(
        [
            "git",
            "diff",
            "--name-only",
            "--diff-filter=A",
            f"{commit}^1",
            commit,
            "--",
            "data/*/*.yaml",
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        check=True,
    ).stdout
    ids = set()
    for line in out.splitlines():
        parts = line.strip().split("/")
        if len(parts) == 3 and parts[2] == parts[1] + ".yaml":
            ids.add(parts[1])
    return ids


def assign(study: str, method: str, v1: set, v2_mini: set) -> dict:
    """The three field values for one dataset."""
    if method == "manual":
        return {"reviewers": REVIEWERS["manual"]}
    if study in v1:
        extraction = "gpt-5-mini" if study in V1_GPT5_MINI else "gpt-4o"
        return {
            "pipeline_version": "v1",
            "models": [f"{extraction} (extraction)"],
            "reviewers": REVIEWERS["ai_assisted"],
        }
    extraction = "gpt-5-mini" if study in v2_mini else "gpt-5.2"
    return {
        "pipeline_version": "v2",
        "models": [f"{extraction} (extraction)", REVIEW],
        "reviewers": REVIEWERS["ai_assisted"],
    }


def render(values: dict) -> str:
    """The lines for the three fields, in schema order."""
    lines = []
    if "pipeline_version" in values:
        lines.append(f"  pipeline_version: {values['pipeline_version']}")
    for key in ("models", "reviewers"):
        if key in values:
            lines.append(f"  {key}:")
            lines += [f"    - {v}" for v in values[key]]
    return "\n".join(lines) + "\n"


# The curation block ends at the next top-level key. The three fields and
# their list items are stripped from it before the new lines are inserted
# after `extracted_on`.
BLOCK = re.compile(r"^curation:\n(?:(?:  |    ).*\n)*", re.M)
OWNED = re.compile(
    r"^  (?:pipeline_version|models|reviewers):.*\n(?:    - .*\n)*", re.M
)


def rewrite(text: str, values: dict, study: str) -> str:
    match = BLOCK.search(text)
    if not match:
        raise SystemExit(f"{study}: no curation block")
    block = OWNED.sub("", match.group(0))
    anchor = re.search(r"^  extracted_on: .*\n", block, re.M)
    if not anchor:
        raise SystemExit(f"{study}: curation block has no extracted_on")
    block = block[: anchor.end()] + render(values) + block[anchor.end() :]
    return text[: match.start()] + block + text[match.end() :]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Report what would change; write nothing.",
    )
    args = parser.parse_args()

    v1 = added_by(V1_MERGE)
    if len(v1) != V1_COUNT:
        raise SystemExit(f"PR #149 added {len(v1)} datasets, expected {V1_COUNT}")
    first_batch = added_by(V2_FIRST_BATCH)
    if len(first_batch) != V2_FIRST_BATCH_COUNT:
        raise SystemExit(
            f"{V2_FIRST_BATCH} added {len(first_batch)} datasets, "
            f"expected {V2_FIRST_BATCH_COUNT}"
        )
    for name, subset, parent in [
        ("V1_GPT5_MINI", V1_GPT5_MINI, v1),
        ("V2_REEXTRACTED", V2_REEXTRACTED, first_batch),
    ]:
        if not subset <= parent:
            raise SystemExit(f"{name} not within its batch: {sorted(subset - parent)}")
    v2_mini = first_batch - V2_REEXTRACTED

    changed, tally = [], {}
    for path in sorted(DATA.glob("*/*.yaml")):
        study = path.parent.name
        if path.name != f"{study}.yaml":
            continue
        text = path.read_text(encoding="utf-8")
        method = yaml.safe_load(text)["curation"]["method"]
        values = assign(study, method, v1, v2_mini)
        new = rewrite(text, values, study)

        parsed = yaml.safe_load(new)["curation"]
        for key, value in values.items():
            if parsed.get(key) != value:
                raise SystemExit(f"{study}: {key} did not round-trip")

        key = (
            method,
            values.get("pipeline_version", "-"),
            tuple(values.get("models", [])),
        )
        tally[key] = tally.get(key, 0) + 1
        if new != text:
            changed.append(study)
            if not args.check:
                path.write_text(new, encoding="utf-8")

    for (method, version, models), n in sorted(tally.items()):
        print(f"{n:4d}  {method:12s} {version:3s} {', '.join(models) or '-'}")
    verb = "would change" if args.check else "changed"
    print(f"datasets {verb}: {len(changed)}")
    return 1 if args.check and changed else 0


if __name__ == "__main__":
    sys.exit(main())
