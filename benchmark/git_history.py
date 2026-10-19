"""Curation throughput from the repository's git history (Analysis C).

Everything is read from the first-parent history of `main` up to the frozen
release, so a pull request counts once, at the moment it was merged.

- A dataset is merged by the first commit on `main` that adds
  `data/<study>/<study>.yaml`. For a pull request that is the merge (or
  squash) commit, so the date is the merge date.
- A batch is one such commit. Its size is the number of datasets it adds.
- Participants and measurements are counted in the file as merged.
- The era is `curation.method` in the release.
- The contributors of a dataset are the people who authored a commit that
  touched its folder, up to and including the merge, with the git identities
  of one person merged (`ALIASES`). An AI-assisted dataset was drafted and
  reviewed in the agents repository and published here, so when that
  repository is given (`--agents-repo`), the people who committed a file
  whose path names the dataset there, on or before the merge date, count
  too. Paths are compared with accents removed.
- A post-merge edit is a later commit on `main`, within 90 days of the merge,
  that changes the dataset's YAML file and touches at most 10 dataset
  folders. A commit that touches more is a bulk commit (a schema migration,
  a label harmonization, a backfill); bulk commits are listed apart. The edit
  rate counts only datasets merged at least 90 days before the release.
  Every commit counted as an edit is listed, so that a reader can see what
  the threshold lets through.

Dates are the committer's local dates, as git records them.

The numbers describe throughput, not person-time: review done outside the
repository leaves no trace in git.

Usage (from the repository root):
    python benchmark/git_history.py --out OUT_DIR [--ref v1.1.0]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import subprocess
import sys
import unicodedata
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "benchmark" / "manifest.csv"
RELEASE = "v1.1.0"
BULK_THRESHOLD = 10  # dataset folders; a commit touching more is a bulk commit
WINDOW_DAYS = 90

# Git identities of the same person, by author name. Names not listed stand
# for themselves.
ALIASES = {
    "Youwei Hu": "Youwei Hu",
    "YouweiHu": "Youwei Hu",
    "Youwei": "Youwei Hu",
    "hww228": "Youwei Hu",
    "WenfeiXiao": "Wenfei Xiao",
    "Wenfei Xiao": "Wenfei Xiao",
    "Liv Z": "Liv Zhai",
    "liv zhai": "Liv Zhai",
    "lllllliv": "Liv Zhai",
}
BOT = re.compile(r"\[bot\]$|^github-actions$|^dependabot", re.I)
PULL_REQUEST = re.compile(r"#(\d+)")
ERAS = ("manual", "ai_assisted")

DATASET_COLUMNS = [
    "study_id",
    "era",
    "in_analysis_set",
    "merge_commit",
    "merged_on",
    "pull_request",
    "batch_size",
    "participants",
    "measurements",
    "contributors",
    "n_contributors",
    "contributors_agents_repo",
    "follow_up_days",
    "edits_90d",
    "bulk_commits_90d",
]
BATCH_COLUMNS = [
    "merge_commit",
    "merged_on",
    "pull_request",
    "era",
    "datasets",
    "measurements",
    "contributors",
]
MONTH_COLUMNS = ["month", "era", "datasets", "measurements"]
BULK_COLUMNS = [
    "commit",
    "date",
    "pull_request",
    "dataset_folders",
    "folders_changed",
    "folders_added",
    "subject",
]
EDIT_COLUMNS = [
    "commit",
    "date",
    "pull_request",
    "dataset_folders",
    "datasets_edited",
    "datasets",
    "subject",
]
IDENTITY_COLUMNS = ["git_name", "person", "commits_touching_data"]


def git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=repo, capture_output=True, check=True
    ).stdout.decode("utf-8")


def person(name: str) -> str | None:
    """The person behind a git author name, or None for a bot."""
    if BOT.search(name):
        return None
    return ALIASES.get(name, name)


@dataclass
class Commit:
    sha: str
    date: datetime
    subject: str
    changes: list[tuple[str, str]] = field(default_factory=list)  # (status, path)

    @property
    def pull_request(self) -> str:
        found = PULL_REQUEST.search(self.subject)
        return found.group(1) if found else ""

    def folders(self, statuses: str = "AMD") -> set[str]:
        """The dataset folders in which this commit adds, changes or deletes files.

        `statuses` narrows it to some of these.
        """
        out = set()
        for status, path in self.changes:
            parts = path.split("/")
            if (
                len(parts) >= 3
                and parts[0] == "data"
                and not parts[1].startswith(".")
                and status[0] in statuses
            ):
                out.add(parts[1])
        return out

    def changed(self) -> set[str]:
        """The dataset folders in which this commit changes or deletes a file."""
        return self.folders("MD")

    def touches(self, path: str, statuses: str) -> bool:
        return any(p == path and s[0] in statuses for s, p in self.changes)


def main_line(repo: Path, ref: str) -> list[Commit]:
    """The first-parent commits up to `ref`, oldest first, with their changes.

    A merge commit's changes are taken against its first parent, which is
    what the pull request brought in.
    """
    log = git(
        repo,
        "log",
        "--first-parent",
        "--reverse",
        "--diff-merges=first-parent",
        "--no-renames",
        "--name-status",
        "--format=%x00%H%x1f%cI%x1f%s",
        ref,
    )
    commits = []
    for chunk in log.split("\x00")[1:]:
        header, *lines = chunk.strip("\n").split("\n")
        sha, date, subject = header.split("\x1f")
        commit = Commit(sha, datetime.fromisoformat(date), subject)
        for line in lines:
            if "\t" in line:
                status, path = line.split("\t", 1)
                commit.changes.append((status, path))
        commits.append(commit)
    return commits


def read_yaml(repo: Path, sha: str, path: str) -> dict:
    text = git(repo, "show", f"{sha}:{path}")
    return yaml.load(text, Loader=yaml.CSafeLoader) or {}


def released_datasets(repo: Path, ref: str) -> dict:
    """{study: era} for every dataset in the release."""
    out = {}
    for line in git(repo, "ls-tree", "--name-only", ref, "data/").split():
        study = line.split("/")[-1]
        if study.startswith("."):
            continue
        try:
            data = read_yaml(repo, ref, f"data/{study}/{study}.yaml")
        except subprocess.CalledProcessError:
            continue
        out[study] = (data.get("curation") or {}).get("method", "")
    return out


def count(dataset: dict) -> tuple[int, int]:
    participants = dataset.get("participants") or []
    return len(participants), sum(
        len(p.get("measurements") or []) for p in participants
    )


def contributors(repo: Path, sha: str, study: str) -> list[str]:
    names = git(repo, "log", "--format=%an", sha, "--", f"data/{study}/").split("\n")
    return sorted({p for p in map(person, filter(None, names)) if p})


def unaccented(text: str) -> str:
    decomposed = unicodedata.normalize("NFKD", text)
    return "".join(c for c in decomposed if not unicodedata.combining(c))


class AgentsRepo:
    """Who committed each dataset's files in the agents repository."""

    def __init__(self, path: Path, ref: str = "main"):
        self.path, self.ref = path, ref
        names = git(
            path,
            "-c",
            "core.quotepath=false",
            "log",
            ref,
            # The repository took in older histories by subtree merges, under
            # other folders, and later moved the files. Without this, a move
            # lists only the new path and the older commits are not found.
            "--no-renames",
            "--name-only",
            "--format=",
        )
        self.paths = sorted({n for n in names.split("\n") if n})
        self.commit = git(path, "rev-parse", f"{ref}^{{commit}}").strip()

    def files_of(self, study: str) -> list[str]:
        """Every path, past or present, with the dataset's name in it."""
        return [p for p in self.paths if study in unaccented(p)]

    def contributors(self, study: str, until: str) -> list[str]:
        files = self.files_of(study)
        if not files:
            return []
        names = git(
            self.path,
            "log",
            self.ref,
            f"--until={until}T23:59:59",
            "--format=%an",
            "--",
            *(f":(literal){f}" for f in files),
        ).split("\n")
        return sorted({p for p in map(person, filter(None, names)) if p})


def history(
    repo: Path,
    ref: str = RELEASE,
    analysis_set: set | None = None,
    bulk_threshold: int = BULK_THRESHOLD,
    window_days: int = WINDOW_DAYS,
    agents: AgentsRepo | None = None,
) -> dict:
    """The tables of Analysis C, as lists of rows."""
    commits = main_line(repo, ref)
    eras = released_datasets(repo, ref)
    end = commits[-1].date
    window = timedelta(days=window_days)
    bulk = [c for c in commits if len(c.folders()) > bulk_threshold]
    bulk_shas = {c.sha for c in bulk}

    datasets, batches = [], []
    seen = set()
    edits: dict = {}  # sha -> (commit, datasets it counts as an edit of)
    for k, commit in enumerate(commits):
        added = sorted(
            s
            for s in eras
            if s not in seen and commit.touches(f"data/{s}/{s}.yaml", "A")
        )
        if not added:
            continue
        seen.update(added)
        rows = []
        for study in added:
            path = f"data/{study}/{study}.yaml"
            participants, measurements = count(read_yaml(repo, commit.sha, path))
            later = [
                c
                for c in commits[k + 1 :]
                if c.date - commit.date <= window and c.touches(path, "M")
            ]
            for c in later:
                if c.sha not in bulk_shas:
                    edits.setdefault(c.sha, (c, set()))[1].add(study)
            merged_on = commit.date.date().isoformat()
            elsewhere = (
                agents.contributors(study, merged_on)
                if agents is not None and eras[study] == "ai_assisted"
                else []
            )
            people = sorted(set(contributors(repo, commit.sha, study)) | set(elsewhere))
            rows.append(
                {
                    "study_id": study,
                    "era": eras[study],
                    "in_analysis_set": (
                        ""
                        if analysis_set is None
                        else "yes" if study in analysis_set else "no"
                    ),
                    "merge_commit": commit.sha[:7],
                    "merged_on": merged_on,
                    "pull_request": commit.pull_request,
                    "batch_size": len(added),
                    "participants": participants,
                    "measurements": measurements,
                    "contributors": ";".join(people),
                    "n_contributors": len(people),
                    "contributors_agents_repo": ";".join(elsewhere),
                    "follow_up_days": (end - commit.date).days,
                    "edits_90d": sum(c.sha not in bulk_shas for c in later),
                    "bulk_commits_90d": sum(c.sha in bulk_shas for c in later),
                }
            )
        datasets += rows
        batches.append(
            {
                "merge_commit": commit.sha[:7],
                "merged_on": rows[0]["merged_on"],
                "pull_request": commit.pull_request,
                "era": ";".join(sorted({r["era"] for r in rows})),
                "datasets": len(rows),
                "measurements": sum(r["measurements"] for r in rows),
                "contributors": ";".join(
                    sorted({p for r in rows for p in r["contributors"].split(";") if p})
                ),
            }
        )

    months: dict = {}
    for row in datasets:
        key = (row["merged_on"][:7], row["era"])
        entry = months.setdefault(key, {"datasets": 0, "measurements": 0})
        entry["datasets"] += 1
        entry["measurements"] += row["measurements"]
    monthly = [
        {"month": month, "era": era, **entry}
        for (month, era), entry in sorted(months.items())
    ]

    identities: dict = {}
    for line in git(repo, "log", "--format=%an", ref, "--", "data/").split("\n"):
        if line:
            identities[line] = identities.get(line, 0) + 1
    return {
        "datasets": datasets,
        "batches": batches,
        "monthly": monthly,
        "bulk": [
            {
                "commit": c.sha[:7],
                "date": c.date.date().isoformat(),
                "pull_request": c.pull_request,
                "dataset_folders": len(c.folders()),
                "folders_changed": len(c.changed()),
                "folders_added": len(c.folders() - c.changed()),
                "subject": c.subject,
            }
            for c in bulk
        ],
        "edits": [
            {
                "commit": c.sha[:7],
                "date": c.date.date().isoformat(),
                "pull_request": c.pull_request,
                "dataset_folders": len(c.folders()),
                "datasets_edited": len(studies),
                "datasets": ";".join(sorted(studies)),
                "subject": c.subject,
            }
            for c, studies in sorted(edits.values(), key=lambda e: e[0].date)
        ],
        "identities": [
            {"git_name": n, "person": person(n) or "(bot)", "commits_touching_data": k}
            for n, k in sorted(identities.items())
        ],
        "not_found": sorted(set(eras) - seen),
        "end": end.date().isoformat(),
        "window_days": window_days,
    }


def summarize(tables: dict) -> dict:
    """Per era: datasets, measurements, batches, contributors and edit rate."""
    out = {"history_ends": tables["end"], "eras": {}}
    for era in ERAS:
        rows = [r for r in tables["datasets"] if r["era"] == era]
        batches = [b for b in tables["batches"] if b["era"] == era]
        followed = [r for r in rows if r["follow_up_days"] >= tables["window_days"]]
        out["eras"][era] = {
            "datasets": len(rows),
            "measurements": sum(r["measurements"] for r in rows),
            "datasets_in_analysis_set": sum(
                r["in_analysis_set"] == "yes" for r in rows
            ),
            "measurements_in_analysis_set": sum(
                r["measurements"] for r in rows if r["in_analysis_set"] == "yes"
            ),
            "first_merge": min((r["merged_on"] for r in rows), default=""),
            "last_merge": max((r["merged_on"] for r in rows), default=""),
            "batches": len(batches),
            "median_batch_size": (
                statistics.median(b["datasets"] for b in batches) if batches else None
            ),
            "largest_batch": max((b["datasets"] for b in batches), default=0),
            "contributors": len(
                {p for r in rows for p in r["contributors"].split(";") if p}
            ),
            "datasets_with_full_follow_up": len(followed),
            "edited_within_window": sum(r["edits_90d"] > 0 for r in followed),
            "edit_rate": (
                round(sum(r["edits_90d"] > 0 for r in followed) / len(followed), 3)
                if followed
                else None
            ),
        }
    out["mixed_batches"] = [
        b["merge_commit"] for b in tables["batches"] if ";" in b["era"]
    ]
    out["bulk_commits"] = len(tables["bulk"])
    out["edit_commits"] = len(tables["edits"])
    out["datasets_not_found_in_history"] = tables["not_found"]
    return out


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, required=True, help="Folder for the tables")
    ap.add_argument("--ref", default=RELEASE, help="Release whose history to read")
    ap.add_argument(
        "--agents-repo", type=Path, help="Local copy of the agents repository"
    )
    ap.add_argument("--agents-ref", default="main", help="Branch to read there")
    args = ap.parse_args(argv)

    with MANIFEST.open(encoding="utf-8", newline="") as fh:
        analysis_set = {row["study_id"] for row in csv.DictReader(fh)}
    agents = AgentsRepo(args.agents_repo, args.agents_ref) if args.agents_repo else None
    tables = history(REPO_ROOT, args.ref, analysis_set, agents=agents)
    summary = summarize(tables)
    summary["ref"] = args.ref
    summary["commit"] = git(REPO_ROOT, "rev-parse", f"{args.ref}^{{commit}}").strip()
    summary["agents_repo"] = (
        {"ref": agents.ref, "commit": agents.commit} if agents else "not read"
    )

    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "datasets.csv", DATASET_COLUMNS, tables["datasets"])
    write_csv(args.out / "batches.csv", BATCH_COLUMNS, tables["batches"])
    write_csv(args.out / "monthly.csv", MONTH_COLUMNS, tables["monthly"])
    write_csv(args.out / "bulk_commits.csv", BULK_COLUMNS, tables["bulk"])
    write_csv(args.out / "edit_commits.csv", EDIT_COLUMNS, tables["edits"])
    write_csv(args.out / "identities.csv", IDENTITY_COLUMNS, tables["identities"])
    (args.out / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"{len(tables['datasets'])} datasets in {len(tables['batches'])} batches, "
        f"{len(tables['bulk'])} bulk commits -> {args.out}"
    )
    if tables["not_found"]:
        print(f"NOT FOUND {', '.join(tables['not_found'])}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
