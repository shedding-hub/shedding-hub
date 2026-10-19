"""Analysis C reads merges, batches, contributors and edits from git correctly.

The tests build a small repository with known dates: two pull requests that
add datasets, a direct edit, a bulk commit, an edit after the window and a
late merge with too little follow-up.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import git_history as gh  # noqa: E402

DAY0 = "2025-01-01T12:00:00-05:00"


def at(days: int) -> str:
    from datetime import datetime, timedelta

    return (datetime.fromisoformat(DAY0) + timedelta(days=days)).isoformat()


class Repo:
    def __init__(self, path: Path):
        self.path = path
        self.run("init", "-q", "-b", "main")

    def run(self, *args, author="Yuke Wang", when=DAY0):
        env = {
            **os.environ,
            "GIT_AUTHOR_NAME": author,
            "GIT_AUTHOR_EMAIL": "a@example.org",
            "GIT_COMMITTER_NAME": author,
            "GIT_COMMITTER_EMAIL": "a@example.org",
            "GIT_AUTHOR_DATE": when,
            "GIT_COMMITTER_DATE": when,
        }
        return subprocess.run(
            ["git", "-c", "commit.gpgsign=false", *args],
            cwd=self.path,
            env=env,
            capture_output=True,
            check=True,
        ).stdout.decode()

    def dataset(self, study, era, participants=2, per_participant=3, note=""):
        folder = self.path / "data" / study
        folder.mkdir(parents=True, exist_ok=True)
        data = {
            "title": study + note,
            "curation": {"method": era},
            "participants": [
                {
                    "measurements": [
                        {"time": t, "value": 1.0} for t in range(per_participant)
                    ]
                }
                for _ in range(participants)
            ],
        }
        (folder / f"{study}.yaml").write_text(yaml.safe_dump(data), encoding="utf-8")

    def commit(self, message, author="Yuke Wang", day=0):
        self.run("add", "-A", author=author, when=at(day))
        self.run("commit", "-q", "-m", message, author=author, when=at(day))

    def merge(self, branch, number, day, author="Yuke Wang"):
        self.run("switch", "-q", "main")
        self.run(
            "merge",
            "-q",
            "--no-ff",
            "-m",
            f"Merge pull request #{number} from org/{branch}",
            branch,
            author=author,
            when=at(day),
        )


@pytest.fixture(scope="module")
def tables(tmp_path_factory):
    repo = Repo(tmp_path_factory.mktemp("history"))
    (repo.path / "README.md").write_text("x", encoding="utf-8")
    repo.commit("start", day=0)

    # PR 1: one manual dataset, drafted by one person and fixed by another.
    repo.run("switch", "-q", "-c", "one")
    repo.dataset("alpha2020one", "manual")
    repo.commit("add alpha", author="Liv Z", day=1)
    repo.dataset("alpha2020one", "manual", per_participant=4)
    repo.commit("fix alpha", author="liv zhai", day=2)
    repo.merge("one", 1, day=3)

    # PR 2: a batch of two AI-assisted datasets.
    repo.run("switch", "-q", "-c", "two")
    repo.dataset("beta2021two", "ai_assisted", participants=1)
    repo.dataset("gamma2022three", "ai_assisted", participants=3)
    repo.commit("add batch", author="hww228", day=10)
    repo.merge("two", 2, day=11)

    # An edit of alpha on main within 90 days.
    repo.dataset("alpha2020one", "manual", per_participant=4, note=" fixed")
    repo.commit("Fix alpha (#3)", day=30)

    # A bulk commit touching all three folders (threshold 2 in the tests).
    for study, era in (
        ("alpha2020one", "manual"),
        ("beta2021two", "ai_assisted"),
        ("gamma2022three", "ai_assisted"),
    ):
        path = repo.path / "data" / study / f"{study}.yaml"
        path.write_text(
            path.read_text(encoding="utf-8") + "# migrated\n", encoding="utf-8"
        )
    repo.commit("Migrate the schema (#4)", day=40)

    # An edit of beta after its window.
    path = repo.path / "data" / "beta2021two" / "beta2021two.yaml"
    path.write_text(path.read_text(encoding="utf-8") + "# late\n", encoding="utf-8")
    repo.commit("Late fix (#5)", day=200)

    # A late merge with less than 90 days of follow-up, edited soon after.
    repo.run("switch", "-q", "-c", "six")
    repo.dataset("delta2023four", "manual", participants=1, per_participant=1)
    repo.commit("add delta", author="Till Hoffmann", day=250)
    repo.merge("six", 6, day=251)
    repo.dataset("delta2023four", "manual", participants=1, per_participant=2)
    repo.commit("Fix delta (#7)", day=260)
    repo.run("tag", "v9", when=at(260))

    out = gh.history(repo.path, "v9", analysis_set={"alpha2020one"}, bulk_threshold=2)
    out["repo"] = str(repo.path)
    return out, {row["study_id"]: row for row in out["datasets"]}


def test_a_dataset_is_merged_by_its_pull_request(tables):
    _, rows = tables
    assert rows["alpha2020one"]["merged_on"] == at(3)[:10]
    assert rows["alpha2020one"]["pull_request"] == "1"
    assert rows["beta2021two"]["pull_request"] == "2"
    assert rows["alpha2020one"]["era"] == "manual"
    assert rows["beta2021two"]["era"] == "ai_assisted"


def test_counts_are_taken_as_merged(tables):
    _, rows = tables
    assert (
        rows["alpha2020one"]["participants"],
        rows["alpha2020one"]["measurements"],
    ) == (2, 8)
    assert (
        rows["delta2023four"]["participants"],
        rows["delta2023four"]["measurements"],
    ) == (1, 1)


def test_batches_and_months(tables):
    out, _ = tables
    sizes = {b["pull_request"]: b["datasets"] for b in out["batches"]}
    assert sizes == {"1": 1, "2": 2, "6": 1}
    two = next(b for b in out["batches"] if b["pull_request"] == "2")
    assert two["measurements"] == 3 + 9 and two["era"] == "ai_assisted"
    assert {(m["month"], m["era"]): m["datasets"] for m in out["monthly"]} == {
        ("2025-01", "manual"): 1,
        ("2025-01", "ai_assisted"): 2,
        ("2025-09", "manual"): 1,
    }


def test_contributors_merge_aliases_and_stop_at_the_merge(tables):
    _, rows = tables
    assert rows["alpha2020one"]["contributors"] == "Liv Zhai"
    assert rows["beta2021two"]["contributors"] == "Youwei Hu"
    # Yuke Wang edits alpha only after the merge.
    assert rows["delta2023four"]["contributors"] == "Till Hoffmann"


def test_edits_skip_bulk_commits_and_stop_at_the_window(tables):
    out, rows = tables
    assert rows["alpha2020one"]["edits_90d"] == 1
    assert rows["alpha2020one"]["bulk_commits_90d"] == 1
    assert rows["beta2021two"]["edits_90d"] == 0  # its edit came on day 200
    assert rows["gamma2022three"]["edits_90d"] == 0
    assert [b["pull_request"] for b in out["bulk"]] == ["4"]
    assert out["bulk"][0]["folders_changed"] == 3
    assert [(e["pull_request"], e["datasets"]) for e in out["edits"]] == [
        ("3", "alpha2020one"),
        ("7", "delta2023four"),
    ]


def test_the_edit_rate_counts_only_full_follow_up(tables):
    out, rows = tables
    assert rows["delta2023four"]["follow_up_days"] == 9
    summary = gh.summarize(out)
    manual = summary["eras"]["manual"]
    assert manual["datasets"] == 2
    assert manual["datasets_with_full_follow_up"] == 1
    assert manual["edited_within_window"] == 1
    assert manual["edit_rate"] == 1.0
    assert manual["datasets_in_analysis_set"] == 1
    ai = summary["eras"]["ai_assisted"]
    assert (ai["batches"], ai["median_batch_size"], ai["edit_rate"]) == (1, 2, 0.0)
    assert summary["datasets_not_found_in_history"] == []


def test_bots_are_not_people():
    assert gh.person("dependabot[bot]") is None
    assert gh.person("github-actions") is None
    assert gh.person("YouweiHu") == "Youwei Hu"
    assert gh.person("Someone New") == "Someone New"


@pytest.fixture(scope="module")
def agents(tmp_path_factory):
    """An agents repository that, like the real one, took in an older history
    by a subtree merge and later moved the files back to the top level."""
    repo = Repo(tmp_path_factory.mktemp("agents"))
    (repo.path / "README.md").write_text("x", encoding="utf-8")
    repo.commit("start", day=0)

    repo.run("switch", "-q", "--orphan", "old")
    (repo.path / "data").mkdir(exist_ok=True)
    (repo.path / "data" / "beta2021two.yaml").write_text("draft", encoding="utf-8")
    repo.commit("draft beta", author="YouweiHu", day=5)
    repo.run("switch", "-q", "main")
    repo.run(
        "merge",
        "-q",
        "-s",
        "ours",
        "--no-commit",
        "--allow-unrelated-histories",
        "old",
        when=at(6),
    )
    repo.run("read-tree", "--prefix=extract/", "-u", "old")
    repo.commit("Add 'extract/' from commit", day=6)
    (repo.path / "data").mkdir(exist_ok=True)
    repo.run("mv", "extract/data/beta2021two.yaml", "data/beta2021two.yaml")
    repo.commit("lift data to the top level", day=7)

    # An accented file name, and an edit after beta's merge on day 11.
    (repo.path / "figure_data").mkdir()
    (repo.path / "figure_data" / "gammá2022three.csv").write_text("t", encoding="utf-8")
    repo.commit("input for gamma", author="Liv Z", day=9)
    (repo.path / "data" / "beta2021two.yaml").write_text("late", encoding="utf-8")
    repo.commit("late change", author="Till Hoffmann", day=20)
    return gh.AgentsRepo(repo.path, "main")


def test_agents_repo_contributors_follow_subtree_merges_and_accents(agents):
    merged = at(11)[:10]
    assert agents.contributors("beta2021two", merged) == ["Youwei Hu", "Yuke Wang"]
    assert agents.contributors("gamma2022three", merged) == ["Liv Zhai"]
    assert agents.contributors("alpha2020one", merged) == []
    assert "Till Hoffmann" in agents.contributors("beta2021two", at(30)[:10])


def test_agents_repo_contributors_count_for_ai_assisted_datasets_only(tables, agents):
    out, _ = tables
    repo = Path(out["repo"])
    again = gh.history(repo, "v9", bulk_threshold=2, agents=agents)
    rows = {row["study_id"]: row for row in again["datasets"]}
    assert rows["beta2021two"]["contributors"] == "Youwei Hu;Yuke Wang"
    assert rows["beta2021two"]["contributors_agents_repo"] == "Youwei Hu;Yuke Wang"
    assert rows["gamma2022three"]["contributors"] == "Liv Zhai;Youwei Hu"
    assert rows["alpha2020one"]["contributors_agents_repo"] == ""
    assert gh.summarize(again)["eras"]["ai_assisted"]["contributors"] == 3
