"""Cost and time are read from run logs and priced as documented."""

import csv
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import cost  # noqa: E402


def log(
    prompt=100_000, completion=10_000, review_in=20_000, review_out=1_000, cached=0
):
    """A log shaped like the runner's, with round numbers."""
    return {
        "started": "2026-10-07T19:47:20+00:00",
        "finished": True,
        "finished_at": "2026-10-07T19:48:50+00:00",
        "extraction": {
            "seconds": 80.0,
            "api_calls": 4,
            "models_returned": ["gpt-5.2-2025-12-11"],
            "usage_totals": {"prompt_tokens": prompt, "completion_tokens": completion},
        },
        "review": {
            "seconds": 9.0,
            "api_calls": 1,
            "models_returned": ["claude-opus-4-8"],
            "usage_totals": {
                "input_tokens": review_in,
                "output_tokens": review_out,
                "cache_creation_input_tokens": 0,
                "cache_read_input_tokens": cached,
            },
        },
    }


def test_tokens_are_priced_per_million():
    row = cost.study_cost(log(), measurements=50)
    assert math.isclose(row["extraction_usd"], 0.1 * 1.75 + 0.01 * 14.0)
    assert math.isclose(row["review_usd"], 0.02 * 5.0 + 0.001 * 25.0)
    assert math.isclose(row["total_usd"], row["extraction_usd"] + row["review_usd"])
    assert row["wall_clock_seconds"] == 90.0
    assert math.isclose(row["usd_per_100_measurements"], 2 * row["total_usd"])
    assert row["seconds_per_100_measurements"] == 180.0
    assert row["models"] == "claude-opus-4-8;gpt-5.2-2025-12-11"


def test_cache_reads_are_cheaper_than_input():
    plain = cost.study_cost(log(review_in=20_000), measurements=50)
    cached = cost.study_cost(log(review_in=0, cached=20_000), measurements=50)
    assert math.isclose(plain["review_usd"] - cached["review_usd"], 0.02 * (5.0 - 0.5))


def test_a_failed_extraction_still_costs_what_it_used():
    broken = log()
    broken["finished"] = False
    broken["review"] = {}
    del broken["finished_at"]
    row = cost.study_cost(broken, measurements=0)
    assert row["review_usd"] == 0 and row["extraction_usd"] > 0
    assert row["wall_clock_seconds"] is None
    assert row["usd_per_100_measurements"] is None


def test_summaries_give_median_quartiles_and_totals():
    rows = [
        {"study_id": f"s{i}", "era": "manual", **cost.study_cost(log(prompt=p), 100)}
        for i, p in enumerate((100_000, 200_000, 300_000, 400_000, 500_000))
    ]
    summary = {r["quantity"]: r for r in cost.summarize("manual", rows)}
    tokens = summary["extraction_input_tokens"]
    assert (tokens["median"], tokens["q1"], tokens["q3"]) == (300_000, 200_000, 400_000)
    assert tokens["total"] == 1_500_000
    assert summary["usd_per_100_measurements"]["total"] == ""


def test_the_command_line_reads_a_run(tmp_path):
    run, cmp_dir, out = tmp_path / "run", tmp_path / "cmp", tmp_path / "out"
    for study in ("tan2021early", "antinori2022epidemiological"):
        (run / study).mkdir(parents=True)
        (run / study / "log.json").write_text(json.dumps(log()), encoding="utf-8")
    (run / "not-a-study").mkdir()
    cmp_dir.mkdir()
    (cmp_dir / "comparison_summary.csv").write_text(
        "study_id,s2_measurements\ntan2021early,82\nantinori2022epidemiological,24\n",
        encoding="utf-8",
    )
    assert cost.main([str(run), "--comparison", str(cmp_dir), "--out", str(out)]) == 0
    rows = list(csv.DictReader((out / "cost_studies.csv").open(encoding="utf-8")))
    assert [(r["study_id"], r["era"]) for r in rows] == [
        ("antinori2022epidemiological", "ai_assisted"),
        ("tan2021early", "manual"),
    ]
    groups = {
        r["group"]
        for r in csv.DictReader((out / "cost_summary.csv").open(encoding="utf-8"))
    }
    assert groups == {"all", "manual", "ai_assisted"}
    assert json.loads((out / "prices.json").read_text(encoding="utf-8")) == cost.PRICES
    # A study folder with no log is reported and fails the run.
    (run / "xu2020characteristics").mkdir()
    assert cost.main([str(run), "--comparison", str(cmp_dir), "--out", str(out)]) == 1
