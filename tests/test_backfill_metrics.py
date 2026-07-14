import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import backfill_metrics as bf  # noqa: E402


def test_read_github_csv_parses_and_skips_blank_rows(tmp_path):
    csv_path = tmp_path / "gh.csv"
    csv_path.write_text(
        "week_start,week_end,stars,forks,open_issues,open_prs,views_this_week,"
        "unique_visitors_this_week,clones_this_week,unique_cloners_this_week,"
        "new_datasets_count\n"
        # a fully-filled row
        "2026-04-06,2026-04-12,16,2,7,3,40,8,30,20,1\n"
        # a template row left blank -> skipped
        "2026-04-13,2026-04-19,,,,,,,,,\n"
        # a row with only some cells filled -> blanks become 0
        "2026-04-20,2026-04-26,16,,,,,,,,\n",
        encoding="utf-8",
    )
    rows = bf.read_github_csv(csv_path)
    assert [r["week_start"] for r in rows] == ["2026-04-06", "2026-04-20"]
    assert rows[0]["github"]["views_this_week"] == 40
    assert rows[0]["github"]["new_datasets_count"] == 1
    # blank numeric cells default to 0
    assert rows[1]["github"]["stars"] == 16
    assert rows[1]["github"]["forks"] == 0


def test_pypi_weekly_sums_correct_windows():
    # one download per day for a long stretch -> week=7, month=30
    daily = {}
    d = __import__("datetime").date(2026, 4, 1)
    for _ in range(60):
        daily[d.isoformat()] = 1
        d += __import__("datetime").timedelta(days=1)
    last_week, last_month = bf.pypi_weekly(daily, "2026-05-04", "2026-05-10")
    assert last_week == 7
    assert last_month == 30


def test_pypi_weekly_handles_missing_days():
    daily = {"2026-05-04": 3, "2026-05-07": 2}  # other days absent -> 0
    last_week, last_month = bf.pypi_weekly(daily, "2026-05-04", "2026-05-10")
    assert last_week == 5
    assert last_month == 5  # 30-day window still only sees these two days


def test_assemble_record_matches_schema():
    gh = {k: i for i, k in enumerate(bf.GH_INT_FIELDS)}
    rec = bf.assemble_record(
        "2026-04-06",
        "2026-04-12",
        gh,
        bf._empty_ga4(),
        {"last_week": 5, "last_month": 20},
        "2026-07-14T00:00:00Z",
    )
    assert set(rec) == {
        "week_start",
        "week_end",
        "collected_at",
        "github",
        "pypi",
        "ga4",
    }
    assert set(rec["github"]) == set(bf.GH_INT_FIELDS)
    assert rec["pypi"] == {"last_week": 5, "last_month": 20}
    assert set(rec["ga4"]) >= {"active_users", "page_views", "traffic_sources"}


def test_merge_records_adds_missing_and_preserves_existing(tmp_path):
    jsonl = tmp_path / "m.jsonl"
    jsonl.write_text(
        json.dumps({"week_start": "2026-06-22", "github": {"stars": 17}}) + "\n",
        encoding="utf-8",
    )
    new = [
        {"week_start": "2026-04-06", "github": {"stars": 16}},
        {"week_start": "2026-06-22", "github": {"stars": 999}},  # existing -> ignored
    ]
    merged, added = bf.merge_records(jsonl, new)
    assert added == ["2026-04-06"]
    assert [r["week_start"] for r in merged] == ["2026-04-06", "2026-06-22"]
    # the existing 2026-06-22 record is preserved, not overwritten
    existing = next(r for r in merged if r["week_start"] == "2026-06-22")
    assert existing["github"]["stars"] == 17


def test_merge_records_missing_file_returns_new_only(tmp_path):
    merged, added = bf.merge_records(
        tmp_path / "nope.jsonl", [{"week_start": "2026-04-06"}]
    )
    assert added == ["2026-04-06"]
    assert len(merged) == 1


def test_load_credentials_info_accepts_raw_and_path(tmp_path):
    info = {"type": "service_account", "project_id": "x"}
    assert bf._load_credentials_info(json.dumps(info)) == info
    p = tmp_path / "sa.json"
    p.write_text(json.dumps(info), encoding="utf-8")
    assert bf._load_credentials_info(str(p)) == info
