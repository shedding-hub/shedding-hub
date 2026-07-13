import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import metrics_report as mr  # noqa: E402


def test_load_records_sorts_and_skips_bad_lines(tmp_path):
    f = tmp_path / "m.jsonl"
    f.write_text(
        '{"week_start":"2026-02-01","github":{"stars":10}}\n'
        "\n"
        "not json\n"
        '{"week_start":"2026-01-01","github":{"stars":9}}\n',
        encoding="utf-8",
    )
    recs = mr.load_records(f)
    assert [r["week_start"] for r in recs] == ["2026-01-01", "2026-02-01"]


def test_load_records_missing_file_returns_empty(tmp_path):
    assert mr.load_records(tmp_path / "nope.jsonl") == []


def test_parse_jsonl_from_text():
    recs = mr._parse_jsonl(
        '{"week_start":"2026-01-01"}\n\n{"week_start":"2026-01-08"}\n'
    )
    assert [r["week_start"] for r in recs] == ["2026-01-01", "2026-01-08"]
