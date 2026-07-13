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


def test_parse_jsonl_skips_non_object_json():
    recs = mr._parse_jsonl('5\n"hello"\nnull\n[1,2]\n{"week_start":"2026-01-01"}\n')
    assert recs == [{"week_start": "2026-01-01"}]


def test_sparkline_returns_svg():
    svg = mr.sparkline([1, 3, 2, 5])
    assert svg.startswith("<svg") and "polyline" in svg


def test_sparkline_empty_is_safe():
    assert "<svg" in mr.sparkline([])
    assert "<svg" in mr.sparkline([None, None])


def test_bar_chart_renders_one_rect_per_item():
    svg = mr.bar_chart([("Direct", 30), ("Search", 5)])
    assert svg.count("<rect") == 2
    assert "Direct" in svg


def test_bar_chart_empty_is_safe():
    assert "No data" in mr.bar_chart([])


def test_line_chart_breaks_line_on_none():
    series = [
        {
            "label": "A",
            "color": "#000000",
            "points": [("2026-01-01", 1), ("2026-01-08", None), ("2026-01-15", 3)],
        }
    ]
    svg = mr.line_chart(series)
    assert svg.count("<polyline") == 2  # the gap splits the line into two segments


def test_line_chart_empty_is_safe():
    assert "No data" in mr.line_chart(
        [{"label": "A", "color": "#000000", "points": []}]
    )
