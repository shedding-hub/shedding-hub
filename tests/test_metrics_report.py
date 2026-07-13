import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import metrics_report as mr  # noqa: E402

SCRIPT = str(Path(__file__).resolve().parent.parent / "scripts" / "metrics_report.py")


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


SAMPLE = [
    {
        "week_start": "2026-06-22",
        "week_end": "2026-06-28",
        "github": {
            "stars": 17,
            "forks": 2,
            "views_this_week": 26,
            "clones_this_week": 29,
        },
        "pypi": {"last_week": 23},
        "ga4": {
            "active_users": 30,
            "new_users": 28,
            "page_views": 42,
            "avg_engagement_seconds": 12.9,
            "traffic_sources": [{"source": "Direct", "sessions": 32}],
            "top_countries": [{"country": "United States", "active_users": 19}],
            "device_types": [{"device": "desktop", "active_users": 28}],
            "page_breakdown": [{"page": "/", "views": 23}],
        },
    },
    {
        "week_start": "2026-06-29",
        "week_end": "2026-07-05",
        "github": {
            "stars": 17,
            "forks": 2,
            "views_this_week": 61,
            "clones_this_week": 59,
        },
        "pypi": {"last_week": 0},
        "ga4": {
            "active_users": 30,
            "new_users": 27,
            "page_views": 42,
            "avg_engagement_seconds": 9.1,
            "traffic_sources": [{"source": "Direct", "sessions": 34}],
            "top_countries": [{"country": "Singapore", "active_users": 16}],
            "device_types": [{"device": "desktop", "active_users": 29}],
            "page_breakdown": [{"page": "/", "views": 14}],
        },
    },
]


def test_build_trends_html_structure():
    html = mr.build_trends_html(SAMPLE)
    assert html.startswith("<!DOCTYPE html>")
    assert "Weekly Trends" in html
    assert "GitHub stars" in html
    assert "PyPI downloads" in html
    assert "Traffic sources" in html
    assert "<svg" in html


def test_build_trends_html_empty():
    html = mr.build_trends_html([])
    assert html.startswith("<!DOCTYPE html>")
    assert "No weekly metrics" in html


def test_build_trends_html_single_record_no_prev():
    html = mr.build_trends_html(SAMPLE[:1])
    assert "<svg" in html  # renders even with no previous week for deltas


def test_build_trends_html_missing_keys_does_not_raise():
    thin = [
        {"week_start": "2026-01-01"},
        {"week_start": "2026-01-08", "github": {"stars": 5}},
    ]
    html = mr.build_trends_html(thin)
    assert "<!DOCTYPE html>" in html


def test_build_trends_html_deterministic_with_fixed_timestamp():
    ts = datetime(2026, 7, 13, 6, 0, 0, tzinfo=timezone.utc)
    a = mr.build_trends_html(SAMPLE, generated_at=ts)
    b = mr.build_trends_html(SAMPLE, generated_at=ts)
    assert a == b
    assert "2026-07-13 06:00 UTC" in a


def test_cli_writes_output_file(tmp_path):
    src = tmp_path / "m.jsonl"
    src.write_text(
        '{"week_start":"2026-06-22","week_end":"2026-06-28","github":{"stars":17}}\n',
        encoding="utf-8",
    )
    out = tmp_path / "r.html"
    subprocess.check_call([sys.executable, SCRIPT, str(src), "-o", str(out)])
    assert out.read_text(encoding="utf-8").startswith("<!DOCTYPE html>")


def test_cli_writes_to_stdout(tmp_path):
    src = tmp_path / "m.jsonl"
    src.write_text(
        '{"week_start":"2026-06-22","github":{"stars":17}}\n', encoding="utf-8"
    )
    res = subprocess.run(
        [sys.executable, SCRIPT, str(src)], capture_output=True, text=True, check=True
    )
    assert "<!DOCTYPE html>" in res.stdout


def test_read_history_from_ref_bad_ref_returns_none():
    assert mr.read_history_from_ref("nonexistent/ref-xyz-123") is None


import weekly_report as wr  # noqa: E402  (scripts dir already on sys.path)


def test_build_email_message_attaches_html():
    msg = wr.build_email_message(
        html="<p>weekly body</p>",
        subject="S",
        sender="a@b.com",
        recipients=["c@d.com"],
        attachments=[("trends.html", "<!DOCTYPE html><p>chart</p>")],
    )
    parts = list(msg.walk())
    dispositions = [p.get("Content-Disposition") for p in parts]
    assert any(d and "attachment" in d and "trends.html" in d for d in dispositions)
    bodies = [
        (p.get_payload(decode=True) or b"").decode()
        for p in parts
        if p.get_content_type() == "text/html"
    ]
    assert any("weekly body" in b for b in bodies)


def test_build_email_message_without_attachments():
    msg = wr.build_email_message(
        html="<p>b</p>", subject="S", sender="a@b.com", recipients=["c@d.com"]
    )
    assert msg["Subject"] == "S"
    assert msg["To"] == "c@d.com"
