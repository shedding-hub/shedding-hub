# Weekly Trends Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone, self-contained HTML trends dashboard from the accumulating `weekly_metrics.jsonl` history and attach it to the existing weekly email.

**Architecture:** A new dependency-free module `scripts/metrics_report.py` turns the JSONL history into an HTML document with hand-built inline SVG charts (`build_trends_html`). `scripts/weekly_report.py` imports it and attaches the output to the email it already sends; the email body is unchanged. The module also has a CLI so the report can be regenerated locally from the `metrics-data` branch.

**Tech Stack:** Python 3.11 standard library only for the report module (`json`, `datetime`, `argparse`, `subprocess`, `pathlib`, `email`). pytest for tests. No new third-party dependencies.

## Global Constraints

- **Python 3.11**; report module uses **standard library only** — do not add any new dependency.
- **Source of truth is the `metrics-data` branch.** History-reading code must read the `metrics-data` version, never `main`'s stale seed. The workflow already restores the full history into `metrics/weekly_metrics.jsonl` before the script runs; the local CLI reads it via `git show origin/metrics-data:metrics/weekly_metrics.jsonl`.
- **Formatting:** all Python must pass `black --check .` (CI enforces this). Run `black` on changed files before every commit.
- **Email body is unchanged.** Only add an attachment; never alter the existing single-week HTML body.
- **Best-effort attachment:** if report generation fails, the email must still send without it.
- **Determinism:** `build_trends_html` output must be identical for identical input except for the injected `generated_at` timestamp.
- **Brand color** matches the email header teal `#1a6b8a`.

---

## File Structure

- **Create `scripts/metrics_report.py`** — report builder: JSONL loading, SVG chart primitives, HTML assembly, CLI. Single responsibility: turn records → HTML. No email or network-collection code.
- **Modify `scripts/weekly_report.py`** — import the report builder, attach its output to the email, make the `anthropic` import lazy so the module is importable in tests.
- **Create `tests/test_metrics_report.py`** — unit tests for loading, charts, HTML assembly, CLI, and the email-attachment builder.
- **Modify `metrics/README.md`** — document the attached report and local regeneration.

---

## Task 1: JSONL loading

**Files:**
- Create: `scripts/metrics_report.py`
- Test: `tests/test_metrics_report.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `load_records(path) -> list[dict]` — parse a JSONL file, skip blank/malformed lines, missing file → `[]`, sorted ascending by `week_start`.
  - `_parse_jsonl(text: str) -> list[dict]` — same parsing from an in-memory string (used later by the CLI's `--ref` path).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_metrics_report.py`:

```python
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
    recs = mr._parse_jsonl('{"week_start":"2026-01-01"}\n\n{"week_start":"2026-01-08"}\n')
    assert [r["week_start"] for r in recs] == ["2026-01-01", "2026-01-08"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_metrics_report.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'metrics_report'`.

- [ ] **Step 3: Write minimal implementation**

Create `scripts/metrics_report.py`:

```python
"""Standalone weekly trends report builder for Shedding Hub.

Reads the accumulating metrics history (``metrics/weekly_metrics.jsonl``, whose
source of truth is the ``metrics-data`` branch) and renders a self-contained HTML
dashboard with inline SVG charts.

Used two ways:
  * imported by ``scripts/weekly_report.py`` to attach the report to the email;
  * run directly to regenerate the report locally, e.g.::

        python scripts/metrics_report.py --ref origin/metrics-data -o report.html
        python scripts/metrics_report.py path/to/weekly_metrics.jsonl -o report.html
"""

import json
from pathlib import Path


def _parse_jsonl(text):
    """Parse JSONL text into records sorted by ``week_start``.

    Blank lines and lines that fail to parse as JSON are skipped.
    """
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    records.sort(key=lambda r: r.get("week_start", ""))
    return records


def load_records(path):
    """Load records from a JSONL file. Missing file returns ``[]``."""
    p = Path(path)
    if not p.exists():
        return []
    return _parse_jsonl(p.read_text(encoding="utf-8"))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_metrics_report.py -v`
Expected: PASS (3 passed).

- [ ] **Step 5: Format and commit**

```bash
black scripts/metrics_report.py tests/test_metrics_report.py
git add scripts/metrics_report.py tests/test_metrics_report.py
git commit -m "feat: add weekly-metrics JSONL loader for trends report"
```

---

## Task 2: SVG chart primitives

**Files:**
- Modify: `scripts/metrics_report.py`
- Test: `tests/test_metrics_report.py`

**Interfaces:**
- Consumes: nothing from other modules.
- Produces (all pure, return SVG/HTML strings):
  - `sparkline(values, *, width=120, height=28, color=INK) -> str` — compact trend line; ignores `None`; empty-safe.
  - `bar_chart(items, *, width=340, bar_color=INK, value_suffix="") -> str` — horizontal bars from `list[(label, value)]`; empty → "No data" paragraph.
  - `line_chart(series, *, width=680, height=240) -> str` — one or more date-indexed series; each series is `{"label": str, "color": "#hex", "points": [(date_str, value_or_None), ...]}`; breaks the line on `None`; empty → "No data" note.
  - Helpers: `_esc(s)`, `_fmt_num(v)`, `_polyline(points, color)`, and palette constants `INK`, `ACCENT`, `MUTED`, `SERIES`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_metrics_report.py`:

```python
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
    assert "No data" in mr.line_chart([{"label": "A", "color": "#000000", "points": []}])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_metrics_report.py -k "sparkline or bar_chart or line_chart" -v`
Expected: FAIL — `AttributeError: module 'metrics_report' has no attribute 'sparkline'`.

- [ ] **Step 3: Write minimal implementation**

Add to `scripts/metrics_report.py` (after the imports, before `_parse_jsonl` add the palette; add the rest below `load_records`):

```python
from datetime import date  # add to the import block at the top of the file

# --- palette (validated against the dataviz skill's contrast rules) ----------
INK = "#1a6b8a"      # brand teal (matches the email header)
ACCENT = "#e07b39"   # warm secondary
MUTED = "#6b7c85"    # gridlines / axis text
SERIES = ["#1a6b8a", "#e07b39", "#4c9a6a", "#9a4c8a"]


def _esc(s):
    return (
        str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )


def _fmt_num(v):
    v = round(v)
    if v >= 1000:
        return f"{v / 1000:.1f}k".replace(".0k", "k")
    return str(int(v))


def _polyline(points, color):
    return (
        f'<polyline fill="none" stroke="{color}" stroke-width="2" '
        f'stroke-linejoin="round" stroke-linecap="round" '
        f'points="{" ".join(points)}"/>'
    )


def sparkline(values, *, width=120, height=28, color=INK):
    vals = [v for v in values if v is not None]
    if not vals:
        return f'<svg width="{width}" height="{height}"></svg>'
    vmax, vmin = max(vals), min(vals)
    span = (vmax - vmin) or 1
    n = len(vals)
    step = width / (n - 1) if n > 1 else 0
    pts = []
    for i, v in enumerate(vals):
        x = i * step
        y = height - 2 - (v - vmin) / span * (height - 4)
        pts.append(f"{x:.1f},{y:.1f}")
    last_x, last_y = pts[-1].split(",")
    return (
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'preserveAspectRatio="none" role="img">'
        f'<polyline fill="none" stroke="{color}" stroke-width="1.5" '
        f'points="{" ".join(pts)}"/>'
        f'<circle cx="{last_x}" cy="{last_y}" r="2.2" fill="{color}"/></svg>'
    )


def bar_chart(items, *, width=340, bar_color=INK, value_suffix=""):
    items = [(str(label), float(v)) for label, v in items if v is not None]
    if not items:
        return '<p style="color:#6b7c85;font-size:13px;">No data.</p>'
    vmax = max(v for _, v in items) or 1
    row_h, label_w = 22, 120
    bar_w = width - label_w - 52
    rows = []
    for i, (label, value) in enumerate(items):
        w = value / vmax * bar_w
        y = i * row_h
        val_txt = (
            f"{int(value):,}{value_suffix}"
            if value == int(value)
            else f"{value:g}{value_suffix}"
        )
        rows.append(
            f'<text x="0" y="{y + 15}" font-size="12" fill="currentColor">{_esc(label)}</text>'
            f'<rect x="{label_w}" y="{y + 4}" width="{w:.1f}" height="14" rx="2" fill="{bar_color}"/>'
            f'<text x="{label_w + w + 6:.1f}" y="{y + 15}" font-size="11" fill="{MUTED}">{val_txt}</text>'
        )
    height = len(items) * row_h + 4
    return (
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'role="img" style="max-width:100%">' + "".join(rows) + "</svg>"
    )


def line_chart(series, *, width=680, height=240):
    all_pts = [(d, v) for s in series for (d, v) in s["points"] if v is not None]
    if not all_pts:
        return '<p style="color:#6b7c85;font-size:13px;">No data yet.</p>'
    dates = sorted({d for d, _ in all_pts})
    vals = [v for _, v in all_pts]
    pad_l, pad_r, pad_t, pad_b = 40, 14, 14, 26
    pw, ph = width - pad_l - pad_r, height - pad_t - pad_b
    ords = [date.fromisoformat(d).toordinal() for d in dates]
    lo, hi = min(ords), max(ords)
    span = (hi - lo) or 1
    vmax = max(vals)
    vmax = vmax * 1.1 if vmax > 0 else 1

    def sx(d):
        return pad_l + (date.fromisoformat(d).toordinal() - lo) / span * pw

    def sy(v):
        return pad_t + ph - (v / vmax) * ph

    parts = []
    for frac in (0.0, 0.5, 1.0):
        yv = vmax * frac
        y = sy(yv)
        parts.append(
            f'<line x1="{pad_l}" y1="{y:.1f}" x2="{pad_l + pw}" y2="{y:.1f}" '
            f'stroke="{MUTED}" stroke-opacity="0.25" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{pad_l - 6}" y="{y + 3:.1f}" text-anchor="end" '
            f'font-size="10" fill="{MUTED}">{_fmt_num(yv)}</text>'
        )

    for i, s in enumerate(series):
        color = s.get("color", SERIES[i % len(SERIES)])
        segment = []
        for d, v in s["points"]:
            if v is None:
                if segment:
                    parts.append(_polyline(segment, color))
                segment = []
            else:
                segment.append(f"{sx(d):.1f},{sy(v):.1f}")
        if segment:
            parts.append(_polyline(segment, color))
        for d, v in s["points"]:
            if v is not None:
                parts.append(
                    f'<circle cx="{sx(d):.1f}" cy="{sy(v):.1f}" r="2.4" fill="{color}"/>'
                )

    for d in (dates[0], dates[len(dates) // 2], dates[-1]):
        parts.append(
            f'<text x="{sx(d):.1f}" y="{height - 8}" text-anchor="middle" '
            f'font-size="10" fill="{MUTED}">{_esc(d[5:])}</text>'
        )

    return (
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'role="img" style="max-width:100%">' + "".join(parts) + "</svg>"
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_metrics_report.py -v`
Expected: PASS (all tests so far).

- [ ] **Step 5: Format and commit**

```bash
black scripts/metrics_report.py tests/test_metrics_report.py
git add scripts/metrics_report.py tests/test_metrics_report.py
git commit -m "feat: add inline SVG chart primitives for trends report"
```

---

## Task 3: HTML assembly (`build_trends_html`)

**Files:**
- Modify: `scripts/metrics_report.py`
- Test: `tests/test_metrics_report.py`

**Interfaces:**
- Consumes: `sparkline`, `bar_chart`, `line_chart`, `_esc`, `SERIES` from Task 2.
- Produces:
  - `build_trends_html(records, *, generated_at=None) -> str` — full standalone HTML document. Empty records → friendly "no data yet" document. `generated_at` (a `datetime`) is injected so output is deterministic in tests.
  - Helpers: `_get(r, *keys)`, `_points(records, *keys)`, `_pairs(items, label_key, value_key)`, `_kpi_tile`, `_line_card`, `_bar_card`, `_shell`, and the `_CSS` constant.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_metrics_report.py`:

```python
from datetime import datetime, timezone  # add near the top imports if not present

SAMPLE = [
    {
        "week_start": "2026-06-22",
        "week_end": "2026-06-28",
        "github": {"stars": 17, "forks": 2, "views_this_week": 26, "clones_this_week": 29},
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
        "github": {"stars": 17, "forks": 2, "views_this_week": 61, "clones_this_week": 59},
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
    thin = [{"week_start": "2026-01-01"}, {"week_start": "2026-01-08", "github": {"stars": 5}}]
    html = mr.build_trends_html(thin)
    assert "<!DOCTYPE html>" in html


def test_build_trends_html_deterministic_with_fixed_timestamp():
    ts = datetime(2026, 7, 13, 6, 0, 0, tzinfo=timezone.utc)
    a = mr.build_trends_html(SAMPLE, generated_at=ts)
    b = mr.build_trends_html(SAMPLE, generated_at=ts)
    assert a == b
    assert "2026-07-13 06:00 UTC" in a
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_metrics_report.py -k build_trends_html -v`
Expected: FAIL — `AttributeError: module 'metrics_report' has no attribute 'build_trends_html'`.

- [ ] **Step 3: Write minimal implementation**

Add `from datetime import datetime, timezone` to the import block (alongside `date`). Append to `scripts/metrics_report.py`:

```python
_CSS = """
:root{--ink:#1a6b8a;--bg:#f4f6f7;--card:#ffffff;--text:#1f2a30;--muted:#6b7c85;--line:#e6ebee}
@media (prefers-color-scheme:dark){:root{--bg:#0f1518;--card:#182126;--text:#e7edf0;--muted:#93a4ad;--line:#26333a}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--text);font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;line-height:1.5}
header{background:var(--ink);color:#fff;padding:22px 28px}
header h1{margin:0;font-size:22px}
.subtitle{margin:4px 0 0;color:#cfe6ef;font-size:13px}
main{max-width:1040px;margin:0 auto;padding:20px}
h2{font-size:16px;margin:26px 0 10px;color:var(--ink)}
h2 .sub{color:var(--muted);font-weight:400;font-size:13px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px;color:var(--text)}
.kpi-label{font-size:12px;color:var(--muted)}
.kpi-value{font-size:26px;font-weight:700;margin:2px 0}
.delta{font-size:12px;font-weight:600}
.delta.up{color:#2e8b57}.delta.down{color:#c0504d}.delta.flat{color:var(--muted)}
.kpi-spark{margin-top:6px;color:var(--ink)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px;overflow-x:auto;color:var(--text)}
.card h3{margin:0 0 8px;font-size:14px}
.legend{margin-bottom:6px;font-size:12px;color:var(--muted)}
.legend .lg{margin-right:12px;white-space:nowrap}
.legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:middle}
.empty{color:var(--muted);font-size:15px;padding:30px 0}
footer{max-width:1040px;margin:0 auto;padding:16px 20px;color:var(--muted);font-size:12px;border-top:1px solid var(--line)}
svg{max-width:100%;height:auto}
"""


def _get(r, *keys):
    node = r
    for k in keys:
        if not isinstance(node, dict) or k not in node:
            return None
        node = node[k]
    return node


def _points(records, *keys):
    return [(r.get("week_start"), _get(r, *keys)) for r in records]


def _pairs(items, label_key, value_key):
    if not isinstance(items, list):
        return []
    out = []
    for it in items:
        if isinstance(it, dict) and it.get(value_key) is not None:
            out.append((it.get(label_key, "?"), it.get(value_key)))
    return out


def _kpi_tile(records, keys, label):
    vals = [_get(r, *keys) for r in records]
    present = [v for v in vals if v is not None]
    latest = present[-1] if present else None
    prev = present[-2] if len(present) >= 2 else None
    if latest is None:
        value_txt, delta_html = "—", ""
    else:
        value_txt = (
            f"{int(latest):,}" if float(latest) == int(latest) else f"{latest:g}"
        )
        if prev is not None:
            d = latest - prev
            sign = "+" if d >= 0 else "−"
            cls = "up" if d >= 0 else "down"
            delta_html = f'<span class="delta {cls}">{sign}{abs(d):g} wk/wk</span>'
        else:
            delta_html = '<span class="delta flat">new</span>'
    return (
        f'<div class="kpi"><div class="kpi-label">{_esc(label)}</div>'
        f'<div class="kpi-value">{value_txt}</div>{delta_html}'
        f'<div class="kpi-spark">{sparkline(vals)}</div></div>'
    )


def _line_card(title, records, spec):
    series = [
        {"label": label, "color": color, "points": _points(records, *keys)}
        for keys, label, color in spec
    ]
    legend = "".join(
        f'<span class="lg"><i style="background:{s["color"]}"></i>{_esc(s["label"])}</span>'
        for s in series
    )
    return (
        f'<div class="card"><h3>{_esc(title)}</h3>'
        f'<div class="legend">{legend}</div>{line_chart(series)}</div>'
    )


def _bar_card(title, pairs):
    return f'<div class="card"><h3>{_esc(title)}</h3>{bar_chart(pairs)}</div>'


def _shell(body, subtitle, gen_label):
    return (
        "<!DOCTYPE html>\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        "<title>Shedding Hub — Weekly Trends</title>"
        f"<style>{_CSS}</style></head><body>"
        f'<header><h1>Shedding Hub — Weekly Trends</h1>'
        f'<p class="subtitle">{_esc(subtitle)}</p></header>'
        f"<main>{body}</main>"
        f"<footer>Generated {_esc(gen_label)} · source: metrics-data branch</footer>"
        "</body></html>"
    )


def build_trends_html(records, *, generated_at=None):
    generated_at = generated_at or datetime.now(timezone.utc)
    gen_label = generated_at.strftime("%Y-%m-%d %H:%M UTC")

    if not records:
        return _shell(
            '<p class="empty">No weekly metrics recorded yet. This report will '
            "populate once the weekly workflow has run.</p>",
            "No data yet",
            gen_label,
        )

    first = records[0].get("week_start", "?")
    last = records[-1].get("week_end") or records[-1].get("week_start", "?")
    subtitle = f"{first} → {last} · {len(records)} week(s)"

    kpis = "".join(
        [
            _kpi_tile(records, ("github", "stars"), "GitHub stars"),
            _kpi_tile(records, ("pypi", "last_week"), "PyPI downloads / wk"),
            _kpi_tile(records, ("ga4", "active_users"), "Weekly active users"),
            _kpi_tile(records, ("ga4", "page_views"), "Page views / wk"),
        ]
    )

    charts = "".join(
        [
            _line_card(
                "GitHub stars & forks",
                records,
                [
                    (("github", "stars"), "Stars", SERIES[0]),
                    (("github", "forks"), "Forks", SERIES[1]),
                ],
            ),
            _line_card(
                "PyPI downloads (per week)",
                records,
                [(("pypi", "last_week"), "Downloads", SERIES[0])],
            ),
            _line_card(
                "Website users",
                records,
                [
                    (("ga4", "active_users"), "Active", SERIES[0]),
                    (("ga4", "new_users"), "New", SERIES[1]),
                ],
            ),
            _line_card(
                "Page views & engagement",
                records,
                [
                    (("ga4", "page_views"), "Page views", SERIES[0]),
                    (("ga4", "avg_engagement_seconds"), "Avg engagement (s)", SERIES[1]),
                ],
            ),
            _line_card(
                "Repo traffic (this-week counts)",
                records,
                [
                    (("github", "views_this_week"), "Views", SERIES[0]),
                    (("github", "clones_this_week"), "Clones", SERIES[1]),
                ],
            ),
        ]
    )

    latest = records[-1]
    comp = "".join(
        [
            _bar_card(
                "Traffic sources",
                _pairs(_get(latest, "ga4", "traffic_sources"), "source", "sessions"),
            ),
            _bar_card(
                "Top countries",
                _pairs(_get(latest, "ga4", "top_countries"), "country", "active_users"),
            ),
            _bar_card(
                "Devices",
                _pairs(_get(latest, "ga4", "device_types"), "device", "active_users"),
            ),
            _bar_card(
                "Top pages",
                _pairs(_get(latest, "ga4", "page_breakdown"), "page", "views"),
            ),
        ]
    )

    body = (
        f'<div class="kpis">{kpis}</div>'
        f"<h2>Trends over time</h2>"
        f'<div class="grid">{charts}</div>'
        f'<h2>Latest week composition '
        f'<span class="sub">({_esc(latest.get("week_start", ""))} – '
        f'{_esc(latest.get("week_end", ""))})</span></h2>'
        f'<div class="grid">{comp}</div>'
    )
    return _shell(body, subtitle, gen_label)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_metrics_report.py -v`
Expected: PASS (all tests).

- [ ] **Step 5: Format and commit**

```bash
black scripts/metrics_report.py tests/test_metrics_report.py
git add scripts/metrics_report.py tests/test_metrics_report.py
git commit -m "feat: assemble standalone HTML trends dashboard"
```

---

## Task 4: CLI + read-from-branch, and docs

**Files:**
- Modify: `scripts/metrics_report.py`
- Modify: `metrics/README.md`
- Test: `tests/test_metrics_report.py`

**Interfaces:**
- Consumes: `load_records`, `_parse_jsonl`, `build_trends_html`.
- Produces:
  - `read_history_from_ref(ref="origin/metrics-data", relpath="metrics/weekly_metrics.jsonl") -> str | None` — returns the file's text at a git ref via `git show`, or `None` if unavailable.
  - `main(argv=None) -> int` — CLI entry: positional `path` (optional), `--ref` (default `origin/metrics-data`), `-o/--output`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_metrics_report.py`:

```python
import subprocess  # add near the top imports if not already present

SCRIPT = str(Path(__file__).resolve().parent.parent / "scripts" / "metrics_report.py")


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
    src.write_text('{"week_start":"2026-06-22","github":{"stars":17}}\n', encoding="utf-8")
    res = subprocess.run(
        [sys.executable, SCRIPT, str(src)], capture_output=True, text=True, check=True
    )
    assert "<!DOCTYPE html>" in res.stdout


def test_read_history_from_ref_bad_ref_returns_none():
    assert mr.read_history_from_ref("nonexistent/ref-xyz-123") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_metrics_report.py -k "cli or read_history" -v`
Expected: FAIL — CLI has no `__main__` handling / `read_history_from_ref` missing.

- [ ] **Step 3: Write minimal implementation**

Add `import argparse`, `import subprocess`, `import sys` to the top import block of `scripts/metrics_report.py`. Append at the end of the file:

```python
def read_history_from_ref(
    ref="origin/metrics-data", relpath="metrics/weekly_metrics.jsonl"
):
    """Return the text of ``relpath`` at git ``ref``, or ``None`` if unavailable."""
    try:
        out = subprocess.run(
            ["git", "show", f"{ref}:{relpath}"],
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build the Shedding Hub weekly trends report (HTML)."
    )
    parser.add_argument(
        "path",
        nargs="?",
        help="Path to weekly_metrics.jsonl. If omitted, reads from --ref.",
    )
    parser.add_argument(
        "--ref",
        default="origin/metrics-data",
        help="Git ref to read history from when no path is given "
        "(default: origin/metrics-data).",
    )
    parser.add_argument("-o", "--output", help="Output HTML file (default: stdout).")
    args = parser.parse_args(argv)

    if args.path:
        records = load_records(args.path)
    else:
        text = read_history_from_ref(args.ref)
        if text is None:
            print(
                f"Could not read metrics/weekly_metrics.jsonl at ref '{args.ref}'.\n"
                "Fetch it first:  git fetch origin metrics-data",
                file=sys.stderr,
            )
            return 1
        records = _parse_jsonl(text)

    html = build_trends_html(records)
    if args.output:
        Path(args.output).write_text(html, encoding="utf-8")
        print(f"Wrote {args.output} ({len(records)} week(s)).")
    else:
        sys.stdout.write(html)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_metrics_report.py -v`
Expected: PASS (all tests).

- [ ] **Step 5: Document the report in `metrics/README.md`**

Append this section to `metrics/README.md`:

```markdown

## Weekly trends report

Each weekly run also builds a standalone HTML trends dashboard
(`scripts/metrics_report.py`) from the full history and **attaches it to the
weekly email** as `shedding-hub-trends_<week_end>.html`. Open the attachment in a
browser to see stars, downloads, active users, and page views over time, plus the
latest week's traffic-source / country / device / page breakdowns.

To regenerate the report locally from the source-of-truth `metrics-data` branch:

```bash
git fetch origin metrics-data
python scripts/metrics_report.py --ref origin/metrics-data -o trends.html
```

Or from an explicit file:

```bash
git show origin/metrics-data:metrics/weekly_metrics.jsonl > hist.jsonl
python scripts/metrics_report.py hist.jsonl -o trends.html
```

Reading `main`'s stale seed instead of `metrics-data` will produce a report from
incomplete history — always use the `metrics-data` version.
```

- [ ] **Step 6: Format and commit**

```bash
black scripts/metrics_report.py tests/test_metrics_report.py
git add scripts/metrics_report.py tests/test_metrics_report.py metrics/README.md
git commit -m "feat: add trends-report CLI and document local regeneration"
```

---

## Task 5: Attach report to the weekly email

**Files:**
- Modify: `scripts/weekly_report.py`
- Test: `tests/test_metrics_report.py`

**Interfaces:**
- Consumes: `metrics_report.load_records`, `metrics_report.build_trends_html`.
- Produces (in `weekly_report.py`):
  - `build_email_message(html, subject, sender, recipients, attachments=None) -> MIMEMultipart` — builds a `mixed` message with the HTML body nested in an `alternative` part plus zero or more HTML attachments (`attachments` is `list[(filename, content)]`).
  - `send_report(html, attachments=None) -> None` — now forwards attachments.
- Behavior change in `main()`: after `save_metrics`, build the trends report (best-effort) and pass it as an attachment.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_metrics_report.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_metrics_report.py -k build_email_message -v`
Expected: FAIL — importing `weekly_report` fails on `import anthropic` (not installed in the test env), or `build_email_message` is undefined.

- [ ] **Step 3a: Make the `anthropic` import lazy**

In `scripts/weekly_report.py`, remove the top-level `import anthropic` line:

```python
import anthropic
import requests
```

becomes:

```python
import requests
```

Then inside `summarize_with_claude`, add the import as the first line of the function body (just before `if not ANTHROPIC_API_KEY:`):

```python
def summarize_with_claude(ga: dict, gh: dict, pypi: dict) -> str:
    import anthropic

    if not ANTHROPIC_API_KEY:
        return "Claude summarization unavailable (no API key)."
```

- [ ] **Step 3b: Add the attachment-aware email builder**

In `scripts/weekly_report.py`, add to the imports near the other `email` imports:

```python
from email.mime.application import MIMEApplication
```

and add `import metrics_report` after the third-party imports:

```python
import requests

import metrics_report
```

Replace the existing `send_report` function:

```python
def send_report(html: str) -> None:
    recipients = [r.strip() for r in REPORT_TO_EMAIL.split(";") if r.strip()]
    if not recipients:
        raise ValueError("REPORT_TO_EMAIL is not set")

    subject = f"Shedding Hub Weekly Report — {WEEK_LABEL}"
    msg = MIMEMultipart("alternative")
    msg["From"] = SMTP_USER
    msg["To"] = ", ".join(recipients)
    msg["Subject"] = subject
    msg.attach(MIMEText(html, "html", "utf-8"))

    with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
        server.starttls()
        server.login(SMTP_USER, SMTP_PASSWORD)
        server.sendmail(SMTP_USER, recipients, msg.as_string())
```

with:

```python
def build_email_message(html, subject, sender, recipients, attachments=None):
    """Build a multipart email: HTML body plus optional HTML attachments.

    ``attachments`` is a list of ``(filename, content)`` tuples.
    """
    msg = MIMEMultipart("mixed")
    msg["From"] = sender
    msg["To"] = ", ".join(recipients)
    msg["Subject"] = subject

    alternative = MIMEMultipart("alternative")
    alternative.attach(MIMEText(html, "html", "utf-8"))
    msg.attach(alternative)

    for filename, content in attachments or []:
        part = MIMEApplication(content.encode("utf-8"), _subtype="html")
        part.add_header("Content-Disposition", "attachment", filename=filename)
        msg.attach(part)
    return msg


def send_report(html: str, attachments=None) -> None:
    recipients = [r.strip() for r in REPORT_TO_EMAIL.split(";") if r.strip()]
    if not recipients:
        raise ValueError("REPORT_TO_EMAIL is not set")

    subject = f"Shedding Hub Weekly Report — {WEEK_LABEL}"
    msg = build_email_message(html, subject, SMTP_USER, recipients, attachments)

    with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
        server.starttls()
        server.login(SMTP_USER, SMTP_PASSWORD)
        server.sendmail(SMTP_USER, recipients, msg.as_string())
```

- [ ] **Step 3c: Wire the report into `main()`**

In `scripts/weekly_report.py`, replace this block in `main()`:

```python
    print("  → Saving metrics snapshot")
    save_metrics(ga, gh, pypi)

    print("  → Claude weekly narrative")
    summary = summarize_with_claude(ga, gh, pypi)

    print("  → Building and sending email")
    html = build_html(summary, ga, gh, pypi)
    send_report(html)
```

with:

```python
    print("  → Saving metrics snapshot")
    save_metrics(ga, gh, pypi)

    print("  → Building trends report attachment")
    attachments = []
    try:
        records = metrics_report.load_records(METRICS_FILE)
        trends_html = metrics_report.build_trends_html(records)
        attachments.append((f"shedding-hub-trends_{END_LABEL}.html", trends_html))
    except Exception as exc:
        print(f"  → Trends report skipped ({exc}); sending email without it.")

    print("  → Claude weekly narrative")
    summary = summarize_with_claude(ga, gh, pypi)

    print("  → Building and sending email")
    html = build_html(summary, ga, gh, pypi)
    send_report(html, attachments=attachments)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_metrics_report.py -v`
Expected: PASS (all tests, including the two email-builder tests).

- [ ] **Step 5: Run the full suite and formatting check**

Run: `pytest -q && black --check scripts/metrics_report.py scripts/weekly_report.py tests/test_metrics_report.py`
Expected: full test suite passes; black reports no changes needed.

- [ ] **Step 6: Format and commit**

```bash
black scripts/weekly_report.py scripts/metrics_report.py tests/test_metrics_report.py
git add scripts/weekly_report.py tests/test_metrics_report.py
git commit -m "feat: attach weekly trends report to the report email"
```

---

## Task 6: End-to-end verification against real history

**Files:** none (verification only).

- [ ] **Step 1: Generate a report from the real `metrics-data` history**

```bash
git fetch origin metrics-data
python scripts/metrics_report.py --ref origin/metrics-data -o /tmp/trends.html
```

Expected: prints `Wrote /tmp/trends.html (N week(s)).` with N ≥ 4.

- [ ] **Step 2: Eyeball the output**

Open `/tmp/trends.html` in a browser. Confirm: KPI tiles show latest values with week-over-week deltas; the five trend charts render with the late-March→late-June gap visible as an actual gap; the four latest-week composition bar charts render; light and dark themes both look right.

- [ ] **Step 3: Confirm the workflow needs no YAML change**

Re-read `.github/workflows/weekly-report.yaml`. Confirm the "Restore metrics history" step still runs before `python scripts/weekly_report.py`, so the report is built from full history. No workflow edit is expected; note this explicitly in the PR description.

---

## Self-Review

**Spec coverage:**
- Standalone HTML dashboard attached to email → Tasks 3, 5. ✅
- Inline SVG charts (line/bar/sparkline) → Task 2. ✅
- KPI header row with deltas → Task 3 (`_kpi_tile`). ✅
- Trend charts (stars/forks, PyPI, active/new users, page views/engagement, repo views/clones) → Task 3. ✅
- Latest-week composition (sources, countries, devices, pages) → Task 3. ✅
- `metrics-data` as source of truth + local `--ref` regeneration → Task 4. ✅
- Best-effort (email sends without attachment on failure) → Task 5 `main()` try/except. ✅
- Determinism via injected `generated_at` → Task 3 test. ✅
- Error handling (empty/single/missing-keys/gap) → Task 3 tests. ✅
- Testing per spec → `tests/test_metrics_report.py` across tasks. ✅
- Docs → Task 4 `metrics/README.md`. ✅

**Placeholder scan:** No TBD/TODO; all steps contain complete code and exact commands. ✅

**Type consistency:** `build_trends_html(records, *, generated_at=None)`, `line_chart(series, ...)`, `bar_chart(items, ...)`, `sparkline(values, ...)`, `load_records(path)`, `_parse_jsonl(text)`, `read_history_from_ref(ref, relpath)`, `main(argv)`, `build_email_message(html, subject, sender, recipients, attachments=None)`, `send_report(html, attachments=None)` — names and signatures are consistent across all tasks. ✅
