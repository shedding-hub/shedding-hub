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
from datetime import date
from pathlib import Path

# --- palette (validated against the dataviz skill's contrast rules) ----------
INK = "#1a6b8a"  # brand teal (matches the email header)
ACCENT = "#e07b39"  # warm secondary
MUTED = "#6b7c85"  # gridlines / axis text
SERIES = ["#1a6b8a", "#e07b39", "#4c9a6a", "#9a4c8a"]


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


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


def _parse_jsonl(text):
    """Parse JSONL text into records sorted by ``week_start``.

    Blank lines and lines that fail to parse as JSON are skipped.
    Non-object JSON values are also skipped.
    """
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(record, dict):
            continue
        records.append(record)
    records.sort(key=lambda r: r.get("week_start", ""))
    return records


def load_records(path):
    """Load records from a JSONL file. Missing file returns ``[]``."""
    p = Path(path)
    if not p.exists():
        return []
    return _parse_jsonl(p.read_text(encoding="utf-8"))


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
