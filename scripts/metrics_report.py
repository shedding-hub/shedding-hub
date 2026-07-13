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
from datetime import date, datetime, timezone
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
        f"<header><h1>Shedding Hub — Weekly Trends</h1>"
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
                    (
                        ("ga4", "avg_engagement_seconds"),
                        "Avg engagement (s)",
                        SERIES[1],
                    ),
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
        f"<h2>Latest week composition "
        f'<span class="sub">({_esc(latest.get("week_start", ""))} – '
        f'{_esc(latest.get("week_end", ""))})</span></h2>'
        f'<div class="grid">{comp}</div>'
    )
    return _shell(body, subtitle, gen_label)
