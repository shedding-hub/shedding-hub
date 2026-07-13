# Weekly Trends Report — Design

**Date:** 2026-07-13
**Status:** Approved (design)
**Author:** brainstormed with Claude Code

## Problem

The weekly email (`scripts/weekly_report.py`, run every Monday by
`.github/workflows/weekly-report.yaml`) is a **single-week snapshot**: it shows
the most recent 7 days of GitHub / PyPI / Google Analytics numbers with small
CSS-bar charts, but it shows no history. Meanwhile, `metrics/weekly_metrics.jsonl`
already accumulates one JSON record per week (stored on the `metrics-data`
branch) and grows every Monday.

We want a **trend report** built from that full history — stars, downloads,
active users, and page views over time, with good visualization — delivered
**alongside** the existing weekly email, and reproducible as new weeks arrive.

## Decision

Build a **standalone, self-contained HTML dashboard** and **attach it to the
existing weekly email**. The email body is unchanged; the attachment opens in a
browser, so it can use real inline **SVG** charts (which Gmail strips from email
bodies but renders fine in an opened HTML file).

Chosen approach (of three considered):

- **A — Standalone module, hand-built inline SVG charts (chosen).** Pure,
  dependency-light, deterministic, and unit-testable. The dataviz skill provides
  the chart design system.
- B — matplotlib PNGs base64-embedded. Adds a workflow dependency; image output
  is harder to diff/test and not crisp on retina. Rejected.
- C — Inline JS charting library (Chart.js). Richest interactivity but a large
  non-deterministic inline blob, harder to test. Rejected for v1.

## Data flow

No new plumbing is required. The workflow already restores the **full history**
into the working copy `metrics/weekly_metrics.jsonl` before the script runs, and
`save_metrics()` appends the current week. So at report-build time the local
JSONL holds every week including the current one.

```
weekly-report.yaml
  └─ restore metrics/weekly_metrics.jsonl from metrics-data branch (full history)
  └─ python scripts/weekly_report.py
        ├─ collect_ga4 / collect_github / collect_pypi   (7-day window, unchanged)
        ├─ save_metrics(...)          → appends current week to the JSONL
        ├─ build_trends_html(records) → reads the JSONL, returns standalone HTML   ◀ NEW
        ├─ summarize_with_claude(...) (unchanged)
        └─ send_report(html, attachment=trends_html)  → email body unchanged,      ◀ MODIFIED
                                                         trends HTML attached
  └─ persist JSONL back to metrics-data branch (unchanged)
```

## Components

### New: `scripts/metrics_report.py`

A self-contained report builder with no dependency on the email/collection code.

- `load_records(path) -> list[dict]` — parse JSONL, skip blank/corrupt lines,
  return records sorted ascending by `week_start`.
- SVG chart helpers (pure functions returning SVG strings):
  - `line_chart(series, ...)` — one or more time series over week dates.
  - `bar_chart(items, ...)` — horizontal bars for latest-week composition.
  - `sparkline(values, ...)` — compact inline trend for KPI tiles.
- `build_trends_html(records) -> str` — assembles a full standalone HTML document:
  inline CSS, theme matched to the email's teal (`#1a6b8a`), responsive layout,
  light/dark aware. Returns `""`-safe placeholder content when history is thin.
- `__main__` CLI: `python scripts/metrics_report.py [JSONL_PATH] [-o OUT.html]`
  (defaults: `metrics/weekly_metrics.jsonl`, stdout) for local regeneration.

### Modified: `scripts/weekly_report.py`

- After `save_metrics(...)`, load the JSONL and call `build_trends_html(records)`.
- `send_report(html, ...)` gains an optional HTML attachment
  (`shedding-hub-trends_<week_end>.html`) via `MIMEApplication` /
  `MIMEText(..., "html")` with a `Content-Disposition: attachment` header. The
  existing `MIMEText(html, "html")` body part is unchanged.
- Report generation is best-effort: if `build_trends_html` raises, log and send
  the email **without** the attachment rather than failing the whole run.

## Report content

Built entirely from existing JSONL fields — no new data collection.

- **KPI header row:** latest `github.stars`, `pypi.last_week`, `ga4.active_users`,
  `ga4.page_views`, each with the week-over-week delta and a sparkline.
- **Trend line charts** (x-axis = real `week_start` dates, so gaps in the data
  render honestly rather than being compressed):
  - GitHub stars & forks
  - PyPI weekly downloads (`pypi.last_week`)
  - GA4 active vs. new users
  - Page views & average engagement seconds
  - GitHub repo views & clones (this-week counts)
- **Latest-week composition** (bar charts from the most recent record):
  traffic sources, top countries, device split, top pages (`page_breakdown`).
- **Header/footer:** date range covered, record count, generated-at timestamp.

## Error handling

- Missing/empty JSONL → report renders a friendly "not enough history yet" state;
  never crashes.
- Missing per-record fields (older records lacking a key) → treated as absent
  data points, not errors; charts skip missing points.
- Single record → KPIs render with no delta; trend charts show a single marker.
- Attachment build failure → email still sends (body unchanged), error logged.

## Testing

`tests/test_metrics_report.py`:

- `build_trends_html` returns valid, non-empty HTML for a multi-record fixture
  and contains the expected section headings and `<svg>` marks.
- Edge cases: empty list, single record, records with missing keys, records with
  a date gap — all build without raising.
- `load_records` skips blank and malformed lines and sorts by `week_start`.

Determinism: the same JSONL input yields byte-identical HTML (aside from the
generated-at timestamp, which is injected via a parameter so tests can pin it),
so local CLI regeneration reproduces the emailed report exactly.

## Scope guards (YAGNI)

- No new metrics collection — uses only fields already in the JSONL.
- No hosting / GitHub Pages, no new workflow, no new secrets.
- No interactivity (hover/tooltips) in v1; SVG is static.

## Out of scope / future

- Interactive charts or a hosted public dashboard.
- Anomaly alerts / thresholds in the narrative.
- Backfilling the pre-June history gap.
