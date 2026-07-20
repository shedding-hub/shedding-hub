# Metrics Backfill — Design

**Date:** 2026-07-14
**Status:** Approved (design)

## Problem

`weekly_metrics.jsonl` (source of truth on the `metrics-data` branch) is missing
12 weeks between `2026-03-30` and `2026-06-21` — the weekly job sent the emails
but did not persist those records (this predates the metrics-data persistence
fix). We want to recover those weeks.

## Data sources per field

- **GitHub** (repo traffic `views/unique_visitors/clones/unique_cloners`, plus
  snapshot `stars/forks/open_issues/open_prs/new_datasets_count`): the GitHub
  traffic API only serves the last 14 days and snapshot counts are current-only,
  so historical values exist **only in the sent emails**. Recovered by manual
  entry into `metrics/backfill/github_backfill.csv`.
- **GA4** (active/new users, page views, engagement, traffic sources, countries,
  devices, page breakdown, daily users): re-queryable from the GA4 Data API by
  date range. Run **locally** with a service-account key the user supplies
  (path + `GA4_PROPERTY_ID`), within the property's data-retention window.
- **PyPI** (`last_week`, `last_month`): re-queryable from the public pypistats
  time-series (`/packages/shedding-hub/overall`), which covers the window
  (verified: data from 2026-01-14). No credentials needed.

## Component

`scripts/backfill_metrics.py` — a **self-contained** one-off recovery tool (does
not import or modify the production `weekly_report.py`, so the live weekly job is
untouched and there is no `anthropic` import coupling). It:

1. Reads `metrics/backfill/github_backfill.csv` (rows the user filled; blank rows
   or blank cells are skipped/zeroed).
2. Per week: builds the GitHub block from the CSV row, queries GA4 for that
   week's date range, and computes PyPI `last_week`/`last_month` from the
   pypistats daily series.
3. Assembles a record in the **exact `save_metrics` schema** and merges it into
   `metrics/weekly_metrics.jsonl` (dedup by `week_start`, sorted ascending).

Pure logic (CSV parsing, PyPI aggregation from a daily series, record assembly,
JSONL merge) is unit-tested with fixtures; the live GA4 call is injected so tests
run without credentials.

## Credentials handling

The GA4 service-account JSON stays in a file **outside the repo** (or a
gitignored path); the script reads `GA4_SERVICE_ACCOUNT_JSON` (raw JSON or a file
path) and `GA4_PROPERTY_ID` from the environment at runtime. The key is never
printed, logged, or committed.

## Output & commit

Produces the merged **16-record** `weekly_metrics.jsonl`. The recovered history
is pushed to the source-of-truth **`metrics-data`** branch via an isolated
worktree (the same pattern the weekly workflow uses). The diff is shown and
confirmed **before** any push.

## PyPI semantics

Match the live report's `/recent` values (real downloads, excluding mirrors):
`last_week` = sum of daily downloads over the 7-day `week_start..week_end`
window; `last_month` = sum over the 30 days ending `week_end`. Validated against a
recent week whose value is already stored.

## Scope guards

- GA4 weeks beyond the property's data-retention window may return empty; those
  weeks are recorded with GitHub + PyPI only and flagged in the run output.
- No change to the weekly workflow. This is a one-off tool, kept in-repo to
  document the recovery and allow reuse if gaps recur.
