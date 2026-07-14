# Metrics backfill (one-off recovery)

Recovers the weekly-metrics records that were emailed but never persisted to the
`metrics-data` branch (weeks `2026-03-30` … `2026-06-21`).

## 1. Fill `github_backfill.csv`

Only the GitHub numbers must be entered by hand — they can't be re-fetched
(GitHub's traffic API only serves the last 14 days, and star/fork/issue/PR counts
are current-only). Open each week's report email and copy the numbers:

| CSV column | Email field |
|---|---|
| `stars`, `forks`, `open_issues`, `open_prs` | "GitHub Repository" section (Stars / Forks / Open issues / Open PRs) |
| `views_this_week` | "Repo views this week" (the number before "(+N vs last week)") |
| `unique_visitors_this_week` | "Unique visitors this week" |
| `clones_this_week` | "Repo clones this week" |
| `unique_cloners_this_week` | "Unique cloners this week" |
| `new_datasets_count` | "New datasets added" (count; `0` if none) |

Notes:
- `week_start` / `week_end` are pre-filled — don't change them.
- Blank cell → treated as `0`. Delete any week's row if you don't have its email.

GA4 (users, page views, sources, countries, devices, pages) and PyPI downloads
are fetched automatically — you do **not** enter those.

## 2. Provide GA4 credentials (kept out of git)

Save the GA4 service-account key **outside the repo** (or anywhere in this folder
— `*.json` here is gitignored), then set:

```bash
export GA4_PROPERTY_ID="<numeric property id>"
export GA4_SERVICE_ACCOUNT_JSON="/absolute/path/to/service-account.json"   # path or raw JSON
```

## 3. Run

```bash
python scripts/backfill_metrics.py --dry-run   # preview, writes nothing
python scripts/backfill_metrics.py             # merge into metrics/weekly_metrics.jsonl
```

The recovered records are merged in (existing weeks are never overwritten). The
updated `weekly_metrics.jsonl` is then pushed to the `metrics-data` branch after
you review the diff.

Run `--no-ga4` to backfill GitHub + PyPI only (leaves GA4 fields empty).
