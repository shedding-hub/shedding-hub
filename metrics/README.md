# Weekly metrics

The Weekly Report workflow (`.github/workflows/weekly-report.yaml`) collects a
snapshot of GitHub, PyPI, and Google Analytics metrics every Monday and appends
one JSON record per week to `weekly_metrics.jsonl`.

Because `main` is a protected branch and rejects direct pushes, the **live,
accumulating history is stored on the dedicated `metrics-data` branch**, not on
`main`. Each weekly run restores the history from `metrics-data`, appends the new
record, and pushes it back to that branch.

To read the full history:

```bash
git fetch origin metrics-data
git show origin/metrics-data:metrics/weekly_metrics.jsonl
```

The copy of `weekly_metrics.jsonl` on `main` is only the initial seed and is not
updated by the workflow.

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
