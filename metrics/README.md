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
