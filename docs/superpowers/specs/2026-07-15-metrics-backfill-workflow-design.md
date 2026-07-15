# Metrics Backfill Workflow — Design

**Date:** 2026-07-15
**Status:** Approved (design)
**Supersedes:** the "run locally with a user-supplied GA4 key" execution path in
[`2026-07-14-metrics-backfill-design.md`](./2026-07-14-metrics-backfill-design.md).
The recovery tool (`scripts/backfill_metrics.py`) is unchanged; only *how it is
run* changes.

## Problem

The backfill tool needs GA4 credentials to recover the GA4 fields for the 12
missing weeks (`2026-03-30`…`2026-06-21`). The user cannot access the GA4
property locally right now. However, the GA4 service account is already stored as
GitHub Actions secrets (`GA4_PROPERTY_ID`, `GA4_SERVICE_ACCOUNT_JSON`) and used by
the Weekly Report workflow. Running the backfill in CI reuses that secret, so no
local GA4 access — and no manual GA4 data entry — is needed.

Key enabling fact: the secret names are the **exact** environment variables
`backfill_metrics.py` already reads, so the script needs **no changes**.

## Component

New workflow `.github/workflows/metrics-backfill.yaml`, triggered by
`workflow_dispatch` **only** (never scheduled, never on push). One input:

- `dry_run` (boolean, **default `true`**): when true, run the tool with
  `--dry-run` (preview in the log, write nothing) and skip the persist step.

### Job steps (reuse the Weekly Report scaffolding)

1. `actions/checkout@v4` — checks out the dispatched ref, so the run reads that
   branch's filled `metrics/backfill/github_backfill.csv`.
2. `actions/setup-python@v5` (3.11) + `pip install google-analytics-data
   google-auth`. PyPI (urllib) and CSV parsing use only the stdlib.
3. **Restore history** — identical to the weekly job: fetch `metrics-data` and
   copy its `weekly_metrics.jsonl` into the working tree. This is the true merge
   base (source of truth), not the dispatched branch's stale copy.
4. **Run backfill** — `env:` supplies `GA4_PROPERTY_ID` /
   `GA4_SERVICE_ACCOUNT_JSON` from secrets. Runs
   `python scripts/backfill_metrics.py --dry-run` when `dry_run` is true, else
   without the flag.
5. **Persist** — gated `if: ${{ !inputs.dry_run }}`. Push the merged JSONL to
   `metrics-data` via the same isolated-worktree pattern as the weekly job, with
   commit message `chore: backfill missing weekly metrics (2026-03-30..2026-06-21)`.

`permissions: contents: write` (needed to push `metrics-data`), mirroring the
weekly job.

## What lands where

`workflow_dispatch` requires the workflow to exist on the default branch to be
runnable, and `main` does not yet have the backfill tool at all (it is stacked on
the in-flight trends-report PR #158, `feature/weekly-trends-report`). The 3
backfill-unique commits touch **only new files** (`scripts/backfill_metrics.py`,
`tests/test_backfill_metrics.py`, `metrics/backfill/*`, the two specs), all
stdlib-only and independent of the trends work, so they land on `main` as a clean
**adds-only PR** with no conflict with #158.

- `main` receives the whole self-contained backfill unit via an **independent PR
  off `main`**: the tool, tests, `metrics/backfill/` (README, `.gitignore`, and
  the **filled** `github_backfill.csv` as the recovery record), the specs, and the
  new workflow file. This is a reusable manual-backfill tool for any future gap.
- Because the filled CSV is committed to `main`, the workflow is dispatched
  **against `main`** — no throwaway branch. `main`'s CSV documents exactly what
  was recovered.

## Execution sequence

1. Author and validate `.github/workflows/metrics-backfill.yaml`; commit it and
   the filled CSV on the working branch.
2. Assemble an adds-only branch off `origin/main` containing only the 7 backfill
   files + the workflow (final content), verify it is additions-only vs `main`
   and that `pytest`/`black` pass, and open a PR to `main`.
3. After the PR merges, dispatch against `main` with `dry_run=true`; confirm the
   log shows the 12 weeks with non-zero GA4/PyPI.
4. Re-dispatch with `dry_run=false`.
5. Verify `metrics-data` lists 16 records.

## Safety

- `merge_records` never overwrites an existing week, so a repeat dispatch is a
  harmless no-op push — real data cannot be corrupted.
- `dry_run` defaults to `true`, so the writing path is never the first action.
- The service-account JSON is passed only through the step `env` from secrets;
  the tool never prints or logs it.

## Verification

`git show origin/metrics-data:metrics/weekly_metrics.jsonl` lists the seed week,
the 12 backfilled weeks (`2026-03-30`…`2026-06-15`), and the 3 recent weeks — 16
total, each backfilled week with non-zero GA4 `active_users`/`page_views`.

## Out of scope

- No change to `scripts/backfill_metrics.py` or to the Weekly Report workflow.
- The abandoned alternative (hand-entering GA4 numbers into the CSV) is dropped:
  CI access to GA4 makes it unnecessary.
