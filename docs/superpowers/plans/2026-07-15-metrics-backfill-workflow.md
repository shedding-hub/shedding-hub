# Metrics Backfill Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover the 12 missing weekly-metrics records (`2026-03-30`…`2026-06-21`) by running the existing `backfill_metrics.py` in CI, using the GA4 service-account already stored as GitHub Actions secrets, and push the result to the `metrics-data` branch.

**Architecture:** Add a `workflow_dispatch`-only workflow `.github/workflows/metrics-backfill.yaml` that reuses the Weekly Report job's proven "restore history from `metrics-data` → run script → push back via isolated worktree" scaffolding, but runs `scripts/backfill_metrics.py` (unchanged) instead of `weekly_report.py`. The backfill tool + workflow + filled CSV land on `main` via one independent, adds-only PR; the workflow is then dispatched against `main` (dry-run first, then real).

**Tech Stack:** GitHub Actions (YAML), Python 3.11, `google-analytics-data` + `google-auth` (GA4), `gh` CLI (PR + dispatch), Git worktrees.

## Global Constraints

- **The backfill script is NOT modified.** It already reads `GA4_PROPERTY_ID` and `GA4_SERVICE_ACCOUNT_JSON` from the environment and merges without overwriting existing weeks. Any temptation to edit it means the plan is wrong — stop and reconsider.
- Secret names are exact: `GA4_PROPERTY_ID`, `GA4_SERVICE_ACCOUNT_JSON` (as used by `.github/workflows/weekly-report.yaml`).
- The workflow trigger is `workflow_dispatch` **only** — never `schedule`, never `push`.
- `dry_run` input defaults to `true`; the persist (write) step is gated on `if: ${{ !inputs.dry_run }}`.
- `main` is a protected branch: reaching it requires a PR, not a direct push.
- The PR to `main` must be **adds-only** (all new files) so it cannot conflict with the in-flight trends-report PR #158.
- Source of truth for metrics history is the `metrics-data` branch, not `main`.
- Target JSONL after the real run: **16 records** (1 seed + 12 backfilled + 3 recent).

---

### Task 1: Author and validate the backfill workflow file

**Files:**
- Create: `.github/workflows/metrics-backfill.yaml`
- Reference (do not modify): `.github/workflows/weekly-report.yaml`, `scripts/backfill_metrics.py`

**Interfaces:**
- Consumes: repo secrets `GA4_PROPERTY_ID`, `GA4_SERVICE_ACCOUNT_JSON`; the `metrics-data` branch; `scripts/backfill_metrics.py` (default `--csv metrics/backfill/github_backfill.csv`, default `--jsonl metrics/weekly_metrics.jsonl`).
- Produces: a dispatchable workflow named `Metrics Backfill (one-off)` with a boolean `dry_run` input (default `true`).

- [ ] **Step 1: Write the workflow file**

Create `.github/workflows/metrics-backfill.yaml` with exactly this content:

```yaml
name: Metrics Backfill (one-off)

# One-off / on-demand recovery of missing weekly-metrics records. Fill
# metrics/backfill/github_backfill.csv with the GitHub numbers from the report
# emails, then run this workflow. GA4 and PyPI are fetched automatically; GA4
# uses the same service-account secret as the Weekly Report workflow.
on:
  workflow_dispatch:
    inputs:
      dry_run:
        description: "Preview only — report the weeks in the log, do NOT write to metrics-data."
        type: boolean
        default: true

jobs:
  backfill:
    name: Backfill weekly metrics
    runs-on: ubuntu-latest
    permissions:
      contents: write

    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"

      - name: Install dependencies
        run: pip install google-analytics-data google-auth

      - name: Restore metrics history from metrics-data branch
        run: |
          # The accumulating snapshot history lives on the unprotected metrics-data
          # branch (main is protected). Seed the working copy with the full history
          # so the backfill merges into it rather than main's stale one-record seed.
          if git ls-remote --exit-code --heads origin metrics-data >/dev/null 2>&1; then
            git fetch --depth=1 origin metrics-data
            if git show FETCH_HEAD:metrics/weekly_metrics.jsonl > "$RUNNER_TEMP/history.jsonl" 2>/dev/null; then
              mkdir -p metrics
              cp "$RUNNER_TEMP/history.jsonl" metrics/weekly_metrics.jsonl
              echo "Restored $(grep -c . metrics/weekly_metrics.jsonl) record(s) from metrics-data."
            else
              echo "metrics-data exists but has no metrics file yet; keeping main's seed."
            fi
          else
            echo "metrics-data branch does not exist yet; using main's file as the seed."
          fi

      - name: Run backfill
        env:
          GA4_PROPERTY_ID: ${{ secrets.GA4_PROPERTY_ID }}
          GA4_SERVICE_ACCOUNT_JSON: ${{ secrets.GA4_SERVICE_ACCOUNT_JSON }}
        run: |
          if [ "${{ inputs.dry_run }}" = "true" ]; then
            echo "DRY RUN — reporting only, no write."
            python scripts/backfill_metrics.py --dry-run
          else
            python scripts/backfill_metrics.py
          fi

      - name: Persist metrics snapshot to metrics-data branch
        if: ${{ !inputs.dry_run }}
        run: |
          git config user.name "github-actions[bot]"
          git config user.email "github-actions[bot]@users.noreply.github.com"

          # The backfill merged the recovered weeks into the history we restored
          # above. Commit it onto the dedicated, unprotected metrics-data branch via
          # an isolated worktree so the protected main branch is never touched.
          cp metrics/weekly_metrics.jsonl "$RUNNER_TEMP/weekly_metrics.jsonl"

          if git ls-remote --exit-code --heads origin metrics-data >/dev/null 2>&1; then
            git fetch --depth=1 origin metrics-data
            git worktree add -B metrics-data "$RUNNER_TEMP/md" FETCH_HEAD
          else
            git worktree add --detach "$RUNNER_TEMP/md"
            git -C "$RUNNER_TEMP/md" checkout --orphan metrics-data
            git -C "$RUNNER_TEMP/md" rm -rf . >/dev/null 2>&1 || true
          fi

          mkdir -p "$RUNNER_TEMP/md/metrics"
          cp "$RUNNER_TEMP/weekly_metrics.jsonl" "$RUNNER_TEMP/md/metrics/weekly_metrics.jsonl"
          git -C "$RUNNER_TEMP/md" add metrics/weekly_metrics.jsonl
          if git -C "$RUNNER_TEMP/md" diff --cached --quiet; then
            echo "No new metrics to commit."
          else
            git -C "$RUNNER_TEMP/md" commit -m "chore: backfill missing weekly metrics (2026-03-30..2026-06-21)"
            git -C "$RUNNER_TEMP/md" push origin metrics-data
            echo "Pushed backfilled metrics to metrics-data."
          fi
```

- [ ] **Step 2: Validate the YAML parses and has the required structure**

Run this from the repo root:

```bash
python - <<'PY'
import yaml
wf = yaml.safe_load(open(".github/workflows/metrics-backfill.yaml", encoding="utf-8"))
# PyYAML parses the bare `on:` key as boolean True — accept either form.
trigger = wf.get("on", wf.get(True))
assert "workflow_dispatch" in trigger, "must be workflow_dispatch"
assert set(trigger) == {"workflow_dispatch"}, "workflow_dispatch ONLY (no schedule/push)"
dr = trigger["workflow_dispatch"]["inputs"]["dry_run"]
assert dr["type"] == "boolean" and dr["default"] is True, "dry_run must default true"
job = wf["jobs"]["backfill"]
assert job["permissions"]["contents"] == "write"
names = [s.get("name", s.get("uses", "")) for s in job["steps"]]
assert "Run backfill" in names and "Persist metrics snapshot to metrics-data branch" in names
persist = next(s for s in job["steps"] if s.get("name", "").startswith("Persist"))
assert persist["if"].strip() == "${{ !inputs.dry_run }}", "persist must be gated on !dry_run"
run = next(s for s in job["steps"] if s.get("name") == "Run backfill")
assert run["env"] == {
    "GA4_PROPERTY_ID": "${{ secrets.GA4_PROPERTY_ID }}",
    "GA4_SERVICE_ACCOUNT_JSON": "${{ secrets.GA4_SERVICE_ACCOUNT_JSON }}",
}, "GA4 secrets must be wired into the run step env"
print("workflow structure OK")
PY
```

Expected output: `workflow structure OK` (no `AssertionError`).

- [ ] **Step 3: Confirm the backfill script still runs unchanged (local sanity)**

Run: `python scripts/backfill_metrics.py --dry-run --no-ga4`
Expected: prints `12 week(s) to backfill …`, one line per week `2026-03-30`…`2026-06-15` with `pypi wk=… mo=…`, and ends `Added 12 new week(s); metrics/weekly_metrics.jsonl would total 16.` / `Dry run: not writing.` (This proves the tool is untouched and the filled CSV parses. `--no-ga4` is used only because there are no local GA4 creds; CI will run the full path.)

- [ ] **Step 4: Commit the workflow and the filled CSV on the working branch**

The filled CSV is currently an uncommitted working-tree change; commit it so it is retrievable by exact path in Task 2.

```bash
git add .github/workflows/metrics-backfill.yaml metrics/backfill/github_backfill.csv
git commit -m "feat: add one-off metrics-backfill CI workflow; fill github_backfill.csv"
```

---

### Task 2: Assemble the adds-only PR branch off `main` and open the PR

**Files:**
- Create branch `feat/metrics-backfill-ci` off `origin/main`, containing only:
  - `scripts/backfill_metrics.py`
  - `tests/test_backfill_metrics.py`
  - `metrics/backfill/.gitignore`
  - `metrics/backfill/README.md`
  - `metrics/backfill/github_backfill.csv` (filled)
  - `docs/superpowers/specs/2026-07-14-metrics-backfill-design.md`
  - `docs/superpowers/specs/2026-07-15-metrics-backfill-workflow-design.md`
  - `.github/workflows/metrics-backfill.yaml`

**Interfaces:**
- Consumes: the committed final content of the above files on `chore/metrics-backfill` (from Task 1 + earlier commits).
- Produces: an open GitHub PR into `main` that is additions-only.

- [ ] **Step 1: Refresh `main` and branch from it**

```bash
git fetch origin main
git switch -c feat/metrics-backfill-ci origin/main
```

Expected: a new branch whose working tree matches `origin/main` (none of the backfill files present yet).

- [ ] **Step 2: Bring in exactly the backfill files from the working branch**

```bash
git checkout chore/metrics-backfill -- \
  scripts/backfill_metrics.py \
  tests/test_backfill_metrics.py \
  metrics/backfill/.gitignore \
  metrics/backfill/README.md \
  metrics/backfill/github_backfill.csv \
  docs/superpowers/specs/2026-07-14-metrics-backfill-design.md \
  docs/superpowers/specs/2026-07-15-metrics-backfill-workflow-design.md \
  .github/workflows/metrics-backfill.yaml
```

- [ ] **Step 3: Verify the change set is additions-only vs `main`**

```bash
git status --short
git diff --cached --name-status origin/main
```

Expected: every line is `A` (added). If any line is `M` or `D`, STOP — a file already exists on `main` and the PR is no longer adds-only; investigate before continuing.

- [ ] **Step 4: Verify tests and formatting pass (build.yaml gates the PR)**

```bash
pytest tests/test_backfill_metrics.py -v
black --check scripts/backfill_metrics.py tests/test_backfill_metrics.py
```

Expected: pytest all-pass; `black --check` reports the two files already formatted (exit 0). If black reports would-reformat, run `black scripts/backfill_metrics.py tests/test_backfill_metrics.py`, re-run `--check`, and re-stage.

- [ ] **Step 5: Commit and push**

```bash
git commit -m "feat: one-off metrics backfill tool + CI workflow

Recover the 12 weekly-metrics records (2026-03-30..2026-06-21) that were emailed
but never persisted. Self-contained tool (GitHub numbers from CSV, GA4 + PyPI
fetched automatically) plus a workflow_dispatch workflow that runs it in CI using
the existing GA4 secret and pushes to the metrics-data branch."
git push -u origin feat/metrics-backfill-ci
```

- [ ] **Step 6: Open the PR into `main`**

```bash
gh pr create --base main --head feat/metrics-backfill-ci \
  --title "feat: one-off metrics backfill tool + CI workflow" \
  --body "Adds \`scripts/backfill_metrics.py\` (recovers the 12 missing weekly-metrics records for 2026-03-30..2026-06-21) and a \`workflow_dispatch\` workflow that runs it in CI with the stored GA4 secret and pushes to \`metrics-data\`. Adds-only; independent of trends-report PR #158. After merge, dispatch the workflow (dry-run first) per docs/superpowers/plans/2026-07-15-metrics-backfill-workflow.md."
```

Expected: prints the new PR URL. Report the URL. **Merging is gated by branch protection / review — this is the checkpoint where a human merges.**

---

### Task 3: Dispatch the backfill and verify `metrics-data`

**Precondition:** the Task 2 PR is merged into `main` (so the workflow exists on the default branch and `main` has the tool + filled CSV). Confirm with `gh workflow list | grep -i "Metrics Backfill"` before dispatching.

**Interfaces:**
- Consumes: the merged workflow on `main`; repo secrets; the `metrics-data` branch.
- Produces: 12 new records on `metrics-data` (16 total).

- [ ] **Step 1: Dry-run dispatch**

```bash
gh workflow run "metrics-backfill.yaml" --ref main -f dry_run=true
sleep 5
gh run list --workflow "metrics-backfill.yaml" --limit 1
```

Then watch the latest run to completion:

```bash
gh run watch "$(gh run list --workflow metrics-backfill.yaml --limit 1 --json databaseId -q '.[0].databaseId')"
```

Expected: run succeeds; the "Run backfill" step log shows `12 week(s) to backfill`, each week `2026-03-30`…`2026-06-15` with **non-zero** `ga4 active=… pageviews=…` and `pypi wk=…`, and `Added 12 new week(s); … would total 16.` The "Persist…" step is **skipped** (dry_run). If any backfilled week shows `ga4 active=0`, note it (may be outside GA4 retention) but this is not a failure.

- [ ] **Step 2: Real dispatch**

```bash
gh workflow run "metrics-backfill.yaml" --ref main -f dry_run=false
sleep 5
gh run watch "$(gh run list --workflow metrics-backfill.yaml --limit 1 --json databaseId -q '.[0].databaseId')"
```

Expected: run succeeds; "Persist…" step runs and prints `Pushed backfilled metrics to metrics-data.`

- [ ] **Step 3: Verify the recovered history**

```bash
git fetch origin metrics-data
git show origin/metrics-data:metrics/weekly_metrics.jsonl | python - <<'PY'
import sys, json
recs = [json.loads(l) for l in sys.stdin if l.strip()]
starts = [r["week_start"] for r in recs]
print("total records:", len(recs))
backfilled = ["2026-03-30","2026-04-06","2026-04-13","2026-04-20","2026-04-27",
              "2026-05-04","2026-05-11","2026-05-18","2026-05-25","2026-06-01",
              "2026-06-08","2026-06-15"]
assert starts == sorted(starts), "records must be sorted by week_start"
assert len(starts) == len(set(starts)), "duplicate week_start(s) present"
assert "2026-03-26" in starts, "seed week missing"
missing = [w for w in backfilled if w not in starts]
assert not missing, f"gap NOT fully filled; missing: {missing}"
# Expect >= 16 (1 seed + 12 backfilled + 3 recent; a later scheduled run may add more).
assert len(recs) >= 16, f"expected at least 16 records, got {len(recs)}"
zero_ga4 = [r["week_start"] for r in recs
            if r["week_start"] in backfilled
            and r.get("ga4", {}).get("active_users", 0) == 0]
print("backfilled weeks with zero GA4 active_users:", zero_ga4 or "none")
print(f"OK: {len(recs)} records, all 12 gap weeks present")
PY
```

Expected: `total records: 16` (or more if a scheduled run has since added a newer week), and `OK: … records, all 12 gap weeks present`. Report any weeks listed under "zero GA4 active_users" (retention-window gaps, acceptable).

- [ ] **Step 4: Idempotency spot-check (optional)**

Re-dispatch with `dry_run=false` once more. Expected: the "Persist…" step prints `No new metrics to commit.` (the 12 weeks already exist; `merge_records` never overwrites). This confirms a repeat run cannot corrupt data.
