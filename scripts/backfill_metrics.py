"""One-off recovery tool: backfill missing weekly-metrics records.

Combines manually-entered GitHub numbers (from a CSV) with GA4 (queried per week
from the Data API using a locally-provided service account) and PyPI (from the
public pypistats time-series), and merges the recovered records into
``metrics/weekly_metrics.jsonl`` in the exact schema ``save_metrics`` writes.

This is a standalone tool: it deliberately does NOT import ``weekly_report.py``,
so the live weekly job is untouched and there is no heavy-import coupling.

Environment (only needed unless ``--no-ga4``):
    GA4_PROPERTY_ID           GA4 numeric property id.
    GA4_SERVICE_ACCOUNT_JSON  Service-account credentials, given either as the
                              raw JSON string or as a path to the JSON file.
                              Keep this file OUTSIDE the repo (or gitignored);
                              it is never printed or committed.

Usage:
    python scripts/backfill_metrics.py            # uses the default paths
    python scripts/backfill_metrics.py --dry-run  # report without writing
    python scripts/backfill_metrics.py --no-ga4   # GitHub + PyPI only
"""

import argparse
import csv
import json
import os
import urllib.request
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

GH_INT_FIELDS = [
    "stars",
    "forks",
    "open_issues",
    "open_prs",
    "views_this_week",
    "unique_visitors_this_week",
    "clones_this_week",
    "unique_cloners_this_week",
    "new_datasets_count",
]

DEFAULT_CSV = "metrics/backfill/github_backfill.csv"
DEFAULT_JSONL = "metrics/weekly_metrics.jsonl"
PACKAGE = "shedding-hub"


# ---------------------------------------------------------------------------
# GitHub CSV (the manually-recovered part)
# ---------------------------------------------------------------------------


def _norm_date(s):
    """Normalize a date to ISO ``YYYY-MM-DD``.

    Accepts ISO already, or the US ``M/D/YYYY`` form Excel tends to write when a
    spreadsheet re-saves the template.
    """
    s = (s or "").strip()
    if not s:
        return s
    if "/" in s:
        m, d, y = s.split("/")
        return f"{int(y):04d}-{int(m):02d}-{int(d):02d}"
    return s


def read_github_csv(path):
    """Read filled rows from the GitHub backfill CSV.

    Rows with a blank ``week_start`` or with every GitHub column blank are
    skipped (unfilled template rows). Blank numeric cells become 0. Dates are
    normalized to ISO ``YYYY-MM-DD`` regardless of how a spreadsheet saved them.
    """
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            week_start = _norm_date(row.get("week_start"))
            week_end = _norm_date(row.get("week_end"))
            if not week_start:
                continue
            if all(not (row.get(k) or "").strip() for k in GH_INT_FIELDS):
                continue
            github = {}
            for k in GH_INT_FIELDS:
                v = (row.get(k) or "").strip()
                github[k] = int(v) if v else 0
            rows.append(
                {"week_start": week_start, "week_end": week_end, "github": github}
            )
    return rows


# ---------------------------------------------------------------------------
# PyPI (public pypistats time-series — no credentials)
# ---------------------------------------------------------------------------


def fetch_pypi_daily(package=PACKAGE, include_mirrors=False, timeout=25):
    """Return {date_str: downloads} of daily PyPI downloads from pypistats.

    ``include_mirrors=False`` matches the live report's ``/recent`` semantics
    (real downloads, excluding mirrors).
    """
    param = "true" if include_mirrors else "false"
    url = f"https://pypistats.org/api/packages/{package}/overall?mirrors={param}"
    req = urllib.request.Request(
        url, headers={"User-Agent": f"shedding-hub-backfill/1.0 ({package})"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = json.load(resp).get("data", [])
    return {row["date"]: row["downloads"] for row in data}


def _range_sum(daily, start, end):
    total = 0
    d, e = date.fromisoformat(start), date.fromisoformat(end)
    while d <= e:
        total += daily.get(d.isoformat(), 0)
        d += timedelta(days=1)
    return total


def pypi_weekly(daily, week_start, week_end):
    """Return (last_week, last_month) sums ending at ``week_end`` (inclusive)."""
    end = date.fromisoformat(week_end)
    month_start = (end - timedelta(days=29)).isoformat()
    return (
        _range_sum(daily, week_start, week_end),
        _range_sum(daily, month_start, week_end),
    )


# ---------------------------------------------------------------------------
# GA4 (queried per week; google libs imported lazily so tests need no creds)
# ---------------------------------------------------------------------------


def _empty_ga4():
    return {
        "active_users": 0,
        "new_users": 0,
        "page_views": 0,
        "avg_engagement_seconds": 0.0,
        "traffic_sources": [],
        "top_countries": [],
        "device_types": [],
        "page_breakdown": [],
        "daily_users": {},
    }


def _load_credentials_info(value):
    """Accept the service-account JSON as a raw string or a path to a file."""
    value = value.strip()
    if os.path.isfile(value):
        with open(value, encoding="utf-8") as f:
            return json.load(f)
    return json.loads(value)


def collect_ga4_for(property_id, credentials_info, start_label, end_label):
    """Query GA4 for a single week (dates inclusive). Returns the stored subset."""
    from google.analytics.data_v1beta import BetaAnalyticsDataClient
    from google.analytics.data_v1beta.types import (
        DateRange,
        Dimension,
        Metric,
        RunReportRequest,
    )
    from google.oauth2 import service_account

    creds = service_account.Credentials.from_service_account_info(
        credentials_info,
        scopes=["https://www.googleapis.com/auth/analytics.readonly"],
    )
    client = BetaAnalyticsDataClient(credentials=creds)
    prop = f"properties/{property_id}"
    date_range = [DateRange(start_date=start_label, end_date=end_label)]
    ga4 = _empty_ga4()

    resp = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            metrics=[
                Metric(name="activeUsers"),
                Metric(name="newUsers"),
                Metric(name="screenPageViews"),
                Metric(name="userEngagementDuration"),
            ],
        )
    )
    if resp.rows:
        v = [mv.value for mv in resp.rows[0].metric_values]
        ga4["active_users"] = int(v[0])
        ga4["new_users"] = int(v[1])
        ga4["page_views"] = int(v[2])
        if ga4["active_users"] > 0:
            ga4["avg_engagement_seconds"] = round(float(v[3]) / ga4["active_users"], 1)

    resp = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            dimensions=[Dimension(name="pagePath")],
            metrics=[
                Metric(name="screenPageViews"),
                Metric(name="activeUsers"),
                Metric(name="userEngagementDuration"),
            ],
            limit=10,
        )
    )
    for row in resp.rows:
        users = int(row.metric_values[1].value)
        eng = float(row.metric_values[2].value)
        ga4["page_breakdown"].append(
            {
                "page": row.dimension_values[0].value,
                "views": int(row.metric_values[0].value),
                "active_users": users,
                "avg_engagement_s": round(eng / users, 1) if users > 0 else 0.0,
            }
        )

    resp = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            dimensions=[Dimension(name="sessionDefaultChannelGroup")],
            metrics=[Metric(name="sessions")],
            limit=10,
        )
    )
    for row in resp.rows:
        ga4["traffic_sources"].append(
            {
                "source": row.dimension_values[0].value,
                "sessions": int(row.metric_values[0].value),
            }
        )

    resp = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            dimensions=[Dimension(name="country")],
            metrics=[Metric(name="activeUsers")],
            limit=10,
        )
    )
    for row in resp.rows:
        ga4["top_countries"].append(
            {
                "country": row.dimension_values[0].value,
                "active_users": int(row.metric_values[0].value),
            }
        )

    resp = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            dimensions=[Dimension(name="deviceCategory")],
            metrics=[Metric(name="activeUsers")],
        )
    )
    for row in resp.rows:
        ga4["device_types"].append(
            {
                "device": row.dimension_values[0].value,
                "active_users": int(row.metric_values[0].value),
            }
        )

    resp = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            dimensions=[Dimension(name="date")],
            metrics=[Metric(name="activeUsers")],
        )
    )
    for row in resp.rows:
        raw = row.dimension_values[0].value  # YYYYMMDD
        ga4["daily_users"][f"{raw[:4]}-{raw[4:6]}-{raw[6:]}"] = int(
            row.metric_values[0].value
        )

    return ga4


# ---------------------------------------------------------------------------
# Record assembly + merge
# ---------------------------------------------------------------------------


def assemble_record(week_start, week_end, github, ga4, pypi, collected_at):
    """Build one record in the exact schema ``save_metrics`` writes."""
    return {
        "week_start": week_start,
        "week_end": week_end,
        "collected_at": collected_at,
        "github": {k: github[k] for k in GH_INT_FIELDS},
        "pypi": {"last_week": pypi["last_week"], "last_month": pypi["last_month"]},
        "ga4": ga4,
    }


def merge_records(existing_path, new_records):
    """Merge new records into the JSONL, keeping existing weeks untouched.

    Returns (merged_sorted, added_week_starts).
    """
    by_week = {}
    p = Path(existing_path)
    if p.exists():
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(rec, dict) and rec.get("week_start"):
                by_week[rec["week_start"]] = rec
    added = []
    for rec in new_records:
        ws = rec["week_start"]
        if ws in by_week:
            continue  # never overwrite an existing (real) record
        by_week[ws] = rec
        added.append(ws)
    merged = sorted(by_week.values(), key=lambda r: r.get("week_start", ""))
    return merged, added


def write_jsonl(path, records):
    with open(path, "w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(description="Backfill missing weekly metrics.")
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument("--jsonl", default=DEFAULT_JSONL)
    parser.add_argument("--package", default=PACKAGE)
    parser.add_argument(
        "--no-ga4", action="store_true", help="Skip GA4 (GitHub + PyPI only)."
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Report without writing the JSONL."
    )
    args = parser.parse_args(argv)

    rows = read_github_csv(args.csv)
    if not rows:
        print(f"No filled rows in {args.csv}; nothing to backfill.")
        return 0
    print(f"{len(rows)} week(s) to backfill from {args.csv}")

    daily = fetch_pypi_daily(args.package)

    cred_info = None
    if not args.no_ga4:
        prop = os.getenv("GA4_PROPERTY_ID", "")
        sa = os.getenv("GA4_SERVICE_ACCOUNT_JSON", "")
        if not prop or not sa:
            print(
                "GA4_PROPERTY_ID / GA4_SERVICE_ACCOUNT_JSON not set. "
                "Set them, or re-run with --no-ga4.",
            )
            return 2
        cred_info = _load_credentials_info(sa)
    else:
        prop = ""

    collected_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    new_records = []
    for row in rows:
        ws, we = row["week_start"], row["week_end"]
        last_week, last_month = pypi_weekly(daily, ws, we)
        if args.no_ga4:
            ga4 = _empty_ga4()
        else:
            try:
                ga4 = collect_ga4_for(prop, cred_info, ws, we)
            except Exception as exc:
                print(f"  {ws}: GA4 query failed ({exc}); recording GitHub+PyPI only.")
                ga4 = _empty_ga4()
        record = assemble_record(
            ws,
            we,
            row["github"],
            ga4,
            {"last_week": last_week, "last_month": last_month},
            collected_at,
        )
        print(
            f"  {ws} -> {we}: pypi wk={last_week} mo={last_month}; "
            f"ga4 active={ga4['active_users']} pageviews={ga4['page_views']}"
        )
        new_records.append(record)

    merged, added = merge_records(args.jsonl, new_records)
    print(f"Added {len(added)} new week(s); {args.jsonl} would total {len(merged)}.")
    if args.dry_run:
        print("Dry run: not writing.")
        return 0
    write_jsonl(args.jsonl, merged)
    print(f"Wrote {args.jsonl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
