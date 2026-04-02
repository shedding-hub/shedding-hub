"""
Shedding Hub Weekly Operations Report

Collects 7-day data from GitHub, PyPI, and Google Analytics,
then uses Claude to write a weekly briefing and sends an HTML email.

Usage:
    python scripts/weekly_report.py

Required environment variables:
    ANTHROPIC_API_KEY        - Claude API key for narrative generation
    SMTP_USER                - Gmail address (sender)
    SMTP_PASSWORD            - Gmail App Password
    REPORT_TO_EMAIL          - Recipient(s), semicolon-separated
    GA4_PROPERTY_ID          - GA4 numeric property ID
    GA4_SERVICE_ACCOUNT_JSON - GCP service account JSON (raw string)
    GH_TOKEN                 - GitHub PAT with repo scope (for traffic API)

Optional environment variables (have defaults):
    SMTP_HOST  - default: smtp.gmail.com
    SMTP_PORT  - default: 587
"""

import json
import os
import re
import smtplib
from datetime import datetime, timedelta, timezone
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path

import anthropic
import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ANTHROPIC_API_KEY = os.getenv("ANTHROPIC_API_KEY", "")
SMTP_HOST = os.getenv("SMTP_HOST", "smtp.gmail.com")
SMTP_PORT = int(os.getenv("SMTP_PORT", "587"))
SMTP_USER = os.getenv("SMTP_USER", "")
SMTP_PASSWORD = os.getenv("SMTP_PASSWORD", "")
GA4_PROPERTY_ID = os.getenv("GA4_PROPERTY_ID", "")
GA4_SERVICE_ACCOUNT_JSON = os.getenv("GA4_SERVICE_ACCOUNT_JSON", "")
GH_TOKEN = os.getenv("GH_TOKEN", "")
REPORT_TO_EMAIL = os.getenv("REPORT_TO_EMAIL", "")

GITHUB_REPO = "shedding-hub/shedding-hub"
PYPI_PACKAGE = "shedding-hub"

_now = datetime.now(timezone.utc)
WEEK_END = _now.replace(hour=0, minute=0, second=0, microsecond=0)
WEEK_START = WEEK_END - timedelta(days=7)
WEEK_END_INCL = WEEK_END - timedelta(seconds=1)

START_LABEL = WEEK_START.strftime("%Y-%m-%d")
END_LABEL = WEEK_END_INCL.strftime("%Y-%m-%d")
WEEK_LABEL = f"{START_LABEL} to {END_LABEL}"

GH_HEADERS = {
    "Authorization": f"token {GH_TOKEN}",
    "Accept": "application/vnd.github+json",
    "X-GitHub-Api-Version": "2022-11-28",
}


# ---------------------------------------------------------------------------
# 1. Google Analytics 4
# ---------------------------------------------------------------------------


def collect_ga4() -> dict:
    result = {
        "active_users": 0,
        "new_users": 0,
        "page_views": 0,
        "avg_engagement_seconds": 0.0,
        "page_breakdown": [],  # [{page, views, active_users, avg_engagement_s}]
        "outbound_clicks": [],  # [{url, clicks}]
        "traffic_sources": [],  # [{source, sessions}]
        "top_countries": [],  # [{country, active_users}]
        "device_types": [],  # [{device, active_users}]
        "daily_users": {},  # date → active_users
        "error": None,
    }

    if not GA4_PROPERTY_ID or not GA4_SERVICE_ACCOUNT_JSON:
        result["error"] = "GA4 credentials not set"
        return result

    try:
        from google.analytics.data_v1beta import BetaAnalyticsDataClient
        from google.analytics.data_v1beta.types import (
            DateRange,
            Dimension,
            FilterExpression,
            Filter,
            Metric,
            RunReportRequest,
        )
        from google.oauth2 import service_account

        credentials = service_account.Credentials.from_service_account_info(
            json.loads(GA4_SERVICE_ACCOUNT_JSON),
            scopes=["https://www.googleapis.com/auth/analytics.readonly"],
        )
        client = BetaAnalyticsDataClient(credentials=credentials)
        prop = f"properties/{GA4_PROPERTY_ID}"
        date_range = [DateRange(start_date=START_LABEL, end_date=END_LABEL)]

        # Overall totals
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
            result["active_users"] = int(v[0])
            result["new_users"] = int(v[1])
            result["page_views"] = int(v[2])
            engagement_total = float(v[3])
            if result["active_users"] > 0:
                result["avg_engagement_seconds"] = round(
                    engagement_total / result["active_users"], 1
                )

        # Per-page breakdown (top 10 pages)
        resp2 = client.run_report(
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
        for row in resp2.rows:
            page = row.dimension_values[0].value
            views = int(row.metric_values[0].value)
            users = int(row.metric_values[1].value)
            eng = float(row.metric_values[2].value)
            avg_eng = round(eng / users, 1) if users > 0 else 0.0
            result["page_breakdown"].append(
                {
                    "page": page,
                    "views": views,
                    "active_users": users,
                    "avg_engagement_s": avg_eng,
                }
            )

        # Outbound link clicks (requires Enhanced Measurement)
        try:
            resp3 = client.run_report(
                RunReportRequest(
                    property=prop,
                    date_ranges=date_range,
                    dimensions=[Dimension(name="linkUrl")],
                    metrics=[Metric(name="eventCount")],
                    dimension_filter=FilterExpression(
                        filter=Filter(
                            field_name="eventName",
                            string_filter=Filter.StringFilter(value="click"),
                        )
                    ),
                    limit=15,
                )
            )
            for row in resp3.rows:
                url = row.dimension_values[0].value
                clicks = int(row.metric_values[0].value)
                result["outbound_clicks"].append({"url": url, "clicks": clicks})
        except Exception:
            pass  # Outbound click tracking may not be configured

        # Traffic sources
        resp4 = client.run_report(
            RunReportRequest(
                property=prop,
                date_ranges=date_range,
                dimensions=[Dimension(name="sessionDefaultChannelGroup")],
                metrics=[Metric(name="sessions")],
                limit=10,
            )
        )
        for row in resp4.rows:
            result["traffic_sources"].append(
                {
                    "source": row.dimension_values[0].value,
                    "sessions": int(row.metric_values[0].value),
                }
            )

        # Geographic reach (top 10 countries)
        resp5 = client.run_report(
            RunReportRequest(
                property=prop,
                date_ranges=date_range,
                dimensions=[Dimension(name="country")],
                metrics=[Metric(name="activeUsers")],
                limit=10,
            )
        )
        for row in resp5.rows:
            result["top_countries"].append(
                {
                    "country": row.dimension_values[0].value,
                    "active_users": int(row.metric_values[0].value),
                }
            )

        # Device type
        resp6 = client.run_report(
            RunReportRequest(
                property=prop,
                date_ranges=date_range,
                dimensions=[Dimension(name="deviceCategory")],
                metrics=[Metric(name="activeUsers")],
            )
        )
        for row in resp6.rows:
            result["device_types"].append(
                {
                    "device": row.dimension_values[0].value,
                    "active_users": int(row.metric_values[0].value),
                }
            )

        # Daily active users sparkline
        resp7 = client.run_report(
            RunReportRequest(
                property=prop,
                date_ranges=date_range,
                dimensions=[Dimension(name="date")],
                metrics=[Metric(name="activeUsers")],
            )
        )
        for row in resp7.rows:
            raw = row.dimension_values[0].value  # YYYYMMDD
            date = f"{raw[:4]}-{raw[4:6]}-{raw[6:]}"
            result["daily_users"][date] = int(row.metric_values[0].value)

    except Exception as exc:
        result["error"] = str(exc)

    return result


# ---------------------------------------------------------------------------
# 2. GitHub repository metrics
# ---------------------------------------------------------------------------


def _gh_get(path: str, params: dict = None) -> requests.Response:
    url = f"https://api.github.com/repos/{GITHUB_REPO}/{path}"
    return requests.get(url, headers=GH_HEADERS, params=params, timeout=20)


def _split_traffic_by_week(items: list) -> tuple[int, int, int, int]:
    """Split 14-day GitHub traffic items into this week vs last week totals.
    Returns (this_count, this_uniques, last_count, last_uniques).
    """
    this_count = this_uniques = last_count = last_uniques = 0
    for item in items:
        ts = datetime.fromisoformat(item["timestamp"].replace("Z", "+00:00"))
        if ts >= WEEK_START:
            this_count += item.get("count", 0)
            this_uniques += item.get("uniques", 0)
        else:
            last_count += item.get("count", 0)
            last_uniques += item.get("uniques", 0)
    return this_count, this_uniques, last_count, last_uniques


def collect_github() -> dict:
    result = {
        "stars": 0,
        "forks": 0,
        "open_issues": 0,
        "open_prs": 0,
        "views_this_week": 0,
        "unique_visitors_this_week": 0,
        "views_last_week": 0,
        "unique_visitors_last_week": 0,
        "clones_this_week": 0,
        "unique_cloners_this_week": 0,
        "clones_last_week": 0,
        "unique_cloners_last_week": 0,
        "top_referrers": [],  # [{referrer, count, uniques}]
        "top_paths": [],  # [{path, count, uniques}]
        "recent_commits": [],  # [{sha, message, author, date}]
        "new_datasets_this_week": [],  # list of new dataset directory names
        "error": None,
    }

    if not GH_TOKEN:
        result["error"] = "GH_TOKEN not set"
        return result

    try:
        # Repository metadata
        resp = requests.get(
            f"https://api.github.com/repos/{GITHUB_REPO}",
            headers=GH_HEADERS,
            timeout=20,
        )
        resp.raise_for_status()
        repo = resp.json()
        result["stars"] = repo.get("stargazers_count", 0)
        result["forks"] = repo.get("forks_count", 0)
        result["open_issues"] = repo.get("open_issues_count", 0)

        # Open PRs count
        pr_resp = _gh_get("pulls", {"state": "open", "per_page": 100})
        if pr_resp.ok:
            prs = pr_resp.json()
            result["open_prs"] = len(prs)
            result["open_issues"] = max(0, result["open_issues"] - result["open_prs"])

        # Traffic: views
        views_resp = _gh_get("traffic/views")
        if views_resp.ok:
            items = views_resp.json().get("views", [])
            (
                result["views_this_week"],
                result["unique_visitors_this_week"],
                result["views_last_week"],
                result["unique_visitors_last_week"],
            ) = _split_traffic_by_week(items)
        else:
            result["traffic_error"] = (
                f"traffic/views {views_resp.status_code}: {views_resp.text[:200]}"
            )

        # Traffic: clones
        clones_resp = _gh_get("traffic/clones")
        if clones_resp.ok:
            items = clones_resp.json().get("clones", [])
            (
                result["clones_this_week"],
                result["unique_cloners_this_week"],
                result["clones_last_week"],
                result["unique_cloners_last_week"],
            ) = _split_traffic_by_week(items)

        # Top referrers
        ref_resp = _gh_get("traffic/popular/referrers")
        if ref_resp.ok:
            result["top_referrers"] = ref_resp.json()[:8]

        # Top paths
        paths_resp = _gh_get("traffic/popular/paths")
        if paths_resp.ok:
            result["top_paths"] = paths_resp.json()[:8]

        # Recent commits
        commits_resp = _gh_get(
            "commits",
            {
                "sha": "main",
                "since": WEEK_START.strftime("%Y-%m-%dT%H:%M:%SZ"),
                "per_page": 10,
            },
        )
        if commits_resp.ok:
            for c in commits_resp.json():
                msg = c.get("commit", {}).get("message", "").split("\n")[0]
                author = c.get("commit", {}).get("author", {})
                result["recent_commits"].append(
                    {
                        "sha": c["sha"][:7],
                        "message": msg[:100],
                        "author": author.get("name", ""),
                        "date": author.get("date", "")[:10],
                    }
                )

        # New datasets added this week
        # Get commits touching data/ directory, collect distinct new subdirs
        data_commits_resp = _gh_get(
            "commits",
            {
                "sha": "main",
                "since": WEEK_START.strftime("%Y-%m-%dT%H:%M:%SZ"),
                "path": "data",
                "per_page": 20,
            },
        )
        if data_commits_resp.ok:
            new_dirs = set()
            for c in data_commits_resp.json():
                detail_resp = requests.get(c["url"], headers=GH_HEADERS, timeout=20)
                if not detail_resp.ok:
                    continue
                for f in detail_resp.json().get("files", []):
                    fname = f.get("filename", "")
                    if fname.startswith("data/") and f.get("status") == "added":
                        parts = fname.split("/")
                        if len(parts) >= 2 and parts[1]:
                            new_dirs.add(parts[1])
            result["new_datasets_this_week"] = sorted(new_dirs)

    except Exception as exc:
        result["error"] = str(exc)

    return result


# ---------------------------------------------------------------------------
# 3. PyPI download statistics
# ---------------------------------------------------------------------------


def collect_pypi() -> dict:
    result = {
        "last_week": 0,
        "last_month": 0,
        "error": None,
    }

    headers = {"User-Agent": f"shedding-hub-weekly-report/1.0 ({GITHUB_REPO})"}

    for attempt in range(3):
        try:
            resp = requests.get(
                f"https://pypistats.org/api/packages/{PYPI_PACKAGE}/recent",
                headers=headers,
                timeout=20,
            )
            if resp.status_code == 429:
                retry_after = int(resp.headers.get("Retry-After", 60))
                import time

                time.sleep(retry_after)
                continue
            resp.raise_for_status()
            data = resp.json().get("data", {})
            result["last_week"] = data.get("last_week", 0)
            result["last_month"] = data.get("last_month", 0)
            break
        except Exception as exc:
            result["error"] = str(exc)

    return result


# ---------------------------------------------------------------------------
# 4. Claude narrative
# ---------------------------------------------------------------------------


def summarize_with_claude(ga: dict, gh: dict, pypi: dict) -> str:
    if not ANTHROPIC_API_KEY:
        return "Claude summarization unavailable (no API key)."

    raw = {
        "week": WEEK_LABEL,
        "website_ga4": ga,
        "github_repository": gh,
        "pypi_downloads": pypi,
    }

    client = anthropic.Anthropic(api_key=ANTHROPIC_API_KEY)
    msg = client.messages.create(
        model="claude-sonnet-4-6",
        max_tokens=1500,
        messages=[
            {
                "role": "user",
                "content": f"""You are the data analyst for Shedding Hub, an open-source scientific data repository
for wastewater-based epidemiology research. It provides standardized biomarker shedding datasets
and statistical models used by public health researchers worldwide.

Below is last week's ({WEEK_LABEL}) raw metrics data. Write a concise weekly briefing with:

1. **Executive Summary** (2-3 sentences): overall reach and engagement trends
2. **Website Highlights**: top pages, traffic sources, geographic reach, engagement quality
3. **Repository Health**: GitHub stars/forks trends, repo traffic, new datasets added, recent development activity
4. **Package Adoption**: PyPI download trends and what they suggest about library usage
5. **Action Items**: 2-3 concrete recommendations for the coming week

Be specific with numbers. Note week-over-week changes where available (views/clones this week vs last week).
Keep it under 400 words. Use markdown headers.

Raw data:
{json.dumps(raw, indent=2, default=str)}""",
            }
        ],
    )
    return msg.content[0].text


# ---------------------------------------------------------------------------
# 5. HTML helpers
# ---------------------------------------------------------------------------


def _mini_chart(daily: dict, label: str) -> str:
    if not daily:
        return ""
    sorted_days = sorted(daily.items())
    max_val = max(v for _, v in sorted_days) or 1
    bars = ""
    for date, val in sorted_days:
        pct = int(val / max_val * 60)
        short = date[5:]
        bars += (
            f"<tr>"
            f"<td style='font-size:11px;color:#666;padding:1px 6px 1px 0;white-space:nowrap;'>{short}</td>"
            f"<td style='padding:1px 0;'>"
            f"<div style='background:#1a6b8a;height:10px;width:{pct}px;border-radius:2px;display:inline-block;'></div>"
            f"</td>"
            f"<td style='font-size:11px;color:#444;padding:1px 0 1px 6px;'>{val}</td>"
            f"</tr>"
        )
    return f"""
    <h4 style='margin:12px 0 4px;font-size:13px;color:#555;'>{label}</h4>
    <table style='border-collapse:collapse;'>{bars}</table>
    """


def _section(title: str, rows: list, err: str = "") -> str:
    err_html = f" <span style='color:red;font-size:12px;'>({err})</span>" if err else ""
    cells = "".join(
        f"<tr><td style='padding:3px 12px 3px 0;color:#555;font-size:13px;'>{k}</td>"
        f"<td style='padding:3px 0;font-weight:600;font-size:13px;'>{v}</td></tr>"
        for k, v in rows
    )
    return (
        f"<h3 style='margin:20px 0 6px;color:#1a6b8a;font-size:15px;"
        f"border-bottom:1px solid #e0e0e0;padding-bottom:4px;'>"
        f"{title}{err_html}</h3>"
        f"<table style='border-collapse:collapse;width:100%;'>{cells}</table>"
    )


def _fmt_seconds(s: float) -> str:
    if s < 60:
        return f"{int(s)}s"
    return f"{int(s // 60)}m {int(s % 60)}s"


# ---------------------------------------------------------------------------
# 6. HTML assembly
# ---------------------------------------------------------------------------


def build_html(summary: str, ga: dict, gh: dict, pypi: dict) -> str:
    # GA4 section
    ga_rows = [
        ("Active users", ga["active_users"]),
        ("New users", ga["new_users"]),
        ("Page views", ga["page_views"]),
        ("Avg engagement / user", _fmt_seconds(ga["avg_engagement_seconds"])),
    ]

    # Per-page breakdown table
    page_table = ""
    if ga["page_breakdown"]:
        header = (
            "<tr style='background:#f0f7fa;'>"
            "<th style='padding:4px 10px 4px 0;text-align:left;font-size:12px;color:#1a6b8a;'>Page</th>"
            "<th style='padding:4px 6px;text-align:right;font-size:12px;color:#1a6b8a;'>Views</th>"
            "<th style='padding:4px 6px;text-align:right;font-size:12px;color:#1a6b8a;'>Users</th>"
            "<th style='padding:4px 0 4px 6px;text-align:right;font-size:12px;color:#1a6b8a;'>Avg time</th>"
            "</tr>"
        )
        body = ""
        for p in ga["page_breakdown"]:
            body += (
                f"<tr>"
                f"<td style='padding:3px 10px 3px 0;font-size:12px;color:#333;'>{p['page']}</td>"
                f"<td style='padding:3px 6px;text-align:right;font-size:12px;'>{p['views']}</td>"
                f"<td style='padding:3px 6px;text-align:right;font-size:12px;'>{p['active_users']}</td>"
                f"<td style='padding:3px 0 3px 6px;text-align:right;font-size:12px;'>{_fmt_seconds(p['avg_engagement_s'])}</td>"
                f"</tr>"
            )
        page_table = (
            f"<h4 style='margin:14px 0 4px;font-size:13px;color:#555;'>Page breakdown</h4>"
            f"<table style='border-collapse:collapse;width:100%;'>{header}{body}</table>"
        )

    # Traffic sources
    if ga["traffic_sources"]:
        total_sess = sum(s["sessions"] for s in ga["traffic_sources"])
        for s in ga["traffic_sources"]:
            pct = round(s["sessions"] / total_sess * 100) if total_sess else 0
            ga_rows.append((f"  ↳ {s['source']}", f"{s['sessions']} sessions ({pct}%)"))

    # Top countries
    if ga["top_countries"]:
        ga_rows.append(("Geographic reach", ""))
        for c in ga["top_countries"][:5]:
            ga_rows.append((f"  ↳ {c['country']}", f"{c['active_users']} users"))

    # Device types
    if ga["device_types"]:
        device_str = " / ".join(
            f"{d['device']}: {d['active_users']}" for d in ga["device_types"]
        )
        ga_rows.append(("Devices", device_str))

    # Outbound clicks
    if ga["outbound_clicks"]:
        ga_rows.append(("Outbound clicks", ""))
        for oc in ga["outbound_clicks"][:5]:
            short_url = oc["url"][:60] + ("…" if len(oc["url"]) > 60 else "")
            ga_rows.append((f"  ↳ {short_url}", str(oc["clicks"])))

    daily_chart = _mini_chart(ga.get("daily_users", {}), "Daily active users")

    # GitHub section
    views_delta = gh["views_this_week"] - gh["views_last_week"]
    delta_sign = "+" if views_delta >= 0 else ""
    traffic_err = gh.get("traffic_error", "")
    views_value = (
        f"unavailable ({traffic_err})"
        if traffic_err
        else f"{gh['views_this_week']} ({delta_sign}{views_delta} vs last week)"
    )
    gh_rows = [
        ("Stars", gh["stars"]),
        ("Forks", gh["forks"]),
        ("Open issues", gh["open_issues"]),
        ("Open PRs", gh["open_prs"]),
        ("Repo views this week", views_value),
        (
            "Unique visitors this week",
            "unavailable" if traffic_err else gh["unique_visitors_this_week"],
        ),
        (
            "Repo clones this week",
            "unavailable" if traffic_err else gh["clones_this_week"],
        ),
        (
            "Unique cloners this week",
            "unavailable" if traffic_err else gh["unique_cloners_this_week"],
        ),
    ]
    if gh["new_datasets_this_week"]:
        gh_rows.append(("New datasets added", ", ".join(gh["new_datasets_this_week"])))
    else:
        gh_rows.append(("New datasets added", "0"))

    if gh["top_referrers"]:
        gh_rows.append(("Top referrers", ""))
        for r in gh["top_referrers"][:5]:
            gh_rows.append(
                (
                    f"  ↳ {r.get('referrer', '')}",
                    f"{r.get('count', 0)} views / {r.get('uniques', 0)} unique",
                )
            )

    commits_html = ""
    if gh["recent_commits"]:
        commit_rows = "".join(
            f"<tr>"
            f"<td style='padding:2px 8px 2px 0;font-size:11px;color:#888;font-family:monospace;'>{c['sha']}</td>"
            f"<td style='padding:2px 0;font-size:12px;color:#333;'>{c['message']}</td>"
            f"<td style='padding:2px 0 2px 8px;font-size:11px;color:#888;white-space:nowrap;'>{c['date']}</td>"
            f"</tr>"
            for c in gh["recent_commits"]
        )
        commits_html = (
            f"<h4 style='margin:14px 0 4px;font-size:13px;color:#555;'>Recent commits</h4>"
            f"<table style='border-collapse:collapse;width:100%;'>{commit_rows}</table>"
        )

    # PyPI section
    pypi_rows = [
        ("Downloads last week", f"{pypi['last_week']:,}"),
        ("Downloads last month", f"{pypi['last_month']:,}"),
    ]
    if pypi["last_month"] > 0:
        weekly_avg = round(pypi["last_month"] / 4.33)
        pypi_rows.append(("Monthly weekly avg", f"{weekly_avg:,}"))

    # Claude narrative → HTML
    summary_html = (
        summary.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )
    summary_html = re.sub(
        r"^#{1,3} (.+)$", r"<strong>\1</strong>", summary_html, flags=re.MULTILINE
    )
    summary_html = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", summary_html)
    summary_html = summary_html.replace("\n", "<br>")

    generated_at = _now.strftime("%Y-%m-%d %H:%M UTC")

    return f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1.0"></head>
<body style="margin:0;padding:0;background:#f4f4f4;font-family:Arial,Helvetica,sans-serif;">
<table width="100%" cellpadding="0" cellspacing="0" style="background:#f4f4f4;">
<tr><td align="center" style="padding:24px 16px;">
<table width="680" cellpadding="0" cellspacing="0"
  style="background:#fff;border-radius:8px;overflow:hidden;box-shadow:0 2px 8px rgba(0,0,0,.08);">

  <!-- Header -->
  <tr><td style="background:#1a6b8a;padding:22px 30px;">
    <h1 style="margin:0;color:#fff;font-size:22px;">Shedding Hub Weekly Report</h1>
    <p style="margin:5px 0 0;color:#b3d9e8;font-size:13px;">{WEEK_LABEL}</p>
  </td></tr>

  <!-- Claude narrative -->
  <tr><td style="padding:22px 30px 10px;">
    <div style="background:#f0f7fa;border-left:4px solid #1a6b8a;padding:16px 18px;border-radius:4px;
                font-size:14px;line-height:1.7;color:#333;">
      {summary_html}
    </div>
  </td></tr>

  <!-- Data sections -->
  <tr><td style="padding:6px 30px 28px;">

    {_section("Website Traffic (GA4)", ga_rows, ga.get("error", ""))}
    {page_table}
    {daily_chart}

    {_section("GitHub Repository", gh_rows, gh.get("error", ""))}
    {commits_html}

    {_section("PyPI Package Downloads", pypi_rows, pypi.get("error", ""))}

  </td></tr>

  <!-- Footer -->
  <tr><td style="background:#f0f0f0;padding:14px 30px;border-top:1px solid #e0e0e0;text-align:center;">
    <p style="margin:0;color:#888;font-size:12px;">
      Shedding Hub automated weekly report &mdash; generated {generated_at}
    </p>
  </td></tr>

</table></td></tr></table>
</body></html>"""


# ---------------------------------------------------------------------------
# 7. Send email
# ---------------------------------------------------------------------------


def send_report(html: str) -> None:
    recipients = [r.strip() for r in REPORT_TO_EMAIL.split(";") if r.strip()]
    if not recipients:
        raise ValueError("REPORT_TO_EMAIL is not set")

    subject = f"Shedding Hub Weekly Report — {WEEK_LABEL}"
    msg = MIMEMultipart("alternative")
    msg["From"] = SMTP_USER
    msg["To"] = ", ".join(recipients)
    msg["Subject"] = subject
    msg.attach(MIMEText(html, "html", "utf-8"))

    with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
        server.starttls()
        server.login(SMTP_USER, SMTP_PASSWORD)
        server.sendmail(SMTP_USER, recipients, msg.as_string())


# ---------------------------------------------------------------------------
# 8. Persist metrics for long-term tracking
# ---------------------------------------------------------------------------

METRICS_FILE = Path(__file__).parent.parent / "metrics" / "weekly_metrics.jsonl"


def save_metrics(ga: dict, gh: dict, pypi: dict) -> None:
    """Append a weekly summary record to metrics/weekly_metrics.jsonl.

    One JSON object per line; appends without modifying prior records.
    Skips writing if a record with the same week_start already exists.
    """
    record = {
        "week_start": START_LABEL,
        "week_end": END_LABEL,
        "collected_at": _now.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "github": {
            "stars": gh["stars"],
            "forks": gh["forks"],
            "open_issues": gh["open_issues"],
            "open_prs": gh["open_prs"],
            "views_this_week": gh["views_this_week"],
            "unique_visitors_this_week": gh["unique_visitors_this_week"],
            "clones_this_week": gh["clones_this_week"],
            "unique_cloners_this_week": gh["unique_cloners_this_week"],
            "new_datasets_count": len(gh.get("new_datasets_this_week", [])),
        },
        "pypi": {
            "last_week": pypi["last_week"],
            "last_month": pypi["last_month"],
        },
        "ga4": {
            "active_users": ga["active_users"],
            "new_users": ga["new_users"],
            "page_views": ga["page_views"],
            "avg_engagement_seconds": ga["avg_engagement_seconds"],
            "traffic_sources": ga["traffic_sources"],
            "top_countries": ga["top_countries"],
            "device_types": ga["device_types"],
            "page_breakdown": ga["page_breakdown"],
            "daily_users": ga["daily_users"],
        },
    }

    METRICS_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Avoid duplicate entries for the same week
    if METRICS_FILE.exists():
        with METRICS_FILE.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    existing = json.loads(line)
                    if existing.get("week_start") == START_LABEL:
                        print(
                            f"  → Metrics for {START_LABEL} already recorded; skipping."
                        )
                        return
                except json.JSONDecodeError:
                    pass

    with METRICS_FILE.open("a") as f:
        f.write(json.dumps(record) + "\n")
    print(f"  → Metrics saved to {METRICS_FILE}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print(f"Collecting weekly data for {WEEK_LABEL}...")

    print("  → Google Analytics GA4")
    ga = collect_ga4()

    print("  → GitHub repository metrics")
    gh = collect_github()

    print("  → PyPI download statistics")
    pypi = collect_pypi()

    print("  → Saving metrics snapshot")
    save_metrics(ga, gh, pypi)

    print("  → Claude weekly narrative")
    summary = summarize_with_claude(ga, gh, pypi)

    print("  → Building and sending email")
    html = build_html(summary, ga, gh, pypi)
    send_report(html)

    print(f"Weekly report sent to {REPORT_TO_EMAIL}.")


if __name__ == "__main__":
    main()
