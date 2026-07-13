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
from pathlib import Path


def _parse_jsonl(text):
    """Parse JSONL text into records sorted by ``week_start``.

    Blank lines and lines that fail to parse as JSON are skipped.
    """
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    records.sort(key=lambda r: r.get("week_start", ""))
    return records


def load_records(path):
    """Load records from a JSONL file. Missing file returns ``[]``."""
    p = Path(path)
    if not p.exists():
        return []
    return _parse_jsonl(p.read_text(encoding="utf-8"))
