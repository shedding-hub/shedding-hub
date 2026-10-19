"""Cost and time of a benchmark run, from its logs (Analysis E, benchmark/PLAN.md).

Each study folder of a run holds `log.json`, written by the runner. From it
this file takes, per agent, the tokens the API reported, the number of API
calls and the wall-clock seconds, and prices the tokens:

    extraction (OpenAI)     prompt tokens at the input price, completion
                            tokens at the output price. The log does not
                            split cached prompt tokens, so all are priced as
                            uncached; the cost is an upper bound.
    review (Anthropic)      input, cache-write, cache-read and output tokens,
                            each at its own price.

Prices are US dollars per million tokens (`PRICES`; another table can be
given with `--prices`). The OpenAI prices are those used for the run's cost
estimate on 7 Oct 2026. The Anthropic prices are the list prices of
claude-opus-4-8 as of 6 Oct 2026, with cache writes (5-minute) at 1.25 times
and cache reads at 0.1 times the input price. Check both against the
providers' pages on the day of the run.

Per study it reports tokens, cost and time, also per 100 reference
measurements, and over studies the median and interquartile range, by era and
overall, with the run's totals. The models each response named are listed, so
that the run can be checked against the pinned identifiers.

Usage (from the repository root):
    python benchmark/cost.py RUN_DIR --comparison CMP_DIR --out OUT_DIR [--prices FILE]
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import discrepancies as disc  # noqa: E402
import metrics  # noqa: E402

PRICES = {
    "extraction": {"input": 1.75, "output": 14.0},
    "review": {"input": 5.0, "cache_write": 6.25, "cache_read": 0.5, "output": 25.0},
}

STUDY_COLUMNS = [
    "study_id",
    "era",
    "finished",
    "reference_measurements",
    "extraction_calls",
    "extraction_input_tokens",
    "extraction_output_tokens",
    "extraction_seconds",
    "extraction_usd",
    "review_calls",
    "review_input_tokens",
    "review_cache_write_tokens",
    "review_cache_read_tokens",
    "review_output_tokens",
    "review_seconds",
    "review_usd",
    "total_usd",
    "wall_clock_seconds",
    "usd_per_100_measurements",
    "seconds_per_100_measurements",
    "models",
]
SUMMARY_COLUMNS = ["group", "quantity", "median", "q1", "q3", "total", "studies"]
SUMMARIZED = [
    "extraction_input_tokens",
    "extraction_output_tokens",
    "review_input_tokens",
    "review_output_tokens",
    "extraction_usd",
    "review_usd",
    "total_usd",
    "extraction_seconds",
    "review_seconds",
    "wall_clock_seconds",
    "usd_per_100_measurements",
    "seconds_per_100_measurements",
]


def usd(tokens: dict, prices: dict) -> float:
    return sum(tokens.get(k, 0) * prices.get(k, 0.0) for k in tokens) / 1e6


def study_cost(log: dict, measurements: int, prices: dict = PRICES) -> dict:
    """Tokens, cost and time of one study, from its log."""
    extraction = log.get("extraction") or {}
    review = log.get("review") or {}
    e_usage = extraction.get("usage_totals") or {}
    r_usage = review.get("usage_totals") or {}
    e_tokens = {
        "input": int(e_usage.get("prompt_tokens") or 0),
        "output": int(e_usage.get("completion_tokens") or 0),
    }
    r_tokens = {
        "input": int(r_usage.get("input_tokens") or 0),
        "cache_write": int(r_usage.get("cache_creation_input_tokens") or 0),
        "cache_read": int(r_usage.get("cache_read_input_tokens") or 0),
        "output": int(r_usage.get("output_tokens") or 0),
    }
    e_usd = usd(e_tokens, prices["extraction"])
    r_usd = usd(r_tokens, prices["review"])
    wall = None
    if log.get("started") and isinstance(log.get("finished_at"), str):
        start = datetime.fromisoformat(log["started"])
        wall = (datetime.fromisoformat(log["finished_at"]) - start).total_seconds()
    per100 = 100 / measurements if measurements else None
    models = sorted(
        set(extraction.get("models_returned") or [])
        | set(review.get("models_returned") or [])
    )
    return {
        "finished": bool(log.get("finished")),
        "reference_measurements": measurements,
        "extraction_calls": int(extraction.get("api_calls") or 0),
        "extraction_input_tokens": e_tokens["input"],
        "extraction_output_tokens": e_tokens["output"],
        "extraction_seconds": float(extraction.get("seconds") or 0.0),
        "extraction_usd": e_usd,
        "review_calls": int(review.get("api_calls") or 0),
        "review_input_tokens": r_tokens["input"],
        "review_cache_write_tokens": r_tokens["cache_write"],
        "review_cache_read_tokens": r_tokens["cache_read"],
        "review_output_tokens": r_tokens["output"],
        "review_seconds": float(review.get("seconds") or 0.0),
        "review_usd": r_usd,
        "total_usd": e_usd + r_usd,
        "wall_clock_seconds": wall,
        "usd_per_100_measurements": (e_usd + r_usd) * per100 if per100 else None,
        "seconds_per_100_measurements": wall * per100 if per100 and wall else None,
        "models": ";".join(models),
    }


def summarize(group: str, rows: list[dict]) -> list[dict]:
    out = []
    for quantity in SUMMARIZED:
        values = [r[quantity] for r in rows if r[quantity] is not None]
        if not values:
            continue
        q1, median, q3 = np.percentile(values, [25, 50, 75])
        out.append(
            {
                "group": group,
                "quantity": quantity,
                "median": float(median),
                "q1": float(q1),
                "q3": float(q3),
                # A sum of per-100-measurement rates means nothing.
                "total": "" if "per_100" in quantity else float(sum(values)),
                "studies": len(values),
            }
        )
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir", type=Path, help="benchmark/runs/<run_id>")
    ap.add_argument("--comparison", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--prices", type=Path, help="JSON table shaped like PRICES")
    args = ap.parse_args(argv)

    prices = (
        json.loads(args.prices.read_text(encoding="utf-8")) if args.prices else PRICES
    )
    manifest = disc.read_manifest()
    measurements = {
        r["study_id"]: int(r["s2_measurements"] or 0)
        for r in metrics.read_csv(args.comparison / "comparison_summary.csv")
    }
    rows, missing = [], []
    for study_dir in sorted(d for d in args.run_dir.iterdir() if d.is_dir()):
        study_id = study_dir.name
        if study_id not in manifest:
            continue
        log = study_dir / "log.json"
        if not log.is_file():
            missing.append(study_id)
            continue
        cost = study_cost(
            json.loads(log.read_text(encoding="utf-8")),
            measurements.get(study_id, 0),
            prices,
        )
        rows.append({"study_id": study_id, "era": manifest[study_id]["era"], **cost})

    summary = summarize("all", rows)
    for era in ("manual", "ai_assisted"):
        summary += summarize(era, [r for r in rows if r["era"] == era])
    args.out.mkdir(parents=True, exist_ok=True)
    disc.write_csv(args.out / "cost_studies.csv", STUDY_COLUMNS, rows)
    disc.write_csv(args.out / "cost_summary.csv", SUMMARY_COLUMNS, summary)
    (args.out / "prices.json").write_text(
        json.dumps(prices, indent=2) + "\n", encoding="utf-8"
    )
    models = sorted({m for r in rows for m in r["models"].split(";") if m})
    total = sum(r["total_usd"] for r in rows)
    print(
        f"{len(rows)} studies, ${total:.2f} in all, models: {', '.join(models)} "
        f"-> {args.out}"
    )
    if missing:
        print(f"NO LOG {', '.join(missing)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
