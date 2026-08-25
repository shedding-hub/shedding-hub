"""Write the catalogue's growth curve for the website's curation page.

The website shows how the catalogue grew: slowly by hand for fifteen months,
then faster once the extraction agents landed. It cannot compute that itself.
The shape lives in *this* repository's commit history, and the site receives
`data/` as a zip with no git attached; its CI is Ruby only besides.

So the series is built here, committed like `figures/`, and copied out of the
archive the site already downloads.

    python scripts/build_curation_growth.py        # -> curation_growth.yaml

Output is normalized SVG coordinates rather than raw counts: the site renders
it in Liquid, which has no arithmetic worth the name, so every number the
figure needs -- axis ticks, label positions, the path itself -- is resolved
here.

REQUIRES FULL HISTORY. Under a shallow clone (`actions/checkout` defaults to
`fetch-depth: 1`) `git log` reports one commit, and the curve silently
collapses to a single point rather than failing. The check below refuses that
instead of publishing a wrong figure.
"""

import argparse
import bisect
import datetime as dt
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "curation_growth.yaml"

# viewBox the figure is drawn in, and the padding that keeps the axis labels
# and the trailing total inside it.
W, H = 1200.0, 380.0
PAD_L, PAD_R, PAD_T, PAD_B = 48.0, 118.0, 26.0, 34.0

# The first agent-extracted batch. Pinned rather than detected: which batch was
# the first the agents produced is a fact about how the work was done, not
# something recoverable from the shape of the curve.
AI_START = dt.date(2026, 2, 5)

# A repository with real history has far more than this; a shallow clone has 1.
MIN_COMMITS = 20


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(ROOT), *args],
        capture_output=True, text=True, check=True,
    ).stdout


def dataset_counts() -> list[tuple[dt.date, int]]:
    """Cumulative dataset count at every commit that touched `data/`."""
    lines = [ln for ln in _git(
        "log", "--format=%H %ad", "--date=short", "--reverse", "--", "data"
    ).split("\n") if ln.strip()]

    if len(lines) < MIN_COMMITS:
        raise SystemExit(
            f"only {len(lines)} commit(s) touch data/ -- this looks like a "
            "shallow clone. The growth curve needs full history; use "
            "`actions/checkout` with `fetch-depth: 0`."
        )

    per_day: dict[dt.date, int] = {}
    for line in lines:
        sha, date = line.split()
        listing = subprocess.run(
            ["git", "-C", str(ROOT), "ls-tree", "-d", "--name-only", f"{sha}:data"],
            capture_output=True, text=True,
        ).stdout
        n = len([x for x in listing.split("\n") if x.strip()])
        day = dt.date.fromisoformat(date)
        # A day can carry several merges; the last state that day is the one.
        per_day[day] = max(per_day.get(day, 0), n)

    return sorted(per_day.items())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=DEST)
    args = ap.parse_args()

    pts = dataset_counts()
    days = [d for d, _ in pts]
    counts = [n for _, n in pts]
    t0, t1 = days[0], days[-1]
    span = (t1 - t0).days or 1

    step = 25
    y_max = max(step * 2, ((counts[-1] // step) + 1) * step)

    def x_of(d: dt.date) -> float:
        return PAD_L + (d - t0).days / span * (W - PAD_L - PAD_R)

    def y_of(n: float) -> float:
        return H - PAD_B - (n / y_max) * (H - PAD_T - PAD_B)

    def count_on(d: dt.date) -> int:
        i = bisect.bisect_right(days, d) - 1
        return counts[i] if i >= 0 else 0

    # A step path: the catalogue holds its value until a batch lands, so
    # interpolating between merges would draw studies that did not exist yet.
    path: list[str] = []
    prev_y = None
    for d, n in pts:
        x, y = x_of(d), y_of(n)
        if prev_y is not None:
            path.append(f"{x:.1f},{prev_y:.1f}")
        path.append(f"{x:.1f},{y:.1f}")
        prev_y = y

    xticks = [
        {"x": round(x_of(dt.date(year, month, 1)), 1),
         "label": dt.date(year, month, 1).strftime("%b %Y")}
        for year in range(t0.year, t1.year + 1)
        for month in (1, 7)
        if t0 <= dt.date(year, month, 1) <= t1
    ]

    doc = {
        "as_of": t1.isoformat(),
        "width": W,
        "height": H,
        "baseline_y": round(y_of(0), 1),
        "plot_top": round(y_of(y_max), 1),
        "plot_left": round(PAD_L, 1),
        "plot_right": round(W - PAD_R, 1),
        "path": " ".join(path),
        "gridlines": [{"y": round(y_of(v), 1), "label": str(v)}
                      for v in range(0, y_max + 1, step)],
        "xticks": xticks,
        "ai_band": {"x": round(x_of(AI_START), 1),
                    "width": round(x_of(t1) - x_of(AI_START), 1)},
        "ai_start": {"x": round(x_of(AI_START), 1),
                     "y": round(y_of(count_on(AI_START)), 1),
                     "date_label": AI_START.strftime("%b %Y"),
                     "count": count_on(AI_START)},
        "end": {"x": round(x_of(t1), 1),
                "y": round(y_of(counts[-1]), 1),
                "count": counts[-1],
                "date_label": t1.strftime("%b %Y")},
        # What the catalogue held the day before the agents landed: the number
        # the page contrasts the current total against.
        "manual_era_end": count_on(AI_START - dt.timedelta(days=1)),
    }

    header = (
        "# Generated by scripts/build_curation_growth.py from this repository's\n"
        "# git history, for the website's curation page. Do not edit by hand.\n"
    )
    args.output.write_text(header + yaml.safe_dump(doc, sort_keys=False),
                           encoding="utf-8")
    print(f"wrote {args.output}")
    print(f"  {len(pts)} change-days, {counts[-1]} studies as of {t1}")
    print(f"  manual era ended at {doc['manual_era_end']}, "
          f"AI era starts {AI_START} at {doc['ai_start']['count']}")


if __name__ == "__main__":
    main()
