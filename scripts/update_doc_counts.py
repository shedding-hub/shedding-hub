"""
Rewrite the catalog counts quoted in the documentation.

Three pages state how large the shipped catalog is, in prose: the fit count,
the study count, the analyte count, the split by model, and the number of
biomarkers and specimen types. Every one of those moves whenever the catalog is
rebuilt, and nothing had been checking them -- the figures in `index.md` and
`modeling-methods.md` are simply what someone typed after the last rebuild.

Run via `make doc_counts`, after `make catalog`.

This substitutes the numbers where they stand rather than regenerating the
sentences around them. Rewriting the prose would reflow it and bury the real
change in a whitespace diff; replacing digits in place leaves a diff that shows
exactly which counts moved. The patterns are anchored on the surrounding words
and each must match exactly once, so a page that has been reworded fails loudly
here instead of silently keeping a stale number.
"""

import argparse
import pathlib
import re
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from shedding_hub import load_shedding_catalog  # noqa: E402

DOCS = REPO_ROOT / "docs"

# `\s+` rather than a literal space throughout: every one of these sentences is
# wrapped, and which words a line break falls between depends on how wide the
# numbers happen to be. Anchoring on single spaces made the patterns miss
# whenever a rebuild changed a count's digit width.
PATTERNS = {
    "getting-started.md": [
        (
            r"(catalog of \*\*)(\d+)(\s+converged fits over\s+)(\d+)( studies\*\*)",
            ("fits", "studies"),
        ),
    ],
    "index.md": [
        (
            r"(catalog of \*\*)(\d+)(\s+converged fits over\s+)(\d+)"
            r"( studies\*\*, spanning\s+)(\d+)(\s+biomarkers and\s+)(\d+)"
            r"(\s+specimen\s+types)",
            ("fits", "studies", "biomarkers", "specimens"),
        ),
    ],
    "modeling-methods.md": [
        (
            r"(The catalog currently holds \*\*)(\d+)(\s+fits over\s+)(\d+)"
            r"(\s+studies and\s+)(\d+)(\s+analytes\*\*:\s+)(\d+)"
            r"(\s+exponential,\s+)(\d+)(\s+gamma,\s+)(\d+)"
            r"(\s+gamma_shifted, across\s+)(\d+)(\s+biomarkers)",
            (
                "fits",
                "studies",
                "analytes",
                "exponential",
                "gamma",
                "gamma_shifted",
                "biomarkers",
            ),
        ),
    ],
}


def catalog_counts(catalog=None) -> dict:
    """Every number the documentation quotes, read off the shipped catalog."""
    catalog = load_shedding_catalog() if catalog is None else catalog
    table = catalog.table
    by_model = table["model"].value_counts()
    return {
        "fits": len(table),
        "studies": table["dataset_id"].nunique(),
        "analytes": table.groupby(["dataset_id", "analyte"]).ngroups,
        "biomarkers": table["biomarker"].nunique(),
        "specimens": table["specimen"].nunique(),
        "exponential": int(by_model.get("exponential", 0)),
        "gamma": int(by_model.get("gamma", 0)),
        "gamma_shifted": int(by_model.get("gamma_shifted", 0)),
    }


def _substitute(text: str, pattern: str, names: tuple, counts: dict) -> str:
    """Replace the captured numbers, leaving every other character alone."""

    def repl(match: re.Match) -> str:
        # Odd groups are the literal text around each number, even groups are
        # the numbers themselves, in the order `names` lists them.
        parts = list(match.groups())
        for index, name in enumerate(names):
            parts[2 * index + 1] = str(counts[name])
        return "".join(parts)

    updated, n = re.subn(pattern, repl, text)
    if n != 1:
        raise SystemExit(
            f"pattern matched {n} time(s), expected exactly 1. The page has been "
            f"reworded and this script needs updating:\n  {pattern}"
        )
    return updated


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Report whether the pages are current; write nothing. Exits 1 if stale.",
    )
    args = parser.parse_args()

    counts = catalog_counts()
    print("catalog counts:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))

    stale = []
    for name, patterns in PATTERNS.items():
        path = DOCS / name
        original = path.read_text(encoding="utf-8")
        updated = original
        for pattern, names in patterns:
            updated = _substitute(updated, pattern, names, counts)
        if updated == original:
            print(f"  {name}: already current")
            continue
        stale.append(name)
        if args.check:
            print(f"  {name}: STALE")
        else:
            path.write_text(updated, encoding="utf-8")
            print(f"  {name}: updated")

    if args.check and stale:
        print(f"\nStale: {', '.join(stale)}. Run `make doc_counts`.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
