"""
Write the fitted parameter table for every dataset, in JSON and CSV.

Reads the shipped catalog; fits nothing. Run via `make parameters` after
`make catalog`.

The JSON is the reusable form: each record carries the population mean and
covariance, the measurement-error SD and the censoring limit, so an estimate
can be simulated from without this package and without refitting. The CSV is
the flat browsing form, one row per fit, and drops the covariance -- a k x k
matrix per fit does not belong in a spreadsheet column.
"""

import argparse
import json
import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from shedding_hub import load_shedding_catalog  # noqa: E402
from shedding_hub.shedding_export import catalog_to_records  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog", default=None, help="Catalog file. Defaults to the shipped one."
    )
    parser.add_argument(
        "--json",
        default=str(REPO_ROOT / "docs" / "shedding_parameters.json"),
        help="JSON file to write.",
    )
    parser.add_argument(
        "--csv",
        default=str(REPO_ROOT / "docs" / "shedding_parameters.csv"),
        help="CSV file to write. Pass an empty string to skip it.",
    )
    args = parser.parse_args()

    catalog = load_shedding_catalog(args.catalog)
    records = catalog_to_records(catalog)
    if not records:
        print("catalog holds no fits, nothing to export")
        return 1

    payload = {
        "n_fits": len(records),
        "n_datasets": len({record["dataset_id"] for record in records}),
        "models": sorted({record["model"] for record in records}),
        "fits": records,
    }
    json_path = pathlib.Path(args.json)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with json_path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=False)
        stream.write("\n")
    print(f"wrote {len(records)} record(s) to {json_path}")

    if args.csv:
        csv_path = pathlib.Path(args.csv)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        catalog.table.to_csv(csv_path, index=False)
        print(f"wrote {len(catalog.table)} row(s) to {csv_path}")

    print(f"{payload['n_datasets']} dataset(s), models {payload['models']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
