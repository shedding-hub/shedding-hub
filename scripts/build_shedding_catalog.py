"""
Regenerate the precomputed shedding-model catalog shipped with the package.

Run via `make catalog`. Fitting every analyte of every dataset takes a while,
which is exactly why the result is precomputed rather than fitted on demand.
"""

import argparse
import pathlib
import sys
import warnings

import yaml

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from shedding_hub import load_dataset  # noqa: E402
from shedding_hub.shedding_catalog import (  # noqa: E402
    CATALOG_PATH,
    fit_shedding_models,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data", default=str(REPO_ROOT / "data"), help="Directory of datasets."
    )
    parser.add_argument(
        "--output", default=str(CATALOG_PATH), help="Catalog file to write."
    )
    args = parser.parse_args()

    data_dir = pathlib.Path(args.data)
    dataset_ids = sorted(
        path.name
        for path in data_dir.iterdir()
        if path.is_dir() and not path.name.startswith(".")
    )

    datasets = []
    for dataset_id in dataset_ids:
        print(f"loading {dataset_id}", flush=True)
        datasets.append(load_dataset(dataset_id, local=str(data_dir)))

    print(f"fitting {len(datasets)} dataset(s)", flush=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        catalog = fit_shedding_models(datasets)

    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(catalog.to_dict(), stream, sort_keys=False)

    print(f"wrote {len(catalog.fits)} fit(s) to {output}")
    print(f"skipped {len(catalog.skipped)} analyte/model combination(s)")
    if not catalog.skipped.empty:
        print(catalog.skipped["reason"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
