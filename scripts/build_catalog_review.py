"""
Render every fit in the shipped catalog against the data behind it.

One page per fit — the observations, the fitted median individual, and the
parameter estimates — collected into a single PDF to page through. Run via
`make review`.

The PDF is deliberately not committed: it is a regenerable binary that would be
rewritten on every catalog rebuild.
"""

import argparse
import pathlib
import sys
import warnings

import matplotlib

matplotlib.use("Agg")

from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from shedding_hub import (  # noqa: E402
    load_dataset,
    load_shedding_catalog,
    plot_fit_diagnostic,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data", default=str(REPO_ROOT / "data"), help="Directory of datasets."
    )
    parser.add_argument(
        "--catalog", default=None, help="Catalog file. Defaults to the shipped one."
    )
    parser.add_argument(
        "--output",
        default=str(REPO_ROOT / "shedding_catalog_review.pdf"),
        help="PDF to write.",
    )
    args = parser.parse_args()

    catalog = load_shedding_catalog(args.catalog)
    if not catalog.fits:
        print("catalog holds no fits, nothing to review")
        return 1

    # Sorted so a study's two models face each other, and each study's analytes
    # stay together.
    fits = sorted(
        catalog.fits, key=lambda fit: (fit.dataset_id, fit.analyte, fit.model)
    )

    datasets: dict = {}
    failures = []
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(output) as pdf:
        for number, fit in enumerate(fits, 1):
            if fit.dataset_id not in datasets:
                print(f"loading {fit.dataset_id}", flush=True)
                datasets[fit.dataset_id] = load_dataset(fit.dataset_id, local=args.data)
            print(
                f"  [{number}/{len(fits)}] {fit.dataset_id} / {fit.analyte} / "
                f"{fit.model}",
                flush=True,
            )
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    figure = plot_fit_diagnostic(fit, datasets[fit.dataset_id])
            except Exception as error:  # noqa: BLE001
                # One unrenderable fit should not cost the other 82 pages, but it
                # is reported rather than passed over: a missing page would
                # otherwise read as a fit that does not exist.
                failures.append((fit.dataset_id, fit.analyte, fit.model, str(error)))
                continue
            pdf.savefig(figure)

    print(f"\nwrote {len(fits) - len(failures)} page(s) to {output}")
    if failures:
        print(f"{len(failures)} fit(s) could not be rendered:")
        for dataset_id, analyte, model, message in failures:
            print(f"  {dataset_id} / {analyte} / {model}: {message}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
