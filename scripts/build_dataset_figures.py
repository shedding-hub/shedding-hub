"""
Render one figure per analyte for the website's dataset pages.

Every analyte in ``data/`` gets at least one figure. Where a fit exists the page
is the fit diagnostic under the gate-2 catalog with a full-range band, which
shows both the fitted curve and the spread a simulation from it would produce.
Where none exists -- 88 of the repository's 216 analytes, mostly cross-sectional
studies with one reading per participant -- the page is observations only, in the
same layout, so a dataset with partial coverage does not render half-blank.

Output is ``figures/<dataset_id>/<analyte>__<model>.png``, plus
``figures/index.json`` naming the default figure and the alternatives for each
analyte. The website's Makefile copies both out of the repository archive it
already downloads for the dataset YAML.

Run via `make figures`, after `make catalog_gate2` and `make catalog_ct_gate2`.
"""

import argparse
import json
import pathlib
import sys
import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from shedding_hub import (  # noqa: E402
    load_dataset,
    load_shedding_catalog,
    plot_analyte_observations,
    plot_fit_diagnostic,
)

# Preference order when an analyte has more than one model. gamma_shifted
# resolves the most structure -- it estimates onset rather than assuming it --
# and exponential resolves the least, so the default page is the most
# informative fit available rather than whichever sorted first.
MODEL_PREFERENCE = ("gamma_shifted", "gamma", "exponential")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default=str(REPO_ROOT / "data"))
    parser.add_argument("--output", default=str(REPO_ROOT / "figures"))
    parser.add_argument(
        "--catalog",
        default=str(REPO_ROOT / "shedding_catalog_gate2.yaml"),
        help="Concentration catalog. Defaults to the gate-2 build.",
    )
    parser.add_argument(
        "--ct-catalog",
        default=str(REPO_ROOT / "shedding_catalog_ct_gate2.yaml"),
        help="Cycle-threshold catalog. Defaults to the gate-2 build.",
    )
    # 100 dpi puts a 9-inch figure at 900px, which is about the width of the
    # site's is-max-desktop container. Below that the browser upscales and the
    # legend -- which carries every fitted parameter -- goes soft.
    parser.add_argument("--dpi", type=int, default=100)
    parser.add_argument(
        "--only", nargs="+", default=None, help="Limit to these dataset ids."
    )
    args = parser.parse_args()

    data_dir = pathlib.Path(args.data)
    output = pathlib.Path(args.output)

    fits: dict[tuple[str, str], dict[str, object]] = {}
    for path in (args.catalog, args.ct_catalog):
        if not pathlib.Path(path).is_file():
            print(f"missing catalog {path}; run the matching make target first")
            return 1
        for fit in load_shedding_catalog(path).fits:
            fits.setdefault((fit.dataset_id, fit.analyte), {})[fit.model] = fit

    dataset_ids = sorted(
        path.name
        for path in data_dir.iterdir()
        if path.is_dir() and not path.name.startswith(".")
    )
    if args.only:
        dataset_ids = [i for i in dataset_ids if i in set(args.only)]

    index: dict[str, dict] = {}
    n_fit = n_obs = 0
    failures: list[tuple[str, str, str]] = []

    for dataset_id in dataset_ids:
        dataset = load_dataset(dataset_id, local=str(data_dir))
        target = output / dataset_id
        target.mkdir(parents=True, exist_ok=True)
        entries = []

        for analyte in dataset.get("analytes") or {}:
            available = fits.get((dataset_id, analyte), {})
            models = [m for m in MODEL_PREFERENCE if m in available]
            figures = []

            if models:
                for model in models:
                    name = f"{analyte}__{model}.png"
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", UserWarning)
                            figure = plot_fit_diagnostic(
                                available[model],
                                dataset,
                                band_quantiles=(0.0, 1.0),
                                band_inner_quantiles=(0.025, 0.975),
                                band_sets_ylim=True,
                                x_from_fitted=True,
                            )
                    except Exception as error:  # noqa: BLE001
                        failures.append((dataset_id, analyte, f"{model}: {error}"))
                        continue
                    figure.savefig(target / name, dpi=args.dpi, bbox_inches="tight")
                    plt.close(figure)
                    figures.append({"model": model, "file": name})
                    n_fit += 1

            if not figures:
                # Either the analyte has no fit, or every one of its fits failed
                # to render. Both leave the page with nothing, so both fall back.
                name = f"{analyte}__observations.png"
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", UserWarning)
                        figure = plot_analyte_observations(dataset, analyte)
                except Exception as error:  # noqa: BLE001
                    failures.append((dataset_id, analyte, f"observations: {error}"))
                    continue
                figure.savefig(target / name, dpi=args.dpi, bbox_inches="tight")
                plt.close(figure)
                figures.append({"model": "observations", "file": name})
                n_obs += 1

            entries.append({"analyte": analyte, "figures": figures})

        index[dataset_id] = {"analytes": entries}
        print(f"  {dataset_id}: {sum(len(e['figures']) for e in entries)}", flush=True)

    (output / "index.json").write_text(
        json.dumps(index, indent=1, sort_keys=True), encoding="utf-8"
    )

    total = sum(len(e["figures"]) for d in index.values() for e in d["analytes"])
    print(f"\nwrote {total} figure(s): {n_fit} fit, {n_obs} observations-only")
    print(f"index at {output / 'index.json'}")
    if failures:
        print(f"{len(failures)} figure(s) could not be rendered:")
        for dataset_id, analyte, message in failures:
            print(f"  {dataset_id} / {analyte} / {message}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
