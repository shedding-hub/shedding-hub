"""
Render one example figure per plotting function, for the documentation.

Run by the docs build, never committed: an image committed once goes stale
silently the moment a plot changes, which is the failure the documentation
design exists to correct. Regenerating means a plot that breaks fails the build.
"""

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import shedding_hub as sh  # noqa: E402

OUTPUT = REPO_ROOT / "docs" / "images"


def _dataset():
    return sh.load_dataset("woelfel2020virological", local=str(REPO_ROOT / "data"))


def _fit(catalog):
    return catalog.select(
        dataset_id="woelfel2020virological", analyte="stool", model="gamma"
    )


def main() -> int:
    import numpy as np

    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = _dataset()
    catalog = sh.load_shedding_catalog()
    fit = _fit(catalog)
    source = sh.shedding_for("SARS-CoV-2", "stool", catalog=catalog)
    traj = sh.simulate_shedding(
        source, n_individuals=200, times=np.arange(1, 31), seed=42
    )

    figures = {
        "plot_time_course": lambda: sh.plot_time_course(data, specimen="sputum"),
        "plot_time_courses": lambda: sh.plot_time_courses([data], specimen="sputum"),
        "plot_shedding_heatmap": lambda: sh.plot_shedding_heatmap(
            data, specimen="sputum", value="concentration"
        ),
        "plot_mean_trajectory": lambda: sh.plot_mean_trajectory(
            data, specimen="sputum", value="concentration"
        ),
        "plot_clearance_curve": lambda: sh.plot_clearance_curve(
            data, specimen="sputum"
        ),
        "plot_detection_probability": lambda: sh.plot_detection_probability(
            data, specimen="sputum"
        ),
        "plot_value_distribution_by_time": lambda: sh.plot_value_distribution_by_time(
            data, specimen="sputum"
        ),
        # These four take the DataFrame the matching calc_* returns, not a
        # dataset. calc_shedding_durations/peaks take dataset *ids* and read
        # from GitHub, so they are driven from the local dataset instead.
        "plot_shedding_duration": lambda: sh.plot_shedding_duration(
            sh.calc_shedding_duration(data)
        ),
        # The plural variants need the *summary* frame -- they read
        # shedding_duration_mean / n_participant, which output='individual'
        # does not produce.
        "plot_shedding_durations": lambda: sh.plot_shedding_durations(
            sh.calc_shedding_duration(data, output="summary")
        ),
        "plot_shedding_peak": lambda: sh.plot_shedding_peak(
            sh.calc_shedding_peak(data)
        ),
        "plot_shedding_peaks": lambda: sh.plot_shedding_peaks(
            sh.calc_shedding_peak(data, output="summary")
        ),
        "plot_fit_diagnostic": lambda: sh.plot_fit_diagnostic(fit, data),
        # The variant the website's dataset pages carry: a stricter
        # extrapolation gate, the full range of a simulated cohort rather than
        # its central 90%, and the axis held to the readings the curve was
        # fitted to. Refitted here rather than taken from the shipped catalog,
        # which is built at the default gate of 3.
        "plot_fit_diagnostic_range": lambda: sh.plot_fit_diagnostic(
            sh.fit_shedding_model(
                data, analyte="stool", model="gamma", max_peak_above_observed=2
            ),
            data,
            band_quantiles=(0.0, 1.0),
            band_inner_quantiles=(0.025, 0.975),
            band_sets_ylim=True,
            x_from_fitted=True,
        ),
        # A cycle-threshold fit, which no catalog ships: the axis carries real
        # Ct numbers and the height is reported as a Ct rather than a log10.
        # wang2020fecal is small enough to refit during a docs build.
        "plot_fit_diagnostic_ct": lambda: sh.plot_fit_diagnostic(
            sh.fit_shedding_model(
                sh.load_dataset("wang2020fecal", local=str(REPO_ROOT / "data")),
                analyte="stool_SARSCoV2_N",
                model="gamma",
            ),
            sh.load_dataset("wang2020fecal", local=str(REPO_ROOT / "data")),
        ),
        # A fitted analyte, drawn without its fit: the reference page is about
        # the layout, and using an unfittable analyte here would need a second
        # dataset loaded for one figure.
        "plot_analyte_observations": lambda: sh.plot_analyte_observations(
            data, "stool"
        ),
        # Scoped to one biomarker and specimen: the whole catalog renders a
        # panel per fit, which came out at 1.4 MB against <=129 KB for every
        # other figure, and grows with every dataset added.
        "plot_catalog_fits": lambda: sh.plot_catalog_fits(
            catalog, biomarker="SARS-CoV-2", specimen="stool"
        ),
        "plot_simulated_shedding": lambda: sh.plot_simulated_shedding(
            traj, source=source
        ),
    }

    for name, build in figures.items():
        # No try/except: a plot that cannot be drawn must fail the docs build
        # rather than leave the previous image in place.
        figure = build()
        figure.savefig(OUTPUT / f"{name}.png", dpi=110, bbox_inches="tight")
        print(f"wrote {name}.png")

    print(f"{len(figures)} figure(s) in {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
