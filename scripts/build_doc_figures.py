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
        "plot_catalog_fits": lambda: sh.plot_catalog_fits(catalog),
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
