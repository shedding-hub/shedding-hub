import shedding_hub as sh
import pandas as pd
from typing import Optional
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import glob
from pathlib import Path


def calc_shedding_duration(
    df_ID: str,
    *,
    plotting: bool = False,
    ref: Optional[str] = None,
    pr: Optional[int] = None,
    local: Optional[str] = None,
) -> pd.DataFrame:
    """
    Calcualte summary statistics for the shedding duration using a loaded dataset by the `load_dataset` function.

    Args:
        df_ID: Dataset identifier, e.g., :code:`woelfel2020virological`.
        plotting: Create a plot for individual level of shedding duration by specimen type.
        ref: Git reference to load. Defaults to the most recent data on the :code:`main`
            branch of https://github.com/shedding-hub/shedding-hub and is automatically
            fetched if a :code:`pr` number is specified.
        pr: Pull request to fetch data from.
        local: Local directory to load data from.

    Returns:
        Summary table of shedding duration (min, max, mean) by biomarker and specimen.
    """
    # load the dataset;
    df = sh.load_dataset(dataset=df_ID, ref=ref, pr=pr, local=local)
    # initialize a pandas dataframe
    df_shedding_duration = pd.DataFrame(
        columns=[
            "study_ID",
            "participant_ID",  # "age","sex",
            "analyte",
            "n_sample",
            "first_sample",
            "last_sample",
            "first_detect",
            "last_detect",
        ]
    )
    df_analyte = pd.DataFrame(
        columns=["analyte", "specimen", "biomarker", "reference_event"]
    )

    # extract analyte data from the standardized shedding data loaded
    for key in df["analytes"]:
        row_new = {
            "analyte": key,
            "specimen": df["analytes"][key]["specimen"],
            "biomarker": df["analytes"][key]["biomarker"],
            "reference_event": df["analytes"][key]["reference_event"],
        }
        df_analyte.loc[len(df_analyte)] = row_new

    # extract participant and measurement data from the standardized shedding data loaded
    participant_counter = 0
    for item in df["participants"]:
        participant_counter += 1
        participant_ID = participant_counter
        # participant_age = item["attributes"]["age"]
        # participant_sex = item["attributes"]["sex"]
        for name, group in pd.DataFrame.from_dict(item["measurements"]).groupby(
            "analyte"
        ):
            row_new = {
                "study_ID": df_ID,
                "participant_ID": participant_ID,
                # "age" : participant_age,
                # "sex" : participant_sex,
                "analyte": name,
                "n_sample": group[group["value"].notna()]["time"].size,
                "first_sample": group[group["value"].notna()]["time"].min(),
                "last_sample": group[group["value"].notna()]["time"].max(),
                "first_detect": group[
                    (group["value"] != "negative") & (group["value"].notna())
                ]["time"].min(),
                "last_detect": group[
                    (group["value"] != "negative") & (group["value"].notna())
                ]["time"].max(),
            }
            df_shedding_duration.loc[len(df_shedding_duration)] = row_new

    # merge analyte information and drop analyte column
    df_shedding_duration = df_shedding_duration.merge(
        df_analyte, how="left", on="analyte"
    ).drop(columns=["analyte"])
    # calculate individual level shedding duration
    df_shedding_duration["first_sample"] = pd.to_numeric(
        df_shedding_duration["first_sample"], errors="coerce"
    )
    df_shedding_duration["last_sample"] = pd.to_numeric(
        df_shedding_duration["last_sample"], errors="coerce"
    )
    df_shedding_duration["first_detect"] = pd.to_numeric(
        df_shedding_duration["first_detect"], errors="coerce"
    )
    df_shedding_duration["last_detect"] = pd.to_numeric(
        df_shedding_duration["last_detect"], errors="coerce"
    )
    df_shedding_duration["shedding_duration"] = (
        df_shedding_duration["last_detect"] - df_shedding_duration["first_detect"] + 1
    )

    if plotting == True:
        plt_shedding = plot_shedding_duration(df_shedding_duration, df_ID=df_ID)
        plt_shedding.show()

    df_return = (
        df_shedding_duration.groupby(
            ["study_ID", "biomarker", "specimen", "reference_event"]
        )
        .agg(
            shedding_duration_min=("shedding_duration", "min"),
            shedding_duration_max=("shedding_duration", "max"),
            shedding_duration_mean=("shedding_duration", "mean"),
            n_sample=("n_sample", "sum"),
            n_participant=("participant_ID", "nunique"),
        )
        .reset_index()
    )

    return df_return


def plot_shedding_duration(df_shedding_duration: pd.DataFrame, df_ID: str):
    """
    Plot shedding duration for each individuals by specimen type.

    Args:
        df_shedding_duration: Shedding duration dataset extracted from the loaded dataset.
        df_ID: Dataset identifier, e.g., :code:`woelfel2020virological`.

    Returns:
        The plot of shedding duration.
    """
    # Plot range bars
    plt.figure(figsize=(8, 6))

    # Sort dataset by specimen group
    df_shedding_duration_sorted = df_shedding_duration.sort_values("specimen")

    specimen_counter = 0
    participant_counter = 0
    color_map = list(mcolors.TABLEAU_COLORS.values())
    specimen_colors = {}

    for name, group in df_shedding_duration_sorted.groupby("specimen"):
        color = color_map[specimen_counter % len(color_map)]
        specimen_colors[name] = color
        for __, row in group.iterrows():
            (line_sampl,) = plt.plot(
                [row["first_sample"], row["last_sample"]],
                [participant_counter, participant_counter],
                linestyle="--",
                marker="o",
                color=color,
            )
            (line_shed,) = plt.plot(
                [row["first_detect"], row["last_detect"]],
                [participant_counter, participant_counter],
                linestyle="-",
                marker="o",
                color=color,
            )
            participant_counter += 1
        specimen_counter += 1

    # Legend for specimen (color)
    specimen_legend = [
        Line2D([0], [0], color=color, lw=2, label=specimen)
        for specimen, color in specimen_colors.items()
    ]

    # Legend for line style (duration type)
    linestyle_legend = [
        Line2D(
            [0], [0], color="black", lw=2, linestyle="--", label="Sampling Duration"
        ),
        Line2D([0], [0], color="black", lw=2, linestyle="-", label="Shedding Duration"),
    ]

    # Add both legends
    plt.legend(
        handles=specimen_legend + linestyle_legend, title="Legend", loc="upper right"
    )

    plt.yticks([])
    plt.xlabel(f"Days after {df_shedding_duration['reference_event'][0]}")
    plt.title(f'Individual Shedding Duration for the Study "{df_ID}"')
    plt.grid(True, axis="x")
    plt.tight_layout()
    return plt


def calc_shedding_durations(
    df_IDs: list[str], *, plotting: bool = False, biomarker: str = "SARS-CoV-2"
) -> pd.DataFrame:
    """
    Calcualte summary statistics for the shedding duration using a loaded dataset by the `load_dataset` function.

    Args:
        df_IDs : A list of dataset identifiers.
        plotting: Create a plot for study level of shedding duration by specimen type.
        biomarker: Filter the data for plotting with a specific biomarker (e.g., "SARS-CoV-2").

    Returns:
        Summary table of shedding duration (min, max, mean, n_sample, n_participant) by study, biomarker, and specimen.
    """
    # Initialize an empty DataFrame with columns
    df_shedding_durations = pd.DataFrame(
        columns=[
            "study_ID",
            "biomarker",
            "specimen",
            "reference_event",
            "shedding_duration_min",
            "shedding_duration_max",
            "shedding_duration_mean",
            "n_sample",
            "n_participant",
        ]
    )

    # Loop to append data
    for (
        filename
    ) in df_IDs:  # [Path(file).stem for file in Path("data").glob("*/*.yaml")]:
        print(f"Load the data: {filename}")
        try:
            if len(df_shedding_durations) == 0:
                df_shedding_durations = calc_shedding_duration(
                    df_ID=filename, plotting=False
                )

            elif len(calc_shedding_duration(df_ID=filename, plotting=False)) == 0:
                pass
            else:
                df_shedding_durations = pd.concat(
                    [
                        df_shedding_durations,
                        calc_shedding_duration(df_ID=filename, plotting=False),
                    ],
                    ignore_index=True,
                )
        except:
            print(f"Cannot load the data {filename}!!!")
    if plotting == True:
        plt_sheddings = plot_shedding_durations(
            df_shedding_durations, biomarker=biomarker
        )
        plt_sheddings.show()
    return df_shedding_durations


def plot_shedding_durations(
    df_shedding_durations: pd.DataFrame, biomarker: str = "SARS-CoV-2"
):
    """
    Plot shedding duration by study and specimen type.

    Args:
        df_shedding_durations: Shedding duration dataset extracted from the multiple loaded datasets.
        biomarker: Filter the data for plotting with a specific biomarker (e.g., "SARS-CoV-2").

    Returns:
        The plot of shedding duration by study and sample type.
    """
    # Plot range bars
    plt.figure(figsize=(10, 8))

    # filter the dataset by biomaker
    df_shedding_durations_filtered = df_shedding_durations.loc[
        (
            (df_shedding_durations["biomarker"] == biomarker)
            & (~df_shedding_durations["shedding_duration_mean"].isnull())
        )
    ]

    # Sort dataset by specimen group
    df_shedding_durations_sorted = df_shedding_durations_filtered.sort_values(
        "specimen"
    )

    specimen_counter = 0
    study_counter = 0
    color_map = list(mcolors.TABLEAU_COLORS.values())

    # Store a handle for each specimen group for the legend
    legend_handles = []

    for name, group in df_shedding_durations_sorted.groupby("specimen"):
        color = color_map[specimen_counter % len(color_map)]
        for __, row in group.iterrows():
            x_vals = [
                row["shedding_duration_min"],
                row["shedding_duration_mean"],
                row["shedding_duration_max"],
            ]
            y_vals = [study_counter] * 3
            # Plot the connecting line without markers
            (line,) = plt.plot(x_vals, y_vals, color=color, linestyle="-")
            # Add custom markers
            plt.plot(
                row["shedding_duration_min"],
                study_counter,
                marker="|",
                color=color,
                markersize=10,
            )
            plt.plot(
                row["shedding_duration_max"],
                study_counter,
                marker="|",
                color=color,
                markersize=10,
            )
            plt.plot(
                row["shedding_duration_mean"], study_counter, marker="o", color=color
            )
            study_counter += 1
        legend_handles.append((line, name))
        specimen_counter += 1

    # Add legend
    handles, labels = zip(*legend_handles)
    plt.legend(handles, labels, title="Specimen", loc="upper right")

    plt.yticks(
        ticks=range(len(df_shedding_durations_sorted["study_ID"])),
        labels=[
            a + " (N=" + str(b) + ")"
            for a, b in zip(
                df_shedding_durations_sorted["study_ID"].values,
                df_shedding_durations_sorted["n_participant"],
            )
        ],
    )
    plt.xlabel("Number of Days")
    plt.title(f"Shedding Duration Plot for {biomarker}")
    plt.grid(True, axis="x")
    plt.tight_layout()
    return plt
