import shedding_hub as sh
import pandas as pd
from typing import List, Dict, Any
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import traceback
import logging

# Constants
DEFAULT_BIOMARKER = "SARS-CoV-2"
DEFAULT_FIGURE_SIZE = (8, 6)
DEFAULT_MULTI_FIGURE_SIZE = (10, 8)
DEFAULT_MARKERSIZE = 10
NEGATIVE_VALUE = "negative"

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def calc_shedding_duration(
    dataset: Dict[str, Any],
    *,
    plotting: bool = False,
) -> pd.DataFrame:
    """
    Calculate summary statistics for the shedding duration using a loaded dataset by the `load_dataset` function.

    Args:
        dataset: Loaded dataset using the `load_dataset` function.
        plotting: Create a plot for individual level of shedding duration by specimen type.

    Returns:
        Summary table of shedding duration (min, max, mean) by biomarker and specimen.

    Raises:
        ValueError: If dataset is missing required keys or is empty.
        KeyError: If required analyte information is missing.
    """
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "study_ID"]
    missing_keys = [key for key in required_keys if key not in dataset]
    if missing_keys:
        raise ValueError(f"Dataset missing required keys: {missing_keys}")

    # extract analyte data from the standardized shedding data loaded
    df_analyte = pd.DataFrame(
        [
            {
                "analyte": key,
                "specimen": analyte["specimen"],
                "biomarker": analyte["biomarker"],
                "reference_event": analyte["reference_event"],
            }
            for key, analyte in dataset["analytes"].items()
        ]
    )

    # extract participant and measurement data from the standardized shedding data loaded
    shedding_duration_data = []
    for participant_counter, item in enumerate(dataset["participants"], 1):
        participant_ID = participant_counter
        # participant_age = item["attributes"]["age"]
        # participant_sex = item["attributes"]["sex"]
        for name, group in pd.DataFrame.from_dict(item["measurements"]).groupby(
            "analyte"
        ):
            # filter the group by value notna and time not unknown
            group = group[group["value"].notna() & (group["time"] != "unknown")].copy()

            # Skip if no valid data after filtering
            if group.empty:
                continue

            # format time to numeric
            group["time"] = pd.to_numeric(group["time"], errors="coerce")

            # Skip if no valid numeric times
            if group["time"].isna().all():
                continue

            # Calculate detection times, handling cases where no positive values exist
            positive_values = group[group["value"] != NEGATIVE_VALUE]
            first_detect = (
                positive_values["time"].min() if not positive_values.empty else pd.NA
            )
            last_detect = (
                positive_values["time"].max() if not positive_values.empty else pd.NA
            )

            row_new = {
                "study_ID": dataset["study_ID"],
                "participant_ID": participant_ID,
                # "age" : participant_age,
                # "sex" : participant_sex,
                "analyte": name,
                "n_sample": group["value"].count(),
                "first_sample": group["time"].min(),
                "last_sample": group["time"].max(),
                "first_detect": first_detect,
                "last_detect": last_detect,
            }
            shedding_duration_data.append(row_new)

    # create dataframe from list
    df_shedding_duration = pd.DataFrame(shedding_duration_data)

    # Return empty DataFrame if no data
    if df_shedding_duration.empty:
        logger.warning(
            f"No valid shedding duration data found for study {dataset['study_ID']}"
        )
        return pd.DataFrame()

    # merge analyte information and drop analyte column
    df_shedding_duration = df_shedding_duration.merge(
        df_analyte, how="left", on="analyte"
    ).drop(columns=["analyte"])

    # concatenate list of specimen types to string;
    df_shedding_duration["specimen"] = df_shedding_duration["specimen"].apply(
        lambda x: ", ".join(map(str, x)) if isinstance(x, list) else str(x)
    )

    # calculate individual level shedding duration
    # Only calculate for rows where both first_detect and last_detect are not NA
    mask = (
        df_shedding_duration["first_detect"].notna()
        & df_shedding_duration["last_detect"].notna()
    )
    df_shedding_duration.loc[mask, "shedding_duration"] = (
        df_shedding_duration.loc[mask, "last_detect"]
        - df_shedding_duration.loc[mask, "first_detect"]
        + 1
    )

    if plotting:
        plt_shedding = plot_shedding_duration(
            df_shedding_duration, dataset_ID=dataset["study_ID"]
        )
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


def plot_shedding_duration(
    df_shedding_duration: pd.DataFrame, dataset_ID: str
) -> plt.Figure:
    """
    Plot shedding duration for each individual by specimen type.

    Args:
        df_shedding_duration: Shedding duration dataset extracted from the loaded dataset.
        dataset_ID: Dataset identifier, e.g., :code:`woelfel2020virological`.

    Returns:
        The plot of shedding duration.

    Raises:
        ValueError: If DataFrame is empty or missing required columns.
    """
    if df_shedding_duration.empty:
        raise ValueError("DataFrame is empty, cannot create plot")

    required_columns = [
        "specimen",
        "first_sample",
        "last_sample",
        "first_detect",
        "last_detect",
        "reference_event",
    ]
    missing_columns = [
        col for col in required_columns if col not in df_shedding_duration.columns
    ]
    if missing_columns:
        raise ValueError(f"DataFrame missing required columns: {missing_columns}")

    # Plot range bars
    plt.figure(figsize=DEFAULT_FIGURE_SIZE)

    # Sort dataset by specimen group
    df_shedding_duration_sorted = df_shedding_duration.sort_values("specimen")

    specimen_counter = 0
    participant_counter = 0
    color_map = list(mcolors.TABLEAU_COLORS.values())
    specimen_colors = {}

    for name, group in df_shedding_duration_sorted.groupby("specimen"):
        color = color_map[specimen_counter % len(color_map)]
        specimen_colors[name] = color
        for _, row in group.iterrows():
            # Only plot if we have valid data
            if pd.notna(row["first_sample"]) and pd.notna(row["last_sample"]):
                plt.plot(
                    [row["first_sample"], row["last_sample"]],
                    [participant_counter, participant_counter],
                    linestyle="--",
                    marker="o",
                    color=color,
                )

            if pd.notna(row["first_detect"]) and pd.notna(row["last_detect"]):
                plt.plot(
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
    plt.xlabel(f"Days after {df_shedding_duration['reference_event'].iloc[0]}")
    plt.title(f'Individual Shedding Duration for the Study "{dataset_ID}"')
    plt.grid(True, axis="x")
    plt.tight_layout()
    return plt


def calc_shedding_durations(
    dataset_IDs: List[str],
    *,
    plotting: bool = False,
    biomarker: str = DEFAULT_BIOMARKER,
) -> pd.DataFrame:
    """
    Calculate summary statistics for the shedding duration using multiple loaded datasets.

    Args:
        dataset_IDs: A list of dataset identifiers.
        plotting: Create a plot for study level of shedding duration by specimen type.
        biomarker: Filter the data for plotting with a specific biomarker (e.g., "SARS-CoV-2").

    Returns:
        Summary table of shedding duration (min, max, mean, n_sample, n_participant) by study, biomarker, and specimen.

    Raises:
        ValueError: If dataset_IDs is empty or contains invalid entries.
    """
    if not dataset_IDs:
        raise ValueError("dataset_IDs cannot be empty")

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
    for filename in dataset_IDs:
        logger.info(f"Loading the data: {filename}")
        try:
            dataset = sh.load_dataset(dataset_ID=filename)
            df_shedding_duration = calc_shedding_duration(
                dataset=dataset, plotting=False
            )

            if df_shedding_duration.empty:
                logger.warning(f"No valid data found for {filename}")
                continue

            if df_shedding_durations.empty:
                df_shedding_durations = df_shedding_duration
            else:
                df_shedding_durations = pd.concat(
                    [df_shedding_durations, df_shedding_duration],
                    ignore_index=True,
                )
        except Exception as e:
            logger.error(f"Cannot load the data {filename}: {e}")
            traceback.print_exc()

    if plotting and not df_shedding_durations.empty:
        plt_sheddings = plot_shedding_durations(
            df_shedding_durations, biomarker=biomarker
        )
        plt_sheddings.show()

    return df_shedding_durations


def plot_shedding_durations(
    df_shedding_durations: pd.DataFrame, biomarker: str = DEFAULT_BIOMARKER
) -> plt.Figure:
    """
    Plot shedding duration by study and specimen type.

    Args:
        df_shedding_durations: Shedding duration dataset extracted from the multiple loaded datasets.
        biomarker: Filter the data for plotting with a specific biomarker (e.g., "SARS-CoV-2").

    Returns:
        The plot of shedding duration by study and sample type.

    Raises:
        ValueError: If DataFrame is empty or missing required columns.
    """
    if df_shedding_durations.empty:
        raise ValueError("DataFrame is empty, cannot create plot")

    required_columns = [
        "biomarker",
        "shedding_duration_mean",
        "specimen",
        "study_ID",
        "n_participant",
    ]
    missing_columns = [
        col for col in required_columns if col not in df_shedding_durations.columns
    ]
    if missing_columns:
        raise ValueError(f"DataFrame missing required columns: {missing_columns}")

    # Plot range bars
    plt.figure(figsize=DEFAULT_MULTI_FIGURE_SIZE)

    # Filter the dataset by biomarker
    df_shedding_durations_filtered = df_shedding_durations.loc[
        (
            (df_shedding_durations["biomarker"] == biomarker)
            & (df_shedding_durations["shedding_duration_mean"].notna())
        )
    ]

    if df_shedding_durations_filtered.empty:
        raise ValueError(f"No data found for biomarker: {biomarker}")

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
        for _, row in group.iterrows():
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
                markersize=DEFAULT_MARKERSIZE,
            )
            plt.plot(
                row["shedding_duration_max"],
                study_counter,
                marker="|",
                color=color,
                markersize=DEFAULT_MARKERSIZE,
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
            f"{a} (N={b})"
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
