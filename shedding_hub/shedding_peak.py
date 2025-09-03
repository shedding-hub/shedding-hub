import shedding_hub as sh
import pandas as pd
from typing import List, Dict, Any, Literal
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.figure as Figure
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


def calc_shedding_peak(
    dataset: Dict[str, Any],
    *,
    output: Literal["summary", "individual"] = "summary",
) -> pd.DataFrame:
    """
    Calculate summary statistics for the shedding peak using a loaded dataset by the `load_dataset` function.

    Args:
        dataset: Loaded dataset using the `load_dataset` function.
        output: Type of dataframe returned.
            summary: Summary table of shedding peak (min, max, mean) by biomarker and specimen.
            individual: Individual entries (Not summarized). This is used for plot_shedding_peak function.

    Returns:
        Summary table of shedding peak (min, max, mean) by biomarker and specimen.

    Raises:
        ValueError: If dataset is missing required keys or is empty.
        KeyError: If required analyte information is missing.
    """
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "dataset_id"]
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
    shedding_peak_data = []
    for participant_id, item in enumerate(dataset["participants"], 1):
        for name, group in pd.DataFrame.from_dict(item["measurements"]).groupby(
            "analyte"
        ):
            # Filter numeric values for shedding peak calculation
            numeric_values = group[
                group["value"].apply(lambda x: isinstance(x, (int, float)))
            ]

            # Skip if no valid data after filtering
            if numeric_values.empty:
                continue

            # format time to numeric
            group["time"] = pd.to_numeric(group["time"], errors="coerce")

            # Skip if no valid numeric times
            if group["time"].isna().all():
                continue

            # Calculate peak time
            shedding_peak = group.loc[numeric_values["value"].idxmax(), "time"]

            row_new = {
                "dataset_id": dataset["dataset_id"],
                "participant_id": participant_id,
                "analyte": name,
                "n_sample": group["value"].count(),
                "first_sample": group[group["time"] != "unknown"]["time"].min(),
                "last_sample": group[group["time"] != "unknown"]["time"].max(),
                "shedding_peak": shedding_peak,
            }
            shedding_peak_data.append(row_new)

    # create dataframe from list
    df_shedding_peak = pd.DataFrame(shedding_peak_data)

    # Return empty DataFrame if no data
    if df_shedding_peak.empty:
        logger.warning(
            f"No valid shedding peak data found for study {dataset['dataset_id']}"
        )
        return pd.DataFrame()

    # merge analyte information and drop analyte column
    df_shedding_peak = df_shedding_peak.merge(
        df_analyte, how="left", on="analyte"
    ).drop(columns=["analyte"])

    # concatenate list of specimen types to string;
    df_shedding_peak["specimen"] = df_shedding_peak["specimen"].apply(
        lambda x: ", ".join(map(str, x)) if isinstance(x, list) else str(x)
    )

    if output == "individual":
        return df_shedding_peak

    df_shedding_peak_summary = (
        df_shedding_peak.groupby(
            ["dataset_id", "biomarker", "specimen", "reference_event"]
        )
        .agg(
            shedding_peak_min=("shedding_peak", "min"),
            shedding_peak_q25=("shedding_peak", lambda x: x.quantile(0.25)),
            shedding_peak_median=("shedding_peak", "median"),
            shedding_peak_q75=("shedding_peak", lambda x: x.quantile(0.75)),
            shedding_peak_max=("shedding_peak", "max"),
            shedding_peak_mean=("shedding_peak", "mean"),
            n_sample=("n_sample", "sum"),
            n_participant=("participant_id", "nunique"),
        )
        .reset_index()
    )

    if output == "summary":
        return df_shedding_peak_summary

    raise ValueError("`output` must be either 'summary' or 'individual'")


def plot_shedding_peak(
    df_shedding_peak: pd.DataFrame,
    *,  # Force keyword arguments for better clarity
    max_nparticipant: int = 30,
    random_seed: int = 12345,
    figsize_width_per_specimen: int = 5,
    figsize_height: int = 6,
    peak_marker: str = "D",
    peak_color: str = "tab:red",
    window_color: str = "gray",
) -> Figure:
    """
    Create a faceted error-bar plot showing each participant's shedding window and peak,
    grouped by specimen type.

    Args:
        df_shedding_peak: Individual level output from calc_shedding_peak containing
            columns: specimen, participant_id, first_sample, last_sample, shedding_peak.
        max_nparticipant: Maximum number of participants to show per specimen type.
            If exceeded, participants are randomly sampled.
        random_seed: Random seed for participant sampling when max_nparticipant is exceeded.
        figsize_width_per_specimen: Width in inches per specimen subplot.
        figsize_height: Height in inches for the entire figure.
        peak_marker: Marker style for peak points. See matplotlib marker styles.
        peak_color: Color for peak markers.
        window_color: Color for the shedding window error bars.

    Returns:
        matplotlib.figure.Figure: The generated figure containing the plot.
            Each specimen type gets its own subplot showing participant shedding windows
            and peaks.

    Note:
        - Participants are labeled as P1, P2, etc. on the y-axis
        - Error bars show the full shedding window from first to last sample
        - Diamond markers show the peak shedding time point
    """
    if df_shedding_peak.empty:
        raise ValueError("DataFrame is empty, cannot create plot")

    # Verify required columns are present for plotting
    required_columns = {
        "specimen",
        "participant_id",
        "first_sample",
        "last_sample",
        "shedding_peak",
    }
    missing_columns = set(required_columns) - set(df_shedding_peak.columns)

    if missing_columns:
        raise ValueError(
            f"DataFrame missing required columns for plotting shedding peak: {', '.join(sorted(missing_columns))}"
        )

    # Drop any rows with NA values in shedding peak
    df = df_shedding_peak.dropna(subset=["shedding_peak"]).copy()

    if df.empty:
        raise ValueError("No valid data remains after dropping NA values.")

    # limit rows per specimen for legibility
    dfs = []
    for specimen, group in df.groupby("specimen"):
        if len(group) > max_nparticipant:
            group = group.sample(n=max_nparticipant, random_state=random_seed)
        dfs.append(group)
    df = pd.concat(dfs, ignore_index=True)

    # error-bar extents
    df["err_plus"] = df["last_sample"] - df["first_sample"]
    df["err_minus"] = 0

    specimens = df["specimen"].unique()
    n = len(specimens)

    # Create figure and axes
    fig, axes = plt.subplots(
        1, n, figsize=(figsize_width_per_specimen * n, figsize_height), sharey=False
    )
    if n == 1:
        axes = [axes]

    for ax, spec in zip(axes, specimens):
        g = (
            df[df["specimen"] == spec]
            .sort_values(["shedding_peak", "participant_id"], ascending=[True, True])
            .reset_index(drop=True)  # reset index after sorting
        ) #Sort by shedding peak day.
        y = list(range(len(g)))[::-1]

        # horizontal error bars (shedding window)
        ax.errorbar(
            x=g["first_sample"],
            y=y,
            xerr=[g["err_minus"], g["err_plus"]],
            fmt="none",
            ecolor=window_color,
            elinewidth=2,
            capsize=0,
            zorder=1,
        )

        # peak markers
        ax.scatter(
            g["shedding_peak"],
            y,
            marker=peak_marker,
            s=35,
            color=peak_color,
            zorder=2,
            label="Peak",
        )

        ax.set_yticks([])
        ax.set_xlabel(f"Days after {df_shedding_peak['reference_event'].iloc[0]}")
        ax.set_title(spec)
        ax.grid(axis="x", alpha=0.3)

    axes[0].set_ylabel("Participant")
    axes[-1].legend(loc="upper right")
    fig.suptitle(
        f"Individual Shedding Peak for Dataset '{df_shedding_peak['dataset_id'].iloc[0]}'",
        y=1.02,
        fontsize=14,
    )
    plt.tight_layout()
    return fig


def calc_shedding_peaks(
    dataset_ids: List[str],
    *,
    biomarker: str = DEFAULT_BIOMARKER,
    reference_event: str = "symptom onset",
    min_nparticipant: int = 5,
) -> pd.DataFrame:
    """
    Calculate summary statistics for the shedding peak using multiple loaded datasets.

    Args:
        dataset_ids: A list of dataset identifiers.
        biomarker: Filter the data for plotting with a specific biomarker (e.g., "SARS-CoV-2").
        reference_event: Filter the data for plotting with a specific reference event (e.g., "symptom onset").
        min_nparticipant: Filter the data for plotting with a minimum number of participants (e.g., 5).

    Returns:
        Summary table of shedding peak (min, q1, median, q3, max, mean, n_sample, n_participant) by study, biomarker, and specimen.

    Raises:
        ValueError: If dataset_ids is empty or contains invalid entries.
    """
    if not dataset_ids:
        raise ValueError("dataset_ids cannot be empty")

    loaded_datasets = []
    for dataset_id in dataset_ids:
        logger.info(f"Loading the data: {dataset_id}")
        loaded_datasets.append(
            calc_shedding_peak(
                dataset=sh.load_dataset(dataset=dataset_id), output="summary"
            )
        )

    df_shedding_peaks = pd.concat(loaded_datasets, ignore_index=True)

    return df_shedding_peaks


def plot_shedding_peaks(
    df_shedding_peaks: pd.DataFrame,
    *,
    biomarker: str = DEFAULT_BIOMARKER,
    reference_event: str = "symptom onset",
    min_nparticipant: int = 5,
    x_axis_upper_limit: int | None = 50,
) -> Figure:
    """
    Plot shedding peaks by study, biomarker, specimen type, and reference event.

    Args:
        df_shedding_peaks: Summary table of shedding peak (min, max, mean, n_sample, n_participant) by study, biomarker, and specimen.
        biomarker: Filter the data for plotting with a specific biomarker (e.g., "SARS-CoV-2").
        reference_event: Filter the data for plotting with a specific reference event (e.g., "symptom onset").
        min_nparticipant: Filter the data for plotting with a minimum number of participants (e.g., 5).
        x_axis_upper_limit: Upper limit for x-axis in days (None for auto-scaling).

    Returns:
        matplotlib.figure.Figure: Box plot showing shedding peaks across studies.
    """
    # Filter summary based on criteria
    df_filtered = df_shedding_peaks.query(
        "biomarker == @biomarker and "
        "reference_event == @reference_event and "
        "n_participant >= @min_nparticipant"
    ).copy()

    if df_filtered.empty:
        raise ValueError(
            "No rows match the chosen biomarker/reference_event "
            "with the minimum number of participants required."
        )

    df_filtered = df_filtered.sort_values(["specimen", "dataset_id"])

    # Build box-plot stats
    bxp_stats = []
    specimen_order = []
    for _, r in df_filtered.iterrows():
        bxp_stats.append(
            {
                "label": f"{r['dataset_id']} (n={r['n_participant']})",
                "whislo": r["shedding_peak_min"],
                "q1": r["shedding_peak_q25"],
                "med": r["shedding_peak_median"],
                "q3": r["shedding_peak_q75"],
                "whishi": r["shedding_peak_max"],
            }
        )
        specimen_order.append(r["specimen"])

    # Create figure and plot
    fig, ax = plt.subplots(figsize=DEFAULT_MULTI_FIGURE_SIZE)

    bp = ax.bxp(
        bxp_stats, orientation="horizontal", patch_artist=True, showfliers=False
    )

    # Style the plot
    for line in bp["medians"]:
        line.set_color("black")
        line.set_linewidth(1.5)

    # Create color map for unique specimens
    unique_specimens = pd.Series(specimen_order).unique()
    cmap = plt.colormaps["tab10"].resampled(len(unique_specimens))
    spec_colors = {spec: cmap(i) for i, spec in enumerate(unique_specimens)}

    for patch, spec in zip(bp["boxes"], specimen_order):
        patch.set_facecolor(spec_colors[spec])
        patch.set_edgecolor("black")
        patch.set_alpha(0.75)

    if x_axis_upper_limit is not None:
        ax.set_xlim(-5, x_axis_upper_limit)

    ax.set_xlabel(f"Days after {reference_event}")
    ax.set_ylabel("")
    ax.set_title(f"Shedding Peak Plot for {biomarker}")
    ax.axvline(0, color="gray", linestyle="--", linewidth=1)

    # Add legend for unique specimens
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markerfacecolor=spec_colors[s],
            markeredgecolor="black",
            label=s,
            markersize=DEFAULT_MARKERSIZE,
        )
        for s in unique_specimens
    ]
    ax.legend(
        handles,
        [h.get_label() for h in handles],
        title="Specimen",
        loc="upper right",
        bbox_to_anchor=(1, 1),
    )

    plt.tight_layout()
    return fig
