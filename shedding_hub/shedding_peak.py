import shedding_hub as sh
import pandas as pd
from typing import List, Dict, Any, Literal
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
from matplotlib.figure import Figure
import logging
import numpy as np

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
    output: Literal["summary", "individual"] = "individual",
) -> pd.DataFrame:
    """
    Calculate summary statistics for the shedding peak using a loaded dataset by the `load_dataset` function.

    Args:
        dataset: Loaded dataset using the `load_dataset` function.
        output: Type of dataframe returned.
            individual: Individual entries (Not summarized). This is used for plot_shedding_peak function.
            summary: Summary table of shedding peak (min, max, mean) by biomarker and specimen.

    Returns:
        DataFrame of shedding peak either individual or summary.

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

    # extract analyte data from the standardized shedding data loaded (include unit)
    df_analyte = pd.DataFrame(
        [
            {
                "analyte": key,
                "specimen": analyte["specimen"],
                "biomarker": analyte["biomarker"],
                "reference_event": analyte["reference_event"],
                "unit": analyte["unit"],
            }
            for key, analyte in dataset["analytes"].items()
        ]
    )

    # map analyte -> unit for per-analyte decision when computing peak
    analyte_unit_map = {
        row["analyte"]: row.get("unit") for _, row in df_analyte.iterrows()
    }

    # extract participant and measurement data from the standardized shedding data loaded
    shedding_peak_data = []
    for participant_id, item in enumerate(dataset["participants"], 1):
        meas = pd.DataFrame.from_dict(item["measurements"])
        if meas.empty:
            continue

        for name, group in meas.groupby("analyte"):
            # coerce values to numeric where possible (handle numeric strings); drop non-numeric rows
            g = group.copy()
            g["value_num"] = pd.to_numeric(g["value"], errors="coerce")
            g = g.dropna(subset=["value_num"])
            if g.empty:
                continue

            # coerce time to numeric for selection and for first/last calculations
            g["time_num"] = pd.to_numeric(g["time"], errors="coerce")
            # if all times are NA, skip
            if g["time_num"].isna().all():
                continue

            # decide whether to pick min (cycle threshold) or max (other units)
            unit = analyte_unit_map.get(name)
            pick_min = (
                isinstance(unit, str) and unit.strip().lower() == "cycle threshold"
            )

            # select index of interest using rows with numeric values only
            if pick_min:
                sel_idx = g["value_num"].idxmin()
            else:
                sel_idx = g["value_num"].idxmax()

            # get shedding_peak time and value from the selected row (use numeric-coerced time)
            shedding_peak_time = g.loc[sel_idx, "time_num"]
            shedding_peak_value = g.loc[sel_idx, "value_num"]
            if pd.isna(shedding_peak_time):
                # skip if selected peak time is NA
                continue

            # compute first/last sample using numeric times (exclude non-numeric / 'unknown')
            times_num = pd.to_numeric(group["time"], errors="coerce")
            valid_times = times_num.loc[times_num.notna()]

            # keep NA (np.nan) if there is no numeric time entry for this group
            first_sample = valid_times.min() if not valid_times.empty else np.nan
            last_sample = valid_times.max() if not valid_times.empty else np.nan

            row_new = {
                "dataset_id": dataset["dataset_id"],
                "participant_id": participant_id,
                "analyte": name,
                "n_sample": group["value"].count(),
                "first_sample": first_sample,
                "last_sample": last_sample,
                "shedding_peak": shedding_peak_time,
                "shedding_peak_value": shedding_peak_value,
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
            ["dataset_id", "biomarker", "specimen", "reference_event", "unit"]
        )
        .agg(
            shedding_peak_min=("shedding_peak", "min"),
            shedding_peak_q25=("shedding_peak", lambda x: x.quantile(0.25)),
            shedding_peak_median=("shedding_peak", "median"),
            shedding_peak_q75=("shedding_peak", lambda x: x.quantile(0.75)),
            shedding_peak_max=("shedding_peak", "max"),
            shedding_peak_mean=("shedding_peak", "mean"),
            shedding_peak_value_min=("shedding_peak_value", "min"),
            shedding_peak_value_q25=("shedding_peak_value", lambda x: x.quantile(0.25)),
            shedding_peak_value_median=("shedding_peak_value", "median"),
            shedding_peak_value_q75=("shedding_peak_value", lambda x: x.quantile(0.75)),
            shedding_peak_value_max=("shedding_peak_value", "max"),
            shedding_peak_value_mean=("shedding_peak_value", "mean"),
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
        )  # Sort by shedding peak day.
        y = list(range(len(g)))[::-1]

        # horizontal dashed lines (shedding window)
        for i, (_, row) in enumerate(g.iterrows()):
            ax.plot(
                [row["first_sample"], row["last_sample"]],
                [y[i], y[i]],
                linestyle="--",
                color=window_color,
                linewidth=2,
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
    # Close the figure in the pyplot state to avoid Jupyter/matplotlib auto-display
    # while still returning the Figure object for the caller to display once.
    plt.close(fig)
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

    bp = ax.bxp(bxp_stats, vert=False, patch_artist=True, showfliers=False)

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
    # Close the figure in the pyplot state to avoid duplicate display in notebooks.
    plt.close(fig)
    return fig


def plot_shedding_peak_values(
    df_shedding_peaks: pd.DataFrame,
    *,
    biomarker: str = DEFAULT_BIOMARKER,
    reference_event: str = "symptom onset",
    min_nparticipant: int = 5,
) -> Figure:
    """
    Plot shedding peak values by study, biomarker, specimen type, and reference event.

    Args:
        df_shedding_peaks: Summary table of shedding peak values (min, q25, median, q75, max, mean)
            by study, biomarker, specimen, and unit.
        biomarker: Filter the data for plotting with a specific biomarker (e.g., "SARS-CoV-2").
        reference_event: Filter the data for plotting with a specific reference event (e.g., "symptom onset").
        min_nparticipant: Filter the data for plotting with a minimum number of participants (e.g., 5).

    Returns:
        matplotlib.figure.Figure: Box plot showing shedding peak values across studies.

    Note:
        - Values are automatically log10-transformed UNLESS the unit contains 'cycle threshold'
        - Different units are plotted in separate subplots
        - X-axis label will show "log10(Peak value)" when transformed
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

    # Get unique units
    unique_units = df_filtered["unit"].unique()
    n_units = len(unique_units)

    # Create figure with subplots for each unit
    fig, axes = plt.subplots(
        1,
        n_units,
        figsize=(DEFAULT_MULTI_FIGURE_SIZE[0] * n_units, DEFAULT_MULTI_FIGURE_SIZE[1]),
        squeeze=False,
    )
    axes = axes.flatten()

    # Plot for each unit
    for unit_idx, unit in enumerate(unique_units):
        ax = axes[unit_idx]

        # Filter data for this unit
        df_unit = df_filtered[df_filtered["unit"] == unit].copy()
        df_unit = df_unit.sort_values(["specimen", "dataset_id"])

        # Determine if log transformation is needed
        is_cycle_threshold = "cycle threshold" in unit.lower()

        # Build box-plot stats
        bxp_stats = []
        specimen_order = []
        for _, r in df_unit.iterrows():
            # Get values and apply transformation if needed
            whislo = r["shedding_peak_value_min"]
            q1 = r["shedding_peak_value_q25"]
            med = r["shedding_peak_value_median"]
            q3 = r["shedding_peak_value_q75"]
            whishi = r["shedding_peak_value_max"]

            if not is_cycle_threshold:
                # Apply log transformation
                whislo = np.log10(whislo)
                q1 = np.log10(q1)
                med = np.log10(med)
                q3 = np.log10(q3)
                whishi = np.log10(whishi)

            bxp_stats.append(
                {
                    "label": f"{r['dataset_id']} (n={r['n_participant']})",
                    "whislo": whislo,
                    "q1": q1,
                    "med": med,
                    "q3": q3,
                    "whishi": whishi,
                }
            )
            specimen_order.append(r["specimen"])

        # Create boxplot
        bp = ax.bxp(bxp_stats, vert=False, patch_artist=True, showfliers=False)

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

        # Set labels based on transformation
        if is_cycle_threshold:
            ax.set_xlabel(f"Peak value ({unit})")
        else:
            ax.set_xlabel(f"log10(Peak value) ({unit})")

        ax.set_ylabel("")
        ax.set_title(f"Unit: {unit}")

        # Add legend for unique specimens (only on first subplot)
        if unit_idx == 0:
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

    # Overall title
    fig.suptitle(
        f"Shedding Peak Values for {biomarker}",
        fontsize=14,
        y=1.00 if n_units == 1 else 1.02,
    )

    plt.tight_layout()
    # Close the figure in the pyplot state to avoid duplicate display in notebooks.
    plt.close(fig)
    return fig


def plot_shedding_peak_value(
    df_shedding_peak: pd.DataFrame,
    *,
    facet_by_specimen: bool = True,
    figsize_width_per_facet: int = 6,
    figsize_height: int = 6,
    marker_size: int = 50,
    alpha: float = 0.6,
    add_trendline: bool = True,
) -> Figure:
    """
    Create an XY scatter plot showing the relationship between shedding peak time (X)
    and shedding peak value (Y). If multiple units are present, they are plotted separately.

    Args:
        df_shedding_peak: Individual level output from calc_shedding_peak containing
            columns: shedding_peak, shedding_peak_value, biomarker, specimen, reference_event, unit.
        facet_by_specimen: If True, create separate subplots for each specimen type (within each unit).
            If False, plot all data on a single plot with biomarker color-coding.
            Note: If multiple units exist, they are always plotted separately regardless of this setting.
        figsize_width_per_facet: Width in inches per facet subplot.
        figsize_height: Height in inches for the entire figure.
        marker_size: Size of scatter plot markers.
        alpha: Transparency of markers (0-1, where 1 is fully opaque).
        add_trendline: If True, add a linear regression trendline for each biomarker.

    Returns:
        matplotlib.figure.Figure: The generated figure containing the scatter plot(s).

    Note:
        - Different biomarkers are shown with different colors
        - Different units are shown in separate subplots
        - Y-axis labels include the unit of measurement
        - Trendlines show the linear relationship between peak time and peak value
        - X-axis: Days after reference event when peak occurs
        - Y-axis: Test value/result at the peak time
        - Log transformation: Values are automatically log10-transformed UNLESS the unit contains
          'cycle threshold'. The Y-axis label will show "log10(Peak value)" when transformed.
    """
    if df_shedding_peak.empty:
        raise ValueError("DataFrame is empty, cannot create plot")

    # Verify required columns are present
    required_columns = {
        "shedding_peak",
        "shedding_peak_value",
        "biomarker",
        "specimen",
        "reference_event",
        "unit",
    }
    missing_columns = required_columns - set(df_shedding_peak.columns)

    if missing_columns:
        raise ValueError(
            f"DataFrame missing required columns for plotting: {', '.join(sorted(missing_columns))}"
        )

    # Drop rows with NA values in key columns
    df = df_shedding_peak.dropna(subset=["shedding_peak", "shedding_peak_value"]).copy()

    if df.empty:
        raise ValueError("No valid data remains after dropping NA values.")

    # Get unique values for faceting
    unique_units = df["unit"].unique()
    unique_biomarkers = df["biomarker"].unique()

    # Create color map for biomarkers
    cmap = plt.colormaps["tab10"].resampled(len(unique_biomarkers))
    biomarker_colors = {
        biomarker: cmap(i) for i, biomarker in enumerate(unique_biomarkers)
    }

    # Determine faceting strategy
    # Priority: Unit > Specimen
    # If multiple units exist, always facet by unit first
    if len(unique_units) > 1:
        # Multiple units: create subplot for each unit
        facet_groups = [(unit, df[df["unit"] == unit]) for unit in unique_units]
        facet_labels = [f"Unit: {unit}" for unit in unique_units]
    elif facet_by_specimen:
        # Single unit, facet by specimen
        unique_specimens = df["specimen"].unique()
        facet_groups = [(spec, df[df["specimen"] == spec]) for spec in unique_specimens]
        facet_labels = [str(spec) for spec in unique_specimens]
    else:
        # Single plot
        facet_groups = [(None, df)]
        facet_labels = ["Shedding Peak Time vs Peak Value"]

    # Create figure
    n_facets = len(facet_groups)
    fig, axes = plt.subplots(
        1,
        n_facets,
        figsize=(figsize_width_per_facet * n_facets, figsize_height),
        squeeze=False,
    )
    axes = axes.flatten()

    # Plot for each facet
    for idx, (facet_key, df_facet) in enumerate(facet_groups):
        ax = axes[idx]
        ax.set_title(facet_labels[idx])

        if df_facet.empty:
            ax.text(
                0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes
            )
            continue

        # Get unit for Y-axis label
        unit = df_facet["unit"].iloc[0] if "unit" in df_facet.columns else "Unknown"

        # Determine if log transformation is needed
        is_cycle_threshold = "cycle threshold" in unit.lower()

        # Plot each biomarker
        for biomarker in unique_biomarkers:
            df_biomarker = df_facet[df_facet["biomarker"] == biomarker]

            if df_biomarker.empty:
                continue

            # Get x and y values
            x_values = df_biomarker["shedding_peak"].values
            y_values = df_biomarker["shedding_peak_value"].values

            # Apply log transformation if not cycle threshold
            if not is_cycle_threshold:
                y_values = np.log10(y_values)

            # Scatter plot
            ax.scatter(
                x_values,
                y_values,
                c=[biomarker_colors[biomarker]],
                s=marker_size,
                alpha=alpha,
                label=biomarker,
                edgecolors="black",
                linewidth=0.5,
            )

            # Add trendline
            if add_trendline and len(df_biomarker) > 1:
                # Calculate correlation coefficient
                correlation = np.corrcoef(x_values, y_values)[0, 1]

                # Fit line
                z = np.polyfit(x_values, y_values, 1)
                p = np.poly1d(z)

                # Plot trendline
                x_line = np.linspace(x_values.min(), x_values.max(), 100)
                ax.plot(
                    x_line,
                    p(x_line),
                    color=biomarker_colors[biomarker],
                    linestyle="--",
                    linewidth=1.5,
                    alpha=0.8,
                    label=f"{biomarker} trend (r={correlation:.2f})",
                )

        # Set labels with unit information
        reference_event = df_facet["reference_event"].iloc[0]
        ax.set_xlabel(f"Days after {reference_event}")

        # Set Y-axis label based on transformation
        if is_cycle_threshold:
            ax.set_ylabel(f"Peak value ({unit})")
        else:
            ax.set_ylabel(f"log10(Peak value) ({unit})")
        ax.grid(True, alpha=0.3)

        # Add legend
        if idx == n_facets - 1:  # Only on last subplot
            ax.legend(loc="best", framealpha=0.9)

    # Overall title
    dataset_id = df["dataset_id"].iloc[0] if "dataset_id" in df.columns else "Unknown"
    fig.suptitle(
        f"Shedding Peak Time vs Peak Value for Dataset '{dataset_id}'",
        y=1.02 if n_facets > 1 else 1.00,
        fontsize=14,
    )

    plt.tight_layout()
    plt.close(fig)
    return fig
