import pandas as pd
import numpy as np
from typing import List, Dict, Any
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
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


def _is_ct_value(unit: str) -> bool:
    """
    Determine if a unit represents cycle threshold (CT) values.

    Args:
        unit: Unit string from analyte metadata.

    Returns:
        True if the unit represents CT values, False otherwise.
    """
    if not unit:
        return False
    unit_lower = str(unit).lower()
    return "cycle threshold" in unit_lower or unit_lower == "ct"


def plot_time_course(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    value: str | None = None,
    max_nparticipant: int = 30,
    random_seed: int = 12345,
    figsize_width_per_specimen: int = 5,
    figsize_height: int = 6,
    line_color: str | None = None,
    marker: str = "o",
    markersize: int = DEFAULT_MARKERSIZE,
    line_alpha: float = 0.7,
    show_negative: bool = False,
) -> Figure:
    """
    Plot individual participant shedding trajectories over time for a single dataset.

    Creates faceted line plots showing biomarker measurements over time for each participant,
    organized by specimen type. Each subplot represents a different specimen type, with
    individual participant trajectories shown as lines with markers.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, plots all biomarkers.
        specimen: Optional filter for a specific specimen type. If None, plots all specimens.
        value: Optional filter for value type. Options are "concentration" or "ct".
            If None, plots all data (may raise error if mixed). Use this to avoid
            mixed CT/concentration errors. Defaults to None.
        max_nparticipant: Maximum number of participants to plot per specimen type.
            If exceeded, randomly samples participants. Defaults to 30.
        random_seed: Random seed for reproducible sampling. Defaults to 12345.
        figsize_width_per_specimen: Width of each specimen subplot in inches. Defaults to 5.
        figsize_height: Height of the figure in inches. Defaults to 6.
        line_color: Color for all lines. If None, uses default color cycle per participant.
        marker: Marker style for data points. Defaults to "o".
        markersize: Size of markers. Defaults to DEFAULT_MARKERSIZE (10).
        line_alpha: Transparency of lines (0-1). Defaults to 0.7.
        show_negative: If True, plots "negative" values at y=0. If False, excludes them.
            Defaults to False.

    Returns:
        matplotlib.figure.Figure: The generated figure containing the time course plots.
            Each subplot shows participant trajectories for one specimen type.

    Raises:
        ValueError: If dataset is missing required keys, is empty, or has no valid data.
        KeyError: If required analyte information is missing.

    Note:
        - Time values of "unknown" are excluded from plots
        - Non-numeric values (except "negative"/"positive") are excluded
        - If more than max_nparticipant per specimen, randomly samples participants
        - Axis labels are derived from the reference_event in analyte metadata
    """
    # Validate input
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "dataset_id"]
    missing_keys = [key for key in required_keys if key not in dataset]
    if missing_keys:
        raise ValueError(f"Dataset missing required keys: {missing_keys}")

    if not dataset["participants"]:
        raise ValueError("Dataset has no participants")

    # Extract time series data from raw dataset
    time_series_data = []
    for participant_id, participant in enumerate(dataset["participants"], 1):
        measurements = participant.get("measurements", [])
        for measurement in measurements:
            analyte_name = measurement.get("analyte")
            time_series_data.append(
                {
                    "participant_id": participant_id,
                    "time": measurement.get("time"),
                    "value": measurement.get("value"),
                    "analyte": analyte_name,
                }
            )

    if not time_series_data:
        raise ValueError("Dataset has no measurements")

    df = pd.DataFrame(time_series_data)

    # Convert time to numeric, filtering out "unknown" values
    df = df[df["time"] != "unknown"].copy()
    df["time_num"] = pd.to_numeric(df["time"], errors="coerce")

    # Join with analyte metadata to get specimen and unit information
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        # Handle both single specimen and array of specimens
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        # Get limit of detection (could be numeric or "unknown")
        lod = analyte_info.get("limit_of_detection")
        lod_numeric = None
        if lod is not None and lod != "unknown":
            try:
                lod_numeric = float(lod)
            except (ValueError, TypeError):
                lod_numeric = None

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "unit": analyte_info.get("unit"),
            "reference_event": analyte_info.get("reference_event"),
            "biomarker": analyte_info.get("biomarker"),
            "limit_of_detection": lod_numeric,
        }

    # Add specimen and metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
    )
    df["unit"] = df["analyte"].map(lambda x: analyte_metadata.get(x, {}).get("unit"))
    df["reference_event"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("reference_event")
    )
    df["biomarker"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("biomarker")
    )
    df["limit_of_detection"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("limit_of_detection")
    )

    # Filter by biomarker if specified
    if biomarker is not None:
        df = df[df["biomarker"] == biomarker]
        if df.empty:
            raise ValueError(f"No measurements found for biomarker '{biomarker}'")

    # Filter by specimen if specified
    if specimen is not None:
        df = df[df["specimen"] == specimen]
        if df.empty:
            raise ValueError(f"No measurements found for specimen '{specimen}'")

    # Determine if each row is CT value or concentration based on unit
    df["is_ct"] = df["unit"].apply(_is_ct_value)

    # Filter by value type if specified
    if value is not None:
        value_lower = value.lower()
        if value_lower == "concentration":
            df = df[~df["is_ct"]].copy()
        elif value_lower == "ct":
            df = df[df["is_ct"]].copy()
        else:
            raise ValueError(
                f"Invalid value '{value}'. " "Must be 'concentration', 'ct', or None."
            )

        if df.empty:
            raise ValueError(f"No {value} data found in dataset after filtering")

    # Handle negative values based on show_negative parameter and data type
    # Track substitution values for caption
    negative_substitution_values = set()

    if show_negative:
        # Substitute negative values: use LOD if available, otherwise 45 for CT or 1 for concentrations
        def substitute_negative(row):
            if row["value"] == NEGATIVE_VALUE:
                # Use limit of detection if available
                if row["limit_of_detection"] is not None and not pd.isna(
                    row["limit_of_detection"]
                ):
                    sub_val = row["limit_of_detection"]
                    negative_substitution_values.add(("LOD", sub_val))
                    return sub_val
                else:
                    # Fallback to default values
                    sub_val = 45 if row["is_ct"] else 1
                    negative_substitution_values.add(("default", sub_val, row["is_ct"]))
                    return sub_val
            return row["value"]

        df["value_num"] = df.apply(substitute_negative, axis=1)
        df["value_num"] = pd.to_numeric(df["value_num"], errors="coerce")
    else:
        # Exclude negative values
        df = df[df["value"] != NEGATIVE_VALUE].copy()
        df["value_num"] = pd.to_numeric(df["value"], errors="coerce")

    # Drop rows with NaN time or value
    df = df.dropna(subset=["time_num", "value_num"])

    if df.empty:
        raise ValueError("No valid numeric measurements found after filtering")

    # Drop rows without specimen info
    df = df.dropna(subset=["specimen"])

    if df.empty:
        raise ValueError("No measurements with valid specimen information found")

    # Sample participants if needed (per specimen)
    dfs = []
    for specimen, group in df.groupby("specimen"):
        unique_participants = group["participant_id"].unique()
        if len(unique_participants) > max_nparticipant:
            sampled_participants = (
                pd.Series(unique_participants)
                .sample(n=max_nparticipant, random_state=random_seed)
                .values
            )
            group = group[group["participant_id"].isin(sampled_participants)]
        dfs.append(group)
    df = pd.concat(dfs, ignore_index=True)

    # Create faceted plot
    specimens = df["specimen"].unique()
    n_specimens = len(specimens)

    fig, axes = plt.subplots(
        1,
        n_specimens,
        figsize=(figsize_width_per_specimen * n_specimens, figsize_height),
        squeeze=False,
    )
    axes = axes.flatten()

    # Get reference event and unit for labeling (use first non-null)
    reference_event = (
        df["reference_event"].dropna().iloc[0]
        if not df["reference_event"].dropna().empty
        else "reference"
    )
    unit = df["unit"].dropna().iloc[0] if not df["unit"].dropna().empty else ""

    # Check if we have mixed CT and concentration data
    if df["is_ct"].nunique() > 1:
        raise ValueError(
            "Dataset contains mixed CT values and concentrations. "
            "Cannot plot on same axes. Use data_type='concentration' or data_type='ct' to filter."
        )

    # Determine plot type based on data
    is_ct = df["is_ct"].iloc[0] if not df.empty else False

    # Compute global axis ranges for consistency across panels
    x_min, x_max = df["time_num"].min(), df["time_num"].max()
    y_min, y_max = df["value_num"].min(), df["value_num"].max()

    # Add some padding to ranges
    x_padding = (x_max - x_min) * 0.05 if x_max > x_min else 1
    x_range = (x_min - x_padding, x_max + x_padding)

    if is_ct:
        # For CT values, add padding but keep linear scale
        y_padding = (y_max - y_min) * 0.05 if y_max > y_min else 1
        y_range = (y_min - y_padding, y_max + y_padding)
    else:
        # For concentrations on log scale, add padding in log space
        if y_min > 0:
            log_y_min, log_y_max = np.log10(y_min), np.log10(y_max)
            log_padding = (
                (log_y_max - log_y_min) * 0.1 if log_y_max > log_y_min else 0.5
            )
            y_range = (10 ** (log_y_min - log_padding), 10 ** (log_y_max + log_padding))
        else:
            y_range = (y_min, y_max)

    # Plot each specimen in its own subplot
    for idx, specimen in enumerate(specimens):
        ax = axes[idx]
        specimen_df = df[df["specimen"] == specimen]

        # Get biomarker for this specimen (first non-null value)
        biomarker = (
            specimen_df["biomarker"].dropna().iloc[0]
            if not specimen_df["biomarker"].dropna().empty
            else ""
        )

        # Plot each participant's trajectory
        for participant_id in specimen_df["participant_id"].unique():
            participant_data = specimen_df[
                specimen_df["participant_id"] == participant_id
            ].sort_values("time_num")
            color = line_color if line_color else None
            ax.plot(
                participant_data["time_num"],
                participant_data["value_num"],
                marker=marker,
                markersize=markersize,
                alpha=line_alpha,
                color=color,
                linewidth=1.5,
            )

        # Apply axis scaling and range
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)

        if is_ct:
            # For CT values, reverse y-axis (smaller on top, larger on bottom)
            ax.invert_yaxis()
        else:
            # For concentrations, use log scale
            ax.set_yscale("log")

        # Styling
        ax.set_xlabel(f"Time after {reference_event} (days)")
        ax.set_ylabel(f"Measurement ({unit})" if unit else "Measurement")

        # Format specimen name (replace underscores with spaces) and include biomarker
        specimen_display = specimen.replace("_", " ")
        title = f"{biomarker} - {specimen_display}" if biomarker else specimen_display
        ax.set_title(title)
        ax.grid(axis="x", alpha=0.3)

    # Remove extra subplots if any
    for idx in range(n_specimens, len(axes)):
        fig.delaxes(axes[idx])

    # Add main title
    dataset_id = dataset.get("dataset_id", "Dataset")
    fig.suptitle(f"Time Course: {dataset_id}", fontsize=14, y=0.98)

    # Add caption for negative value substitution if applicable
    if show_negative and negative_substitution_values:
        caption_parts = []
        # Check for LOD substitutions
        lod_values = [
            item[1] for item in negative_substitution_values if item[0] == "LOD"
        ]
        if lod_values:
            lod_str = ", ".join([f"{v:.2g}" for v in sorted(set(lod_values))])
            caption_parts.append(f"limit of detection ({lod_str})")

        # Check for default substitutions
        default_items = [
            item for item in negative_substitution_values if item[0] == "default"
        ]
        if default_items:
            # Group by CT vs concentration
            ct_vals = [item[1] for item in default_items if item[2]]
            conc_vals = [item[1] for item in default_items if not item[2]]
            if ct_vals:
                caption_parts.append(f"CT value of {ct_vals[0]}")
            if conc_vals:
                caption_parts.append(f"concentration of {conc_vals[0]}")

        if caption_parts:
            caption = "Note: Negative values were substituted with " + " or ".join(
                caption_parts
            )
            fig.text(
                0.5, -0.05, caption, ha="center", fontsize=9, style="italic", wrap=True
            )

    plt.tight_layout()
    plt.close(fig)
    return fig


def plot_time_courses(
    datasets: List[Dict[str, Any]],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    value: str | None = None,
    max_nparticipant: int = 10,
    random_seed: int = 12345,
    figsize_width_per_study: int = 4,
    figsize_height_per_specimen: int = 3,
    marker: str = "o",
    markersize: int = DEFAULT_MARKERSIZE,
    line_alpha: float = 0.5,
    show_negative: bool = False,
) -> Figure:
    """
    Plot individual participant shedding trajectories across multiple datasets.

    Creates a grid of faceted plots comparing time courses across different studies.
    Each column represents a different dataset/study, and rows represent different
    specimen types. Individual participant trajectories are shown as lines with markers,
    with different colors for different studies.

    Args:
        datasets: List of raw dataset dictionaries from load_dataset().
        biomarker: Optional filter for a specific biomarker. If None, plots all biomarkers.
        specimen: Optional filter for a specific specimen type. If None, plots all specimens.
        value: Optional filter for value type. Options are "concentration" or "ct".
            If None, plots all data (may raise error if mixed). Use this to avoid
            mixed CT/concentration errors. Defaults to None.
        max_nparticipant: Maximum number of participants per dataset per specimen.
            Defaults to 10 (lower than single dataset due to multiple studies).
        random_seed: Random seed for reproducible sampling. Defaults to 12345.
        figsize_width_per_study: Width of each study subplot in inches. Defaults to 4.
        figsize_height_per_specimen: Height of each specimen row in inches. Defaults to 3.
        marker: Marker style for data points. Defaults to "o".
        markersize: Size of markers. Defaults to DEFAULT_MARKERSIZE (10).
        line_alpha: Transparency of lines (0-1). Defaults to 0.5 for better overlay visibility.
        show_negative: If True, plots "negative" values at y=0. If False, excludes them.
            Defaults to False.

    Returns:
        matplotlib.figure.Figure: The generated figure containing the multi-study comparison.
            Grid layout with rows=specimens, columns=datasets.

    Raises:
        ValueError: If datasets list is empty, or datasets are invalid.
        KeyError: If required analyte information is missing.

    Note:
        - Datasets are distinguished by color using TABLEAU_COLORS
        - All datasets should ideally measure the same analyte for comparison
        - Time axis and measurement units should be consistent across datasets
    """
    # Validate input
    if not datasets or not isinstance(datasets, list):
        raise ValueError("Datasets must be a non-empty list")

    if len(datasets) == 0:
        raise ValueError("Datasets list is empty")

    # Collect data from all datasets
    all_data = []
    for dataset in datasets:
        # Validate each dataset
        if not dataset or not isinstance(dataset, dict):
            raise ValueError("Each dataset must be a non-empty dictionary")

        required_keys = ["analytes", "participants", "dataset_id"]
        missing_keys = [key for key in required_keys if key not in dataset]
        if missing_keys:
            raise ValueError(f"Dataset missing required keys: {missing_keys}")

        # Extract time series data
        time_series_data = []
        for participant_id, participant in enumerate(dataset["participants"], 1):
            measurements = participant.get("measurements", [])
            for measurement in measurements:
                analyte_name = measurement.get("analyte")
                time_series_data.append(
                    {
                        "participant_id": f"{dataset['dataset_id']}_P{participant_id}",
                        "time": measurement.get("time"),
                        "value": measurement.get("value"),
                        "analyte": analyte_name,
                        "dataset_id": dataset["dataset_id"],
                    }
                )

        df_dataset = pd.DataFrame(time_series_data)

        # Add analyte metadata
        for analyte_name, analyte_info in dataset["analytes"].items():
            specimen_value = analyte_info.get("specimen")
            if isinstance(specimen_value, list):
                specimen_value = "+".join(specimen_value)

            # Get limit of detection (could be numeric or "unknown")
            lod = analyte_info.get("limit_of_detection")
            lod_numeric = None
            if lod is not None and lod != "unknown":
                try:
                    lod_numeric = float(lod)
                except (ValueError, TypeError):
                    lod_numeric = None

            df_dataset.loc[df_dataset["analyte"] == analyte_name, "specimen"] = (
                specimen_value
            )
            df_dataset.loc[df_dataset["analyte"] == analyte_name, "unit"] = (
                analyte_info.get("unit")
            )
            df_dataset.loc[df_dataset["analyte"] == analyte_name, "reference_event"] = (
                analyte_info.get("reference_event")
            )
            df_dataset.loc[df_dataset["analyte"] == analyte_name, "biomarker"] = (
                analyte_info.get("biomarker")
            )
            df_dataset.loc[
                df_dataset["analyte"] == analyte_name, "limit_of_detection"
            ] = lod_numeric

        all_data.append(df_dataset)

    # Combine all datasets
    df = pd.concat(all_data, ignore_index=True)

    if df.empty:
        raise ValueError("No measurements found across all datasets")

    # Filter by biomarker if specified
    if biomarker is not None:
        df = df[df["biomarker"] == biomarker]
        if df.empty:
            raise ValueError(f"No measurements found for biomarker '{biomarker}'")

    # Filter by specimen if specified
    if specimen is not None:
        df = df[df["specimen"] == specimen]
        if df.empty:
            raise ValueError(f"No measurements found for specimen '{specimen}'")

    # Convert time to numeric
    df = df[df["time"] != "unknown"].copy()
    df["time_num"] = pd.to_numeric(df["time"], errors="coerce")

    # Determine if each row is CT value or concentration based on unit
    df["is_ct"] = df["unit"].apply(_is_ct_value)

    # Filter by value type if specified
    if value is not None:
        value_lower = value.lower()
        if value_lower == "concentration":
            df = df[~df["is_ct"]].copy()
        elif value_lower == "ct":
            df = df[df["is_ct"]].copy()
        else:
            raise ValueError(
                f"Invalid value '{value}'. " "Must be 'concentration', 'ct', or None."
            )

        if df.empty:
            raise ValueError(f"No {value} data found across datasets after filtering")

    # Handle negative values based on show_negative parameter and data type
    # Track substitution values for caption
    negative_substitution_values = set()

    if show_negative:
        # Substitute negative values: use LOD if available, otherwise 45 for CT or 1 for concentrations
        def substitute_negative(row):
            if row["value"] == NEGATIVE_VALUE:
                # Use limit of detection if available
                if row["limit_of_detection"] is not None and not pd.isna(
                    row["limit_of_detection"]
                ):
                    sub_val = row["limit_of_detection"]
                    negative_substitution_values.add(("LOD", sub_val))
                    return sub_val
                else:
                    # Fallback to default values
                    sub_val = 45 if row["is_ct"] else 1
                    negative_substitution_values.add(("default", sub_val, row["is_ct"]))
                    return sub_val
            return row["value"]

        df["value_num"] = df.apply(substitute_negative, axis=1)
        df["value_num"] = pd.to_numeric(df["value_num"], errors="coerce")
    else:
        df = df[df["value"] != NEGATIVE_VALUE].copy()
        df["value_num"] = pd.to_numeric(df["value"], errors="coerce")

    df = df.dropna(subset=["time_num", "value_num", "specimen"])

    if df.empty:
        raise ValueError("No valid numeric measurements found after filtering")

    # Sample participants per dataset per specimen
    dfs = []
    for (dataset_id, specimen), group in df.groupby(["dataset_id", "specimen"]):
        unique_participants = group["participant_id"].unique()
        if len(unique_participants) > max_nparticipant:
            sampled_participants = (
                pd.Series(unique_participants)
                .sample(n=max_nparticipant, random_state=random_seed)
                .values
            )
            group = group[group["participant_id"].isin(sampled_participants)]
        dfs.append(group)
    df = pd.concat(dfs, ignore_index=True)

    # Create faceted plot: rows=specimens, cols=datasets
    specimens = sorted(df["specimen"].unique())
    dataset_ids = [d["dataset_id"] for d in datasets]
    n_specimens = len(specimens)
    n_datasets = len(dataset_ids)

    fig, axes = plt.subplots(
        n_specimens,
        n_datasets,
        figsize=(
            figsize_width_per_study * n_datasets,
            figsize_height_per_specimen * n_specimens,
        ),
        squeeze=False,
    )

    # Color map for datasets
    color_map = list(mcolors.TABLEAU_COLORS.values())

    # Get reference event and unit for labeling
    reference_event = (
        df["reference_event"].dropna().iloc[0]
        if not df["reference_event"].dropna().empty
        else "reference"
    )
    unit = df["unit"].dropna().iloc[0] if not df["unit"].dropna().empty else ""

    # Check if we have mixed CT and concentration data
    if df["is_ct"].nunique() > 1:
        raise ValueError(
            "Datasets contain mixed CT values and concentrations. "
            "Cannot plot on same axes. Use data_type='concentration' or data_type='ct' to filter."
        )

    # Determine plot type based on data
    is_ct = df["is_ct"].iloc[0] if not df.empty else False

    # Compute global axis ranges for consistency across all panels
    x_min, x_max = df["time_num"].min(), df["time_num"].max()
    y_min, y_max = df["value_num"].min(), df["value_num"].max()

    # Add some padding to ranges
    x_padding = (x_max - x_min) * 0.05 if x_max > x_min else 1
    x_range = (x_min - x_padding, x_max + x_padding)

    if is_ct:
        # For CT values, add padding but keep linear scale
        y_padding = (y_max - y_min) * 0.05 if y_max > y_min else 1
        y_range = (y_min - y_padding, y_max + y_padding)
    else:
        # For concentrations on log scale, add padding in log space
        if y_min > 0:
            log_y_min, log_y_max = np.log10(y_min), np.log10(y_max)
            log_padding = (
                (log_y_max - log_y_min) * 0.1 if log_y_max > log_y_min else 0.5
            )
            y_range = (10 ** (log_y_min - log_padding), 10 ** (log_y_max + log_padding))
        else:
            y_range = (y_min, y_max)

    # Plot each specimen × dataset combination
    for spec_idx, specimen in enumerate(specimens):
        # Get biomarker for this specimen (from first available dataset with this specimen)
        specimen_df = df[df["specimen"] == specimen]
        biomarker = (
            specimen_df["biomarker"].dropna().iloc[0]
            if not specimen_df["biomarker"].dropna().empty
            else ""
        )

        for ds_idx, dataset_id in enumerate(dataset_ids):
            ax = axes[spec_idx, ds_idx]
            subset = df[(df["specimen"] == specimen) & (df["dataset_id"] == dataset_id)]

            if not subset.empty:
                # Plot each participant's trajectory
                color = color_map[ds_idx % len(color_map)]
                for participant_id in subset["participant_id"].unique():
                    participant_data = subset[
                        subset["participant_id"] == participant_id
                    ].sort_values("time_num")
                    ax.plot(
                        participant_data["time_num"],
                        participant_data["value_num"],
                        marker=marker,
                        markersize=markersize,
                        alpha=line_alpha,
                        color=color,
                        linewidth=1.5,
                    )

            # Apply axis scaling and range
            ax.set_xlim(x_range)
            ax.set_ylim(y_range)

            if is_ct:
                # For CT values, reverse y-axis (smaller on top, larger on bottom)
                ax.invert_yaxis()
            else:
                # For concentrations, use log scale
                ax.set_yscale("log")

            # Styling
            if spec_idx == n_specimens - 1:
                ax.set_xlabel(f"Time after {reference_event} (days)")
            if ds_idx == 0:
                # Format specimen name (replace underscores with spaces) and include biomarker
                specimen_display = specimen.replace("_", " ")
                ylabel = (
                    f"{biomarker} - {specimen_display}\n({unit})"
                    if biomarker and unit
                    else (f"{specimen_display}\n({unit})" if unit else specimen_display)
                )
                ax.set_ylabel(ylabel)
            if spec_idx == 0:
                ax.set_title(dataset_id, fontsize=10)
            ax.grid(axis="x", alpha=0.3)

    fig.suptitle("Time Course Comparison Across Studies", fontsize=14, y=0.995)

    # Add caption for negative value substitution if applicable
    if show_negative and negative_substitution_values:
        caption_parts = []
        # Check for LOD substitutions
        lod_values = [
            item[1] for item in negative_substitution_values if item[0] == "LOD"
        ]
        if lod_values:
            lod_str = ", ".join([f"{v:.2g}" for v in sorted(set(lod_values))])
            caption_parts.append(f"limit of detection ({lod_str})")

        # Check for default substitutions
        default_items = [
            item for item in negative_substitution_values if item[0] == "default"
        ]
        if default_items:
            # Group by CT vs concentration
            ct_vals = [item[1] for item in default_items if item[2]]
            conc_vals = [item[1] for item in default_items if not item[2]]
            if ct_vals:
                caption_parts.append(f"CT value of {ct_vals[0]}")
            if conc_vals:
                caption_parts.append(f"concentration of {conc_vals[0]}")

        if caption_parts:
            caption = "Note: Negative values were substituted with " + " or ".join(
                caption_parts
            )
            fig.text(
                0.5, -0.05, caption, ha="center", fontsize=9, style="italic", wrap=True
            )

    plt.tight_layout()
    plt.close(fig)
    return fig


def plot_shedding_heatmap(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    value: str | None = None,
    time_bin_size: float = 1.0,
    time_range: tuple[float, float] | None = None,
    sort_by: str = "first_positive",
    max_nparticipant: int = 50,
    random_seed: int = 12345,
    figsize: tuple[int, int] | None = None,
    cmap: str | None = None,
    show_negative: bool = True,
    show_colorbar: bool = True,
    show_participant_labels: bool = False,
) -> Figure:
    """
    Plot a heatmap of shedding intensity over time across participants.

    Creates a heatmap where rows represent participants, columns represent time bins,
    and color intensity represents measurement values. Provides a compact overview
    of shedding patterns across all participants.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, plots all biomarkers.
        specimen: Optional filter for a specific specimen type. If None, uses first specimen found.
        value: Optional filter for value type. Options are "concentration" or "ct".
            If None, plots all data (may raise error if mixed). Defaults to None.
        time_bin_size: Size of time bins in days. Defaults to 1.0.
        time_range: Optional tuple (min_time, max_time) to limit the time axis.
            If None, uses the full range of data.
        sort_by: How to sort participants. Options:
            - "first_positive": Sort by time of first positive measurement (default)
            - "peak_time": Sort by time of peak value
            - "peak_value": Sort by peak measurement value
            - "participant_id": Sort by participant ID (original order)
        max_nparticipant: Maximum number of participants to display. If exceeded,
            randomly samples participants. Defaults to 50.
        random_seed: Random seed for reproducible sampling. Defaults to 12345.
        figsize: Figure size as (width, height). If None, automatically calculated.
        cmap: Colormap name. If None, uses "YlOrRd" for concentrations, "YlOrRd_r" for CT.
        show_negative: If True, shows negative values using limit of detection or default.
            If False, negative values appear as missing (white). Defaults to True.
        show_colorbar: If True, displays a colorbar. Defaults to True.
        show_participant_labels: If True, shows participant IDs on y-axis. Defaults to False.

    Returns:
        matplotlib.figure.Figure: The generated heatmap figure.

    Raises:
        ValueError: If dataset is missing required keys, is empty, or has no valid data.
    """
    # Validate input
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "dataset_id"]
    missing_keys = [key for key in required_keys if key not in dataset]
    if missing_keys:
        raise ValueError(f"Dataset missing required keys: {missing_keys}")

    if not dataset["participants"]:
        raise ValueError("Dataset has no participants")

    # Extract time series data from raw dataset
    time_series_data = []
    for participant_id, participant in enumerate(dataset["participants"], 1):
        measurements = participant.get("measurements", [])
        for measurement in measurements:
            analyte_name = measurement.get("analyte")
            time_series_data.append(
                {
                    "participant_id": participant_id,
                    "time": measurement.get("time"),
                    "value": measurement.get("value"),
                    "analyte": analyte_name,
                }
            )

    if not time_series_data:
        raise ValueError("Dataset has no measurements")

    df = pd.DataFrame(time_series_data)

    # Convert time to numeric, filtering out "unknown" values
    df = df[df["time"] != "unknown"].copy()
    df["time_num"] = pd.to_numeric(df["time"], errors="coerce")

    # Join with analyte metadata to get specimen and unit information
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_val = analyte_info.get("specimen")
        if isinstance(specimen_val, list):
            specimen_val = "+".join(specimen_val)

        lod = analyte_info.get("limit_of_detection")
        lod_numeric = None
        if lod is not None and lod != "unknown":
            try:
                lod_numeric = float(lod)
            except (ValueError, TypeError):
                lod_numeric = None

        analyte_metadata[analyte_name] = {
            "specimen": specimen_val,
            "unit": analyte_info.get("unit"),
            "reference_event": analyte_info.get("reference_event"),
            "biomarker": analyte_info.get("biomarker"),
            "limit_of_detection": lod_numeric,
        }

    # Add metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
    )
    df["unit"] = df["analyte"].map(lambda x: analyte_metadata.get(x, {}).get("unit"))
    df["reference_event"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("reference_event")
    )
    df["biomarker"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("biomarker")
    )
    df["limit_of_detection"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("limit_of_detection")
    )

    # Filter by biomarker if specified
    if biomarker is not None:
        df = df[df["biomarker"] == biomarker]
        if df.empty:
            raise ValueError(f"No measurements found for biomarker '{biomarker}'")

    # Filter by specimen if specified
    if specimen is not None:
        df = df[df["specimen"] == specimen]
        if df.empty:
            raise ValueError(f"No measurements found for specimen '{specimen}'")

    # Determine if each row is CT value or concentration based on unit
    df["is_ct"] = df["unit"].apply(_is_ct_value)

    # Filter by value type if specified
    if value is not None:
        value_lower = value.lower()
        if value_lower == "concentration":
            df = df[~df["is_ct"]].copy()
        elif value_lower == "ct":
            df = df[df["is_ct"]].copy()
        else:
            raise ValueError(
                f"Invalid value '{value}'. " "Must be 'concentration', 'ct', or None."
            )
        if df.empty:
            raise ValueError(f"No {value} data found in dataset after filtering")

    # Check for mixed CT and concentration data
    if df["is_ct"].nunique() > 1:
        raise ValueError(
            "Dataset contains mixed CT values and concentrations. "
            "Use value='concentration' or value='ct' to filter."
        )

    is_ct = df["is_ct"].iloc[0] if not df.empty else False

    # Handle negative values - track them separately for distinct coloring
    df["is_negative"] = df["value"] == NEGATIVE_VALUE

    if show_negative:
        # For negative values, we'll show them with a distinct color (skyblue)
        # Set numeric value to NaN for negatives, they'll be overlaid later
        df["value_num"] = df["value"].apply(
            lambda x: np.nan if x == NEGATIVE_VALUE else x
        )
        df["value_num"] = pd.to_numeric(df["value_num"], errors="coerce")
    else:
        df["value_num"] = df["value"].apply(
            lambda x: np.nan if x == NEGATIVE_VALUE else x
        )
        df["value_num"] = pd.to_numeric(df["value_num"], errors="coerce")

    # Drop rows with NaN time
    df = df.dropna(subset=["time_num"])

    if df.empty:
        raise ValueError("No valid measurements found after filtering")

    # Apply time range filter if specified
    if time_range is not None:
        df = df[(df["time_num"] >= time_range[0]) & (df["time_num"] <= time_range[1])]
        if df.empty:
            raise ValueError(f"No measurements found in time range {time_range}")

    # Sample participants if needed
    unique_participants = df["participant_id"].unique()
    if len(unique_participants) > max_nparticipant:
        sampled_participants = (
            pd.Series(unique_participants)
            .sample(n=max_nparticipant, random_state=random_seed)
            .values
        )
        df = df[df["participant_id"].isin(sampled_participants)]

    # Create time bins
    time_min = df["time_num"].min()
    time_max = df["time_num"].max()
    if time_range is not None:
        time_min, time_max = time_range

    # Floor time_min and ceil time_max to bin_size boundaries for clean tick labels
    time_min = np.floor(time_min / time_bin_size) * time_bin_size
    time_max = np.ceil(time_max / time_bin_size) * time_bin_size

    # Create bin edges
    bins = np.arange(time_min - time_bin_size, time_max + time_bin_size, time_bin_size)
    bin_labels = bins[:-1] + time_bin_size / 2  # Center of each bin

    df["time_bin"] = pd.cut(
        df["time_num"], bins=bins, labels=bin_labels, include_lowest=True
    )

    # Aggregate values within each bin (take mean if multiple measurements)
    pivot_df = (
        df.groupby(["participant_id", "time_bin"], observed=False)["value_num"]
        .mean()
        .unstack()
    )

    # Create a pivot table for negative values (True if any measurement in bin was negative)
    negative_pivot = (
        df.groupby(["participant_id", "time_bin"], observed=False)["is_negative"]
        .any()
        .unstack()
    )

    # Reindex to ensure all bins are present (fills missing bins with NaN/False)
    pivot_df = pivot_df.reindex(columns=bin_labels, fill_value=np.nan)
    negative_pivot = negative_pivot.reindex(columns=bin_labels, fill_value=False)

    # Sort participants based on sort_by parameter
    if sort_by == "first_positive":
        # Sort by time of first non-NaN value
        first_positive_time = pivot_df.apply(
            lambda row: row.first_valid_index(), axis=1
        )
        sort_order = first_positive_time.sort_values().index
    elif sort_by == "peak_time":
        # Sort by time of peak value
        if is_ct:
            # For CT, peak is minimum value
            peak_time = pivot_df.apply(
                lambda row: row.idxmin() if row.notna().any() else np.nan, axis=1
            )
        else:
            # For concentration, peak is maximum value
            peak_time = pivot_df.apply(
                lambda row: row.idxmax() if row.notna().any() else np.nan, axis=1
            )
        sort_order = peak_time.sort_values().index
    elif sort_by == "peak_value":
        # Sort by peak measurement value
        if is_ct:
            # For CT, lower is higher viral load, so sort ascending
            peak_value = pivot_df.min(axis=1)
            sort_order = peak_value.sort_values(ascending=True).index
        else:
            # For concentration, sort descending
            peak_value = pivot_df.max(axis=1)
            sort_order = peak_value.sort_values(ascending=False).index
    elif sort_by == "participant_id":
        sort_order = pivot_df.index.sort_values()
    else:
        raise ValueError(
            f"Invalid sort_by '{sort_by}'. "
            "Must be 'first_positive', 'peak_time', 'peak_value', or 'participant_id'."
        )

    pivot_df = pivot_df.loc[sort_order]
    negative_pivot = negative_pivot.loc[sort_order]

    # Determine figure size
    if figsize is None:
        n_participants = len(pivot_df)
        n_time_bins = len(pivot_df.columns)
        width = max(8, n_time_bins * 0.3 + 2)
        height = max(6, n_participants * 0.2 + 2)
        figsize = (width, height)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Determine colormap
    # For CT: use reversed colormap so smaller values = darker color (higher viral load)
    # For concentration: use log10 scale, higher values = darker color
    if cmap is None:
        cmap = "YlOrRd_r" if is_ct else "YlOrRd"

    # Prepare heatmap data with appropriate scaling
    heatmap_data = pivot_df.values.copy()

    if is_ct:
        # For CT values: keep original values, reversed colormap handles color direction
        pass
    else:
        # For concentration: apply log10 scale
        # Replace zeros and negatives with NaN to avoid log errors
        heatmap_data = np.where(heatmap_data > 0, np.log10(heatmap_data), np.nan)

    # Create heatmap
    im = ax.imshow(
        heatmap_data,
        aspect="auto",
        cmap=cmap,
        interpolation="nearest",
    )

    # Overlay negative values with skyblue color
    if show_negative:
        negative_mask = negative_pivot.values
        for i in range(negative_mask.shape[0]):
            for j in range(negative_mask.shape[1]):
                if negative_mask[i, j] == True:
                    # Draw a rectangle for negative values
                    rect = plt.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1, facecolor="skyblue", edgecolor="none"
                    )
                    ax.add_patch(rect)

    # Set axis labels with larger font
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.set_ylabel("Participants", fontsize=12)

    # Set x-axis ticks
    # Derive tick labels from pivot_df.columns to ensure consistency
    # Columns are bin centers, so subtract half bin_size to get bin start times
    n_cols = len(pivot_df.columns)
    n_xticks = min(5, n_cols)
    # Select evenly spaced tick positions
    xtick_indices = np.round(np.linspace(0, n_cols - 1, n_xticks)).astype(int)
    # Convert bin centers to bin start times for more intuitive labels
    xtick_labels = [
        f"{round(float(pivot_df.columns[i]) + time_bin_size / 2)}"
        for i in xtick_indices
    ]
    ax.set_xticks(xtick_indices)
    ax.set_xticklabels(xtick_labels, fontsize=10)

    # Set y-axis ticks
    if show_participant_labels:
        ax.set_yticks(range(len(pivot_df)))
        ax.set_yticklabels([f"P{pid}" for pid in pivot_df.index], fontsize=10)
    else:
        n_yticks = min(10, len(pivot_df))
        ytick_positions = np.linspace(0, len(pivot_df) - 1, n_yticks, dtype=int)
        ax.set_yticks(ytick_positions)
        ax.set_yticklabels([f"{i+1}" for i in ytick_positions], fontsize=10)

    # Add colorbar
    if show_colorbar:
        unit = df["unit"].dropna().iloc[0] if not df["unit"].dropna().empty else ""
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        if is_ct:
            cbar_label = f"CT value ({unit})\n(darker = higher viral load)"
        else:
            cbar_label = f"log₁₀ Concentration ({unit})"
        cbar.set_label(cbar_label, fontsize=11)
        cbar.ax.tick_params(labelsize=10)

    # Add legend for negative values if shown
    if show_negative and negative_pivot.values.any():
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor="skyblue", edgecolor="none", label="Negative")
        ]
        ax.legend(handles=legend_elements, loc="upper right", fontsize=10)

    # Add title
    dataset_id = dataset.get("dataset_id", "Dataset")
    specimen_display = (
        df["specimen"].dropna().iloc[0] if not df["specimen"].dropna().empty else ""
    )
    biomarker_display = (
        df["biomarker"].dropna().iloc[0] if not df["biomarker"].dropna().empty else ""
    )

    title_parts = [f"Shedding Heatmap: {dataset_id}"]
    if biomarker_display:
        title_parts.append(f"({biomarker_display}")
        if specimen_display:
            title_parts[-1] += f" - {specimen_display})"
        else:
            title_parts[-1] += ")"
    elif specimen_display:
        title_parts.append(f"({specimen_display})")

    ax.set_title(" ".join(title_parts), fontsize=14)

    plt.tight_layout()
    plt.close(fig)
    return fig


def plot_mean_trajectory(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    value: str | None = None,
    central_tendency: str = "mean",
    uncertainty: str = "95ci",
    time_bin_size: float = 1.0,
    time_range: tuple[float, float] | None = None,
    min_observations: int = 3,
    figsize: tuple[int, int] = (10, 6),
    line_color: str = "steelblue",
    fill_alpha: float = 0.3,
    show_individual: bool = False,
    individual_alpha: float = 0.5,
    show_n: bool = True,
) -> Figure:
    """
    Plot mean/median trajectory with confidence bands across participants.

    Creates a line plot showing the central tendency (mean or median) of measurements
    over time, with a shaded band showing the uncertainty range. Useful for visualizing
    the "typical" shedding pattern across all participants.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, uses first biomarker found.
        specimen: Optional filter for a specific specimen type. If None, uses first specimen found.
        value: Optional filter for value type. Options are "concentration" or "ct".
            If None, plots all data (may raise error if mixed). Defaults to None.
        central_tendency: Measure of central tendency. Options:
            - "mean": Arithmetic mean (default)
            - "median": Median value
        uncertainty: Type of uncertainty band. Options:
            - "95ci": 95% confidence interval (default)
            - "iqr": Interquartile range (25th-75th percentile)
            - "sd": Standard deviation (mean ± 1 SD)
            - "range": Full range (min-max)
        time_bin_size: Size of time bins in days for aggregating measurements. Defaults to 1.0.
        time_range: Optional tuple (min_time, max_time) to limit the time axis.
            If None, uses the full range of data.
        min_observations: Minimum number of observations required per time bin to include
            in the plot. Bins with fewer observations are excluded. Defaults to 3.
        figsize: Figure size as (width, height). Defaults to (10, 6).
        line_color: Color for the central tendency line and fill. Defaults to "steelblue".
        fill_alpha: Transparency of the uncertainty band (0-1). Defaults to 0.3.
        show_individual: If True, shows individual participant trajectories in background.
            Defaults to False.
        individual_alpha: Transparency of individual trajectories (0-1). Defaults to 0.5.
        show_n: If True, shows the number of observations at each time point. Defaults to True.

    Returns:
        matplotlib.figure.Figure: The generated figure containing the mean trajectory plot.

    Raises:
        ValueError: If dataset is missing required keys, is empty, or has no valid data.
    """
    # Validate input
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "dataset_id"]
    missing_keys = [key for key in required_keys if key not in dataset]
    if missing_keys:
        raise ValueError(f"Dataset missing required keys: {missing_keys}")

    if not dataset["participants"]:
        raise ValueError("Dataset has no participants")

    # Validate parameters
    if central_tendency not in ["mean", "median"]:
        raise ValueError(
            f"Invalid central_tendency '{central_tendency}'. "
            "Must be 'mean' or 'median'."
        )

    if uncertainty not in ["95ci", "iqr", "sd", "range"]:
        raise ValueError(
            f"Invalid uncertainty '{uncertainty}'. "
            "Must be '95ci', 'iqr', 'sd', or 'range'."
        )

    # Extract time series data from raw dataset
    time_series_data = []
    for participant_id, participant in enumerate(dataset["participants"], 1):
        measurements = participant.get("measurements", [])
        for measurement in measurements:
            analyte_name = measurement.get("analyte")
            time_series_data.append(
                {
                    "participant_id": participant_id,
                    "time": measurement.get("time"),
                    "value": measurement.get("value"),
                    "analyte": analyte_name,
                }
            )

    if not time_series_data:
        raise ValueError("Dataset has no measurements")

    df = pd.DataFrame(time_series_data)

    # Convert time to numeric, filtering out "unknown" values
    df = df[df["time"] != "unknown"].copy()
    df["time_num"] = pd.to_numeric(df["time"], errors="coerce")

    # Join with analyte metadata to get specimen and unit information
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        lod = analyte_info.get("limit_of_detection")
        lod_numeric = None
        if lod is not None and lod != "unknown":
            try:
                lod_numeric = float(lod)
            except (ValueError, TypeError):
                lod_numeric = None

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "unit": analyte_info.get("unit"),
            "reference_event": analyte_info.get("reference_event"),
            "biomarker": analyte_info.get("biomarker"),
            "limit_of_detection": lod_numeric,
        }

    # Add metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
    )
    df["unit"] = df["analyte"].map(lambda x: analyte_metadata.get(x, {}).get("unit"))
    df["reference_event"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("reference_event")
    )
    df["biomarker"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("biomarker")
    )
    df["limit_of_detection"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("limit_of_detection")
    )

    # Filter by biomarker if specified
    if biomarker is not None:
        df = df[df["biomarker"] == biomarker]
        if df.empty:
            raise ValueError(f"No measurements found for biomarker '{biomarker}'")

    # Filter by specimen if specified
    if specimen is not None:
        df = df[df["specimen"] == specimen]
        if df.empty:
            raise ValueError(f"No measurements found for specimen '{specimen}'")

    # Determine if each row is CT value or concentration based on unit
    df["is_ct"] = df["unit"].apply(_is_ct_value)

    # Filter by value type if specified
    if value is not None:
        value_lower = value.lower()
        if value_lower == "concentration":
            df = df[~df["is_ct"]].copy()
        elif value_lower == "ct":
            df = df[df["is_ct"]].copy()
        else:
            raise ValueError(
                f"Invalid value '{value}'. " "Must be 'concentration', 'ct', or None."
            )
        if df.empty:
            raise ValueError(f"No {value} data found in dataset after filtering")

    # Check for mixed CT and concentration data
    if df["is_ct"].nunique() > 1:
        raise ValueError(
            "Dataset contains mixed CT values and concentrations. "
            "Use value='concentration' or value='ct' to filter."
        )

    is_ct = df["is_ct"].iloc[0] if not df.empty else False

    # Exclude negative values for trajectory calculation
    df = df[df["value"] != NEGATIVE_VALUE].copy()
    df["value_num"] = pd.to_numeric(df["value"], errors="coerce")

    # Drop rows with NaN time or value
    df = df.dropna(subset=["time_num", "value_num"])

    if df.empty:
        raise ValueError("No valid numeric measurements found after filtering")

    # Apply time range filter if specified
    if time_range is not None:
        df = df[(df["time_num"] >= time_range[0]) & (df["time_num"] <= time_range[1])]
        if df.empty:
            raise ValueError(f"No measurements found in time range {time_range}")

    # Create time bins
    time_min = df["time_num"].min()
    time_max = df["time_num"].max()
    if time_range is not None:
        time_min, time_max = time_range

    # Floor time_min and ceil time_max to bin_size boundaries
    time_min = np.floor(time_min / time_bin_size) * time_bin_size
    time_max = np.ceil(time_max / time_bin_size) * time_bin_size

    bins = np.arange(time_min, time_max + time_bin_size, time_bin_size)
    bin_centers = bins[:-1] + time_bin_size / 2

    df["time_bin"] = pd.cut(
        df["time_num"], bins=bins, labels=bin_centers, include_lowest=True
    )

    # Calculate statistics per time bin
    def calc_stats(group):
        values = group["value_num"].values
        n = len(values)

        if n < min_observations:
            return pd.Series(
                {
                    "n": n,
                    "center": np.nan,
                    "lower": np.nan,
                    "upper": np.nan,
                }
            )

        if central_tendency == "mean":
            center = np.mean(values)
        else:  # median
            center = np.median(values)

        if uncertainty == "95ci":
            # 95% confidence interval
            sem = np.std(values, ddof=1) / np.sqrt(n)
            lower = center - 1.96 * sem
            upper = center + 1.96 * sem
        elif uncertainty == "iqr":
            # Interquartile range
            lower = np.percentile(values, 25)
            upper = np.percentile(values, 75)
        elif uncertainty == "sd":
            # Standard deviation
            sd = np.std(values, ddof=1)
            lower = center - sd
            upper = center + sd
        elif uncertainty == "range":
            # Full range
            lower = np.min(values)
            upper = np.max(values)

        return pd.Series(
            {
                "n": n,
                "center": center,
                "lower": lower,
                "upper": upper,
            }
        )

    stats_df = df.groupby("time_bin", observed=False).apply(calc_stats).reset_index()
    stats_df["time"] = stats_df["time_bin"].astype(float)

    # Remove bins with insufficient observations
    stats_df = stats_df.dropna(subset=["center"])

    if stats_df.empty:
        raise ValueError(
            f"No time bins have at least {min_observations} observations. "
            "Try reducing min_observations or using a larger time_bin_size."
        )

    # Sort by time
    stats_df = stats_df.sort_values("time")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot individual trajectories in background if requested
    if show_individual:
        for participant_id in df["participant_id"].unique():
            participant_data = df[df["participant_id"] == participant_id].sort_values(
                "time_num"
            )
            ax.plot(
                participant_data["time_num"],
                participant_data["value_num"],
                color="gray",
                alpha=individual_alpha,
                linewidth=0.5,
            )

    # Plot uncertainty band
    ax.fill_between(
        stats_df["time"],
        stats_df["lower"],
        stats_df["upper"],
        color=line_color,
        alpha=fill_alpha,
        label=f"{uncertainty.upper()}" if uncertainty != "95ci" else "95% CI",
    )

    # Plot central tendency line
    ax.plot(
        stats_df["time"],
        stats_df["center"],
        color=line_color,
        linewidth=2,
        label=central_tendency.capitalize(),
    )

    # Set axis scaling
    if is_ct:
        ax.invert_yaxis()
    else:
        ax.set_yscale("log")

    # Labels and title
    reference_event = (
        df["reference_event"].dropna().iloc[0]
        if not df["reference_event"].dropna().empty
        else "reference"
    )
    unit = df["unit"].dropna().iloc[0] if not df["unit"].dropna().empty else ""
    specimen_display = (
        df["specimen"].dropna().iloc[0] if not df["specimen"].dropna().empty else ""
    )
    biomarker_display = (
        df["biomarker"].dropna().iloc[0] if not df["biomarker"].dropna().empty else ""
    )

    ax.set_xlabel(f"Time after {reference_event} (days)", fontsize=12)
    if is_ct:
        ax.set_ylabel(f"CT value ({unit})", fontsize=12)
    else:
        ax.set_ylabel(f"Concentration ({unit})", fontsize=12)

    # Title
    dataset_id = dataset.get("dataset_id", "Dataset")
    title_parts = [f"{central_tendency.capitalize()} Trajectory: {dataset_id}"]
    if biomarker_display:
        title_parts.append(f"({biomarker_display}")
        if specimen_display:
            title_parts[-1] += f" - {specimen_display})"
        else:
            title_parts[-1] += ")"
    elif specimen_display:
        title_parts.append(f"({specimen_display})")

    ax.set_title(" ".join(title_parts), fontsize=14)

    # Add sample size annotations if requested
    if show_n:
        # Add n values at the top of the plot (vertical text to avoid overlap)
        # Add offset to avoid clipping at axis boundary
        y_lim = ax.get_ylim()
        if is_ct:
            # For inverted CT axis: y_lim[1] is smaller value (top of plot)
            # Add offset towards larger values (move down visually)
            y_range = abs(y_lim[1] - y_lim[0])
            y_pos = y_lim[1] + y_range * 0.08
        else:
            # For log scale concentrations, apply offset in log space
            log_range = np.log10(y_lim[1]) - np.log10(y_lim[0])
            y_pos = 10 ** (np.log10(y_lim[1]) - log_range * 0.08)
        for _, row in stats_df.iterrows():
            ax.annotate(
                f"n={int(row['n'])}",
                xy=(row["time"], y_pos),
                fontsize=8,
                ha="center",
                va="bottom",
                alpha=0.7,
                rotation=90,
            )

    ax.legend(loc="best", fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.close(fig)
    return fig


def plot_value_distribution_by_time(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    value: str | None = None,
    plot_type: str = "box",
    time_bin_size: float = 1.0,
    time_range: tuple[float, float] | None = None,
    min_observations: int = 1,
    figsize: tuple[int, int] = (12, 6),
    color: str = "steelblue",
    show_points: bool = False,
    point_alpha: float = 0.5,
    show_n: bool = True,
    show_mean: bool = False,
    widths: float | None = None,
) -> Figure:
    """
    Plot distribution of measurement values at each time bin.

    Creates box plots or violin plots showing how measurement values are distributed
    at different time points. Useful for understanding variability in shedding
    patterns across participants at each stage of infection.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, uses first biomarker found.
        specimen: Optional filter for a specific specimen type. If None, uses first specimen found.
        value: Optional filter for value type. Options are "concentration" or "ct".
            If None, plots all data (may raise error if mixed). Defaults to None.
        plot_type: Type of distribution plot. Options:
            - "box": Box plot showing median, quartiles, and outliers (default)
            - "violin": Violin plot showing full distribution shape
        time_bin_size: Size of time bins in days for grouping measurements. Defaults to 1.0.
        time_range: Optional tuple (min_time, max_time) to limit the time axis.
            If None, uses the full range of data.
        min_observations: Minimum number of observations required per time bin to include
            in the plot. Bins with fewer observations are excluded. Defaults to 1.
        figsize: Figure size as (width, height). Defaults to (12, 6).
        color: Color for the box/violin plots. Defaults to "steelblue".
        show_points: If True, overlays individual data points on the distribution plots.
            Defaults to False.
        point_alpha: Transparency of individual points (0-1). Only used if show_points=True.
            Defaults to 0.5.
        show_n: If True, shows the number of observations at each time bin. Defaults to True.
        show_mean: If True, shows mean markers on box plots. Defaults to False.
        widths: Width of box/violin plots. If None, automatically calculated based on
            time_bin_size. Defaults to None.

    Returns:
        matplotlib.figure.Figure: The generated figure containing the distribution plots.

    Raises:
        ValueError: If dataset is missing required keys, is empty, has no valid data,
            or if invalid parameters are provided.

    Note:
        - Negative values are excluded from the distribution calculations
        - For CT values, the y-axis is inverted (lower CT = higher viral load)
        - For concentration values, a log scale is used on the y-axis
    """
    # Validate input
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "dataset_id"]
    missing_keys = [key for key in required_keys if key not in dataset]
    if missing_keys:
        raise ValueError(f"Dataset missing required keys: {missing_keys}")

    if not dataset["participants"]:
        raise ValueError("Dataset has no participants")

    # Validate plot_type
    if plot_type not in ["box", "violin"]:
        raise ValueError(
            f"Invalid plot_type '{plot_type}'. " "Must be 'box' or 'violin'."
        )

    # Extract time series data from raw dataset
    time_series_data = []
    for participant_id, participant in enumerate(dataset["participants"], 1):
        measurements = participant.get("measurements", [])
        for measurement in measurements:
            analyte_name = measurement.get("analyte")
            time_series_data.append(
                {
                    "participant_id": participant_id,
                    "time": measurement.get("time"),
                    "value": measurement.get("value"),
                    "analyte": analyte_name,
                }
            )

    if not time_series_data:
        raise ValueError("Dataset has no measurements")

    df = pd.DataFrame(time_series_data)

    # Convert time to numeric, filtering out "unknown" values
    df = df[df["time"] != "unknown"].copy()
    df["time_num"] = pd.to_numeric(df["time"], errors="coerce")

    # Join with analyte metadata to get specimen and unit information
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        lod = analyte_info.get("limit_of_detection")
        lod_numeric = None
        if lod is not None and lod != "unknown":
            try:
                lod_numeric = float(lod)
            except (ValueError, TypeError):
                lod_numeric = None

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "unit": analyte_info.get("unit"),
            "reference_event": analyte_info.get("reference_event"),
            "biomarker": analyte_info.get("biomarker"),
            "limit_of_detection": lod_numeric,
        }

    # Add metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
    )
    df["unit"] = df["analyte"].map(lambda x: analyte_metadata.get(x, {}).get("unit"))
    df["reference_event"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("reference_event")
    )
    df["biomarker"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("biomarker")
    )
    df["limit_of_detection"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("limit_of_detection")
    )

    # Filter by biomarker if specified
    if biomarker is not None:
        df = df[df["biomarker"] == biomarker]
        if df.empty:
            raise ValueError(f"No measurements found for biomarker '{biomarker}'")

    # Filter by specimen if specified
    if specimen is not None:
        df = df[df["specimen"] == specimen]
        if df.empty:
            raise ValueError(f"No measurements found for specimen '{specimen}'")

    # Determine if each row is CT value or concentration based on unit
    df["is_ct"] = df["unit"].apply(_is_ct_value)

    # Filter by value type if specified
    if value is not None:
        value_lower = value.lower()
        if value_lower == "concentration":
            df = df[~df["is_ct"]].copy()
        elif value_lower == "ct":
            df = df[df["is_ct"]].copy()
        else:
            raise ValueError(
                f"Invalid value '{value}'. " "Must be 'concentration', 'ct', or None."
            )
        if df.empty:
            raise ValueError(f"No {value} data found in dataset after filtering")

    # Check for mixed CT and concentration data
    if df["is_ct"].nunique() > 1:
        raise ValueError(
            "Dataset contains mixed CT values and concentrations. "
            "Use value='concentration' or value='ct' to filter."
        )

    is_ct = df["is_ct"].iloc[0] if not df.empty else False

    # Exclude negative values for distribution calculation
    df = df[df["value"] != NEGATIVE_VALUE].copy()
    df["value_num"] = pd.to_numeric(df["value"], errors="coerce")

    # Drop rows with NaN time or value
    df = df.dropna(subset=["time_num", "value_num"])

    if df.empty:
        raise ValueError("No valid numeric measurements found after filtering")

    # Apply time range filter if specified
    if time_range is not None:
        df = df[(df["time_num"] >= time_range[0]) & (df["time_num"] <= time_range[1])]
        if df.empty:
            raise ValueError(f"No measurements found in time range {time_range}")

    # Create time bins centered at integers (or multiples of bin_size)
    time_min = df["time_num"].min()
    time_max = df["time_num"].max()
    if time_range is not None:
        time_min, time_max = time_range

    # Round time_min down and time_max up to nearest multiple of bin_size to get bin centers
    center_min = np.floor(time_min / time_bin_size) * time_bin_size
    center_max = np.ceil(time_max / time_bin_size) * time_bin_size

    # Create bin edges at ±bin_size/2 around centers (e.g., for center=1 with bin_size=1: edges are 0.5 and 1.5)
    bins = np.arange(
        center_min - time_bin_size / 2, center_max + time_bin_size, time_bin_size
    )
    bin_centers = np.arange(center_min, center_max + time_bin_size / 2, time_bin_size)

    df["time_bin"] = pd.cut(
        df["time_num"], bins=bins, labels=bin_centers.astype(int), include_lowest=True
    )
    df["time_bin_num"] = df["time_bin"].astype(float)

    # Group by time bin and filter by min_observations
    bin_counts = df.groupby("time_bin_num").size()
    valid_bins = bin_counts[bin_counts >= min_observations].index
    df = df[df["time_bin_num"].isin(valid_bins)]

    if df.empty:
        raise ValueError(
            f"No time bins have at least {min_observations} observations. "
            "Try reducing min_observations or using a larger time_bin_size."
        )

    # Get unique time bins and sort
    unique_bins = sorted(df["time_bin_num"].unique())

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data for plotting
    plot_data = [
        df[df["time_bin_num"] == bin_val]["value_num"].values for bin_val in unique_bins
    ]

    # Calculate widths
    if widths is None:
        widths = time_bin_size * 0.7

    # Create distribution plots
    if plot_type == "box":
        bp = ax.boxplot(
            plot_data,
            positions=unique_bins,
            widths=widths,
            patch_artist=True,
            showmeans=show_mean,
            meanprops=dict(
                marker="D",
                markerfacecolor="white",
                markeredgecolor="black",
                markersize=6,
            ),
        )
        # Style box plots
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(2)
        for whisker in bp["whiskers"]:
            whisker.set_color("gray")
        for cap in bp["caps"]:
            cap.set_color("gray")
        for flier in bp["fliers"]:
            flier.set_markerfacecolor(color)
            flier.set_markeredgecolor(color)
            flier.set_alpha(0.5)

    elif plot_type == "violin":
        vp = ax.violinplot(
            plot_data,
            positions=unique_bins,
            widths=widths,
            showmeans=show_mean,
            showmedians=True,
            showextrema=True,
        )
        # Style violin plots
        for body in vp["bodies"]:
            body.set_facecolor(color)
            body.set_alpha(0.7)
        if "cmedians" in vp:
            vp["cmedians"].set_color("white")
            vp["cmedians"].set_linewidth(2)
        if "cmeans" in vp:
            vp["cmeans"].set_color("black")
        if "cbars" in vp:
            vp["cbars"].set_color("gray")
        if "cmins" in vp:
            vp["cmins"].set_color("gray")
        if "cmaxes" in vp:
            vp["cmaxes"].set_color("gray")

    # Overlay individual points if requested
    if show_points:
        for bin_val, data in zip(unique_bins, plot_data):
            # Add jitter to x positions for visibility
            jitter = np.random.RandomState(12345).uniform(
                -widths / 3, widths / 3, size=len(data)
            )
            ax.scatter(
                [bin_val] * len(data) + jitter,
                data,
                color="black",
                alpha=point_alpha,
                s=15,
                zorder=3,
            )

    # Set axis scaling
    if is_ct:
        ax.invert_yaxis()
    else:
        ax.set_yscale("log")

    # Set x-axis ticks to show actual day values with reasonable spacing
    n_bins = len(unique_bins)
    if n_bins <= 10:
        # Show all ticks if there aren't too many
        ax.set_xticks(unique_bins)
        ax.set_xticklabels([int(b) for b in unique_bins])
    else:
        # Select a subset of ticks for readability
        n_ticks = min(10, n_bins)
        tick_indices = np.round(np.linspace(0, n_bins - 1, n_ticks)).astype(int)
        tick_positions = [unique_bins[i] for i in tick_indices]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([int(p) for p in tick_positions])

    # Labels and title
    reference_event = (
        df["reference_event"].dropna().iloc[0]
        if not df["reference_event"].dropna().empty
        else "reference"
    )
    unit = df["unit"].dropna().iloc[0] if not df["unit"].dropna().empty else ""
    specimen_display = (
        df["specimen"].dropna().iloc[0] if not df["specimen"].dropna().empty else ""
    )
    biomarker_display = (
        df["biomarker"].dropna().iloc[0] if not df["biomarker"].dropna().empty else ""
    )

    ax.set_xlabel(f"Time after {reference_event} (days)", fontsize=12)
    if is_ct:
        ax.set_ylabel(f"CT value ({unit})", fontsize=12)
    else:
        ax.set_ylabel(f"Concentration ({unit})", fontsize=12)

    # Title
    dataset_id = dataset.get("dataset_id", "Dataset")
    plot_type_display = "Box Plot" if plot_type == "box" else "Violin Plot"
    title_parts = [f"Value Distribution ({plot_type_display}): {dataset_id}"]
    if biomarker_display:
        title_parts.append(f"({biomarker_display}")
        if specimen_display:
            title_parts[-1] += f" - {specimen_display})"
        else:
            title_parts[-1] += ")"
    elif specimen_display:
        title_parts.append(f"({specimen_display})")

    ax.set_title(" ".join(title_parts), fontsize=14)

    # Add sample size annotations if requested
    if show_n:
        y_lim = ax.get_ylim()
        if is_ct:
            # For inverted CT axis
            y_range = abs(y_lim[1] - y_lim[0])
            y_pos = y_lim[1] + y_range * 0.08
        else:
            # For log scale
            log_range = np.log10(y_lim[1]) - np.log10(y_lim[0])
            y_pos = 10 ** (np.log10(y_lim[1]) - log_range * 0.08)

        for bin_val, data in zip(unique_bins, plot_data):
            ax.annotate(
                f"n={len(data)}",
                xy=(bin_val, y_pos),
                fontsize=8,
                ha="center",
                va="bottom",
                alpha=0.7,
                rotation=90,
            )

    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    plt.close(fig)
    return fig


def plot_detection_probability(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    time_bin_size: float = 1.0,
    time_range: tuple[float, float] | None = None,
    min_observations: int = 3,
    figsize: tuple[int, int] = (10, 6),
    line_color: str = "steelblue",
    show_ci: bool = True,
    ci_alpha: float = 0.3,
    show_n: bool = True,
) -> Figure:
    """
    Plot detection probability (proportion of positive measurements) over time.

    Creates a line plot showing the probability of detecting a positive measurement
    at each time bin, with optional 95% confidence intervals. Useful for understanding
    when during infection a biomarker is most likely to be detected.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, uses first biomarker found.
        specimen: Optional filter for a specific specimen type. If None, uses first specimen found.
        time_bin_size: Size of time bins in days for aggregating measurements. Defaults to 1.0.
        time_range: Optional tuple (min_time, max_time) to limit the time axis.
            If None, uses the full range of data.
        min_observations: Minimum number of observations required per time bin to include
            in the plot. Bins with fewer observations are excluded. Defaults to 3.
        figsize: Figure size as (width, height). Defaults to (10, 6).
        line_color: Color for the probability line and confidence band. Defaults to "steelblue".
        show_ci: If True, shows 95% confidence interval band using Wilson score interval.
            Defaults to True.
        ci_alpha: Transparency of the confidence interval band (0-1). Defaults to 0.3.
        show_n: If True, shows the number of observations at each time point. Defaults to True.

    Returns:
        matplotlib.figure.Figure: The generated figure containing the detection probability plot.

    Raises:
        ValueError: If dataset is missing required keys, is empty, or has no valid data.

    Note:
        - Detection probability is calculated as the proportion of non-negative measurements
        - Confidence intervals use the Wilson score interval for binomial proportions
        - Y-axis is fixed to 0-1 range
    """
    # Validate input
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "dataset_id"]
    missing_keys = [key for key in required_keys if key not in dataset]
    if missing_keys:
        raise ValueError(f"Dataset missing required keys: {missing_keys}")

    if not dataset["participants"]:
        raise ValueError("Dataset has no participants")

    # Extract time series data from raw dataset
    time_series_data = []
    for participant_id, participant in enumerate(dataset["participants"], 1):
        measurements = participant.get("measurements", [])
        for measurement in measurements:
            analyte_name = measurement.get("analyte")
            time_series_data.append(
                {
                    "participant_id": participant_id,
                    "time": measurement.get("time"),
                    "value": measurement.get("value"),
                    "analyte": analyte_name,
                }
            )

    if not time_series_data:
        raise ValueError("Dataset has no measurements")

    df = pd.DataFrame(time_series_data)

    # Convert time to numeric, filtering out "unknown" values
    df = df[df["time"] != "unknown"].copy()
    df["time_num"] = pd.to_numeric(df["time"], errors="coerce")

    # Join with analyte metadata to get specimen and unit information
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "reference_event": analyte_info.get("reference_event"),
            "biomarker": analyte_info.get("biomarker"),
        }

    # Add metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
    )
    df["reference_event"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("reference_event")
    )
    df["biomarker"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("biomarker")
    )

    # Filter by biomarker if specified
    if biomarker is not None:
        df = df[df["biomarker"] == biomarker]
        if df.empty:
            raise ValueError(f"No measurements found for biomarker '{biomarker}'")

    # Filter by specimen if specified
    if specimen is not None:
        df = df[df["specimen"] == specimen]
        if df.empty:
            raise ValueError(f"No measurements found for specimen '{specimen}'")

    # Drop rows with NaN time
    df = df.dropna(subset=["time_num"])

    if df.empty:
        raise ValueError("No valid measurements found after filtering")

    # Apply time range filter if specified
    if time_range is not None:
        df = df[(df["time_num"] >= time_range[0]) & (df["time_num"] <= time_range[1])]
        if df.empty:
            raise ValueError(f"No measurements found in time range {time_range}")

    # Mark positive vs negative measurements
    df["is_positive"] = df["value"] != NEGATIVE_VALUE

    # Create time bins centered at integers (or multiples of bin_size)
    time_min = df["time_num"].min()
    time_max = df["time_num"].max()
    if time_range is not None:
        time_min, time_max = time_range

    # Round time_min down and time_max up to nearest multiple of bin_size to get bin centers
    center_min = np.floor(time_min / time_bin_size) * time_bin_size
    center_max = np.ceil(time_max / time_bin_size) * time_bin_size

    # Create bin edges at ±bin_size/2 around centers
    bins = np.arange(
        center_min - time_bin_size / 2, center_max + time_bin_size, time_bin_size
    )
    bin_centers = np.arange(center_min, center_max + time_bin_size / 2, time_bin_size)

    df["time_bin"] = pd.cut(
        df["time_num"], bins=bins, labels=bin_centers.astype(int), include_lowest=True
    )
    df["time_bin_num"] = df["time_bin"].astype(float)

    # Calculate detection probability per time bin
    def calc_detection_stats(group):
        n = len(group)
        n_positive = group["is_positive"].sum()
        p = n_positive / n if n > 0 else 0

        # Wilson score interval for 95% CI
        z = 1.96
        if n > 0:
            denominator = 1 + z**2 / n
            center = (p + z**2 / (2 * n)) / denominator
            margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denominator
            ci_lower = max(0, center - margin)
            ci_upper = min(1, center + margin)
        else:
            ci_lower = 0
            ci_upper = 0

        return pd.Series(
            {
                "n": n,
                "n_positive": n_positive,
                "probability": p,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
            }
        )

    stats_df = (
        df.groupby("time_bin_num", observed=False)
        .apply(calc_detection_stats)
        .reset_index()
    )

    # Filter by min_observations
    stats_df = stats_df[stats_df["n"] >= min_observations]

    if stats_df.empty:
        raise ValueError(
            f"No time bins have at least {min_observations} observations. "
            "Try reducing min_observations or using a larger time_bin_size."
        )

    # Sort by time
    stats_df = stats_df.sort_values("time_bin_num")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot confidence interval band
    if show_ci:
        ax.fill_between(
            stats_df["time_bin_num"],
            stats_df["ci_lower"],
            stats_df["ci_upper"],
            color=line_color,
            alpha=ci_alpha,
            label="95% CI",
        )

    # Plot probability line
    ax.plot(
        stats_df["time_bin_num"],
        stats_df["probability"],
        color=line_color,
        linewidth=2,
        marker="o",
        markersize=5,
        label="Detection probability",
    )

    # Set y-axis to 0-1 range with buffer for annotations
    ax.set_ylim(-0.05, 1.15 if show_n else 1.05)

    # Set x-axis ticks to show actual day values with reasonable spacing
    unique_bins = stats_df["time_bin_num"].tolist()
    n_bins = len(unique_bins)
    if n_bins <= 10:
        ax.set_xticks(unique_bins)
        ax.set_xticklabels([int(b) for b in unique_bins])
    else:
        n_ticks = min(10, n_bins)
        tick_indices = np.round(np.linspace(0, n_bins - 1, n_ticks)).astype(int)
        tick_positions = [unique_bins[i] for i in tick_indices]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([int(p) for p in tick_positions])

    # Labels and title
    reference_event = (
        df["reference_event"].dropna().iloc[0]
        if not df["reference_event"].dropna().empty
        else "reference"
    )
    specimen_display = (
        df["specimen"].dropna().iloc[0] if not df["specimen"].dropna().empty else ""
    )
    biomarker_display = (
        df["biomarker"].dropna().iloc[0] if not df["biomarker"].dropna().empty else ""
    )

    ax.set_xlabel(f"Time after {reference_event} (days)", fontsize=12)
    ax.set_ylabel("Detection probability", fontsize=12)

    # Title
    dataset_id = dataset.get("dataset_id", "Dataset")
    title_parts = [f"Detection Probability: {dataset_id}"]
    if biomarker_display:
        title_parts.append(f"({biomarker_display}")
        if specimen_display:
            title_parts[-1] += f" - {specimen_display})"
        else:
            title_parts[-1] += ")"
    elif specimen_display:
        title_parts.append(f"({specimen_display})")

    ax.set_title(" ".join(title_parts), fontsize=14)

    # Add sample size annotations if requested
    if show_n:
        y_pos = 1.02
        for _, row in stats_df.iterrows():
            ax.annotate(
                f"n={int(row['n'])}",
                xy=(row["time_bin_num"], y_pos),
                fontsize=8,
                ha="center",
                va="bottom",
                alpha=0.7,
                rotation=90,
            )

    ax.legend(loc="best", fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.close(fig)
    return fig


def plot_clearance_curve(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    figsize: tuple[int, int] = (10, 6),
    line_color: str = "steelblue",
    show_ci: bool = True,
    ci_alpha: float = 0.3,
    show_censored: bool = True,
    show_n_at_risk: bool = True,
) -> Figure:
    """
    Plot Kaplan-Meier style clearance curve showing proportion still shedding over time.

    Creates a survival curve showing the proportion of participants who are still
    shedding at each time point. Clearance is defined as the time of the last
    positive measurement for each participant.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, uses first biomarker found.
        specimen: Optional filter for a specific specimen type. If None, uses first specimen found.
        figsize: Figure size as (width, height). Defaults to (10, 6).
        line_color: Color for the survival curve and confidence band. Defaults to "steelblue".
        show_ci: If True, shows 95% confidence interval band using Greenwood's formula.
            Defaults to True.
        ci_alpha: Transparency of the confidence interval band (0-1). Defaults to 0.3.
        show_censored: If True, shows tick marks for censored observations (participants
            whose last measurement was positive). Defaults to True.
        show_n_at_risk: If True, shows the number of participants at risk at key time points.
            Defaults to True.

    Returns:
        matplotlib.figure.Figure: The generated figure containing the clearance curve.

    Raises:
        ValueError: If dataset is missing required keys, is empty, or has no valid data.

    Note:
        - Clearance time is the time of last positive measurement for each participant
        - Participants whose last measurement is positive are treated as censored
        - Confidence intervals use Greenwood's formula for variance estimation
    """
    # Validate input
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "dataset_id"]
    missing_keys = [key for key in required_keys if key not in dataset]
    if missing_keys:
        raise ValueError(f"Dataset missing required keys: {missing_keys}")

    if not dataset["participants"]:
        raise ValueError("Dataset has no participants")

    # Extract time series data from raw dataset
    time_series_data = []
    for participant_id, participant in enumerate(dataset["participants"], 1):
        measurements = participant.get("measurements", [])
        for measurement in measurements:
            analyte_name = measurement.get("analyte")
            time_series_data.append(
                {
                    "participant_id": participant_id,
                    "time": measurement.get("time"),
                    "value": measurement.get("value"),
                    "analyte": analyte_name,
                }
            )

    if not time_series_data:
        raise ValueError("Dataset has no measurements")

    df = pd.DataFrame(time_series_data)

    # Convert time to numeric, filtering out "unknown" values
    df = df[df["time"] != "unknown"].copy()
    df["time_num"] = pd.to_numeric(df["time"], errors="coerce")

    # Join with analyte metadata to get specimen information
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "reference_event": analyte_info.get("reference_event"),
            "biomarker": analyte_info.get("biomarker"),
        }

    # Add metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
    )
    df["reference_event"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("reference_event")
    )
    df["biomarker"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("biomarker")
    )

    # Filter by biomarker if specified
    if biomarker is not None:
        df = df[df["biomarker"] == biomarker]
        if df.empty:
            raise ValueError(f"No measurements found for biomarker '{biomarker}'")

    # Filter by specimen if specified
    if specimen is not None:
        df = df[df["specimen"] == specimen]
        if df.empty:
            raise ValueError(f"No measurements found for specimen '{specimen}'")

    # Drop rows with NaN time
    df = df.dropna(subset=["time_num"])

    if df.empty:
        raise ValueError("No valid measurements found after filtering")

    # Mark positive vs negative measurements
    df["is_positive"] = df["value"] != NEGATIVE_VALUE

    # Calculate clearance time and censoring status for each participant
    clearance_data = []
    for participant_id, participant_df in df.groupby("participant_id"):
        participant_df = participant_df.sort_values("time_num")

        # Get positive measurements only
        positive_df = participant_df[participant_df["is_positive"]]

        if positive_df.empty:
            # No positive measurements - skip this participant
            continue

        # Check if censored: last measurement overall is positive
        last_measurement = participant_df.iloc[-1]
        censored = last_measurement["is_positive"]

        if censored:
            # Censored: use time of last positive measurement (no future measurements)
            event_time = positive_df["time_num"].max()
        else:
            # Cleared: use time of negative observation (first negative after last positive)
            last_positive_time = positive_df["time_num"].max()
            negative_df = participant_df[~participant_df["is_positive"]]
            # Get the first negative measurement after the last positive
            negatives_after = negative_df[negative_df["time_num"] > last_positive_time]
            if not negatives_after.empty:
                event_time = negatives_after["time_num"].min()
            else:
                # Fallback: use the last negative measurement
                event_time = negative_df["time_num"].max()

        clearance_data.append(
            {
                "participant_id": participant_id,
                "clearance_time": event_time,
                "censored": censored,
            }
        )

    if not clearance_data:
        raise ValueError("No participants with positive measurements found")

    clearance_df = pd.DataFrame(clearance_data)

    # Sort by clearance time
    clearance_df = clearance_df.sort_values("clearance_time")

    # Calculate Kaplan-Meier survival estimates
    n_total = len(clearance_df)
    times = [0]  # Start at time 0 with S(0) = 1
    survival = [1.0]
    ci_lower = [1.0]
    ci_upper = [1.0]
    n_at_risk_times = [0]
    n_at_risk_values = [n_total]

    # Censored times and survival values for plotting
    censored_times = []
    censored_survival = []

    # Group events by time
    event_times = (
        clearance_df.groupby("clearance_time")
        .agg(
            {
                "censored": lambda x: (~x).sum(),  # Number of events (cleared)
                "participant_id": "count",  # Total at this time
            }
        )
        .rename(columns={"censored": "events", "participant_id": "total"})
    )

    n_at_risk = n_total
    cumulative_var = 0.0
    current_survival = 1.0

    for time, row in event_times.iterrows():
        d = row["events"]  # Number who cleared (not censored)
        c = row["total"] - d  # Number censored at this time

        if d > 0 and n_at_risk > 0:
            # Kaplan-Meier update
            current_survival = current_survival * (1 - d / n_at_risk)

            # Greenwood's formula for variance
            if n_at_risk > d:
                cumulative_var += d / (n_at_risk * (n_at_risk - d))

            times.append(time)
            survival.append(current_survival)

            # Calculate CI
            se = current_survival * np.sqrt(cumulative_var) if cumulative_var > 0 else 0
            ci_lower.append(max(0, current_survival - 1.96 * se))
            ci_upper.append(min(1, current_survival + 1.96 * se))

            n_at_risk_times.append(time)
            n_at_risk_values.append(n_at_risk - d - c)

        # Record censored observations
        if c > 0:
            # Get the survival value at this time for plotting censored marks
            censored_times.extend([time] * c)
            censored_survival.extend([current_survival] * c)

        n_at_risk -= d + c

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot confidence interval band
    if show_ci and len(times) > 1:
        # Create step-filled CI band
        ax.fill_between(
            times,
            ci_lower,
            ci_upper,
            step="post",
            color=line_color,
            alpha=ci_alpha,
            label="95% CI",
        )

    # Plot survival curve as step function
    ax.step(
        times,
        survival,
        where="post",
        color=line_color,
        linewidth=2,
        label="Clearance curve",
    )

    # Plot censored observations
    if show_censored and censored_times:
        ax.scatter(
            censored_times,
            censored_survival,
            marker="|",
            s=100,
            color=line_color,
            zorder=3,
            label="Censored",
        )

    # Set y-axis to 0-1 range with buffer
    ax.set_ylim(-0.05, 1.05)

    # Labels and title
    reference_event = (
        df["reference_event"].dropna().iloc[0]
        if not df["reference_event"].dropna().empty
        else "reference"
    )
    specimen_display = (
        df["specimen"].dropna().iloc[0] if not df["specimen"].dropna().empty else ""
    )
    biomarker_display = (
        df["biomarker"].dropna().iloc[0] if not df["biomarker"].dropna().empty else ""
    )

    ax.set_xlabel(f"Time after {reference_event} (days)", fontsize=12)
    ax.set_ylabel("Proportion still shedding", fontsize=12)

    # Title
    dataset_id = dataset.get("dataset_id", "Dataset")
    title_parts = [f"Clearance Curve: {dataset_id}"]
    if biomarker_display:
        title_parts.append(f"({biomarker_display}")
        if specimen_display:
            title_parts[-1] += f" - {specimen_display})"
        else:
            title_parts[-1] += ")"
    elif specimen_display:
        title_parts.append(f"({specimen_display})")

    ax.set_title(" ".join(title_parts), fontsize=14)

    # Add number at risk annotations
    if show_n_at_risk and n_at_risk_times:
        # Select a subset of time points for display
        n_points = min(6, len(n_at_risk_times))
        indices = np.round(np.linspace(0, len(n_at_risk_times) - 1, n_points)).astype(
            int
        )

        # Add text below the plot
        risk_text = "At risk: "
        for i in indices:
            t = n_at_risk_times[i]
            n = n_at_risk_values[i]
            risk_text += f"  t={int(t)}: {n}"

        ax.text(
            0.5,
            -0.12,
            risk_text,
            transform=ax.transAxes,
            fontsize=9,
            ha="center",
            va="top",
        )

    ax.legend(loc="best", fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.close(fig)
    return fig
