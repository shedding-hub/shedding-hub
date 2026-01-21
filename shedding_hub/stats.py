import warnings
import pandas as pd
import numpy as np
from typing import Dict, Any, Literal

# Constants
NEGATIVE_VALUE = "negative"


def calc_shedding_summary(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
) -> pd.DataFrame:
    """
    Calculate per-participant shedding summary statistics.

    Creates a summary table with key shedding metrics for each participant,
    including timing, duration, peak values, and clearance status.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, uses all biomarkers.
        specimen: Optional filter for a specific specimen type. If None, uses all specimens.

    Returns:
        pandas.DataFrame with columns:
            - participant_id: Participant identifier
            - biomarker: Biomarker name
            - specimen: Specimen type
            - value_type: Type of value ('concentration' or 'ct')
            - reference_event: Reference event for time measurements
            - first_positive_time: Time of first positive measurement
            - last_positive_time: Time of last positive measurement
            - shedding_duration: Duration from first to last positive (days)
            - peak_value: Maximum measurement value
            - peak_time: Time of peak measurement
            - n_positive: Number of positive measurements
            - n_negative: Number of negative measurements
            - n_total: Total number of measurements
            - clearance_status: 'cleared' if last measurement is negative, 'censored' otherwise
            - clearance_time: Time of clearance (first negative after last positive) or censoring

    Raises:
        ValueError: If dataset is missing required keys, is empty, or has no valid data.

    Example:
        >>> import shedding_hub as sh
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> from shedding_hub.stats import calc_shedding_summary
        >>> summary = calc_shedding_summary(data, specimen='sputum')
        >>> print(summary[['participant_id', 'shedding_duration', 'peak_value', 'clearance_status']])
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

    # Helper function to check if unit is CT value
    def _is_ct_value(unit: str | None) -> bool:
        if unit is None:
            return False
        unit_lower = str(unit).lower()
        return "ct" in unit_lower or "cycle" in unit_lower

    # Join with analyte metadata
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        unit = analyte_info.get("unit")
        value_type = "ct" if _is_ct_value(unit) else "concentration"

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "reference_event": analyte_info.get("reference_event"),
            "biomarker": analyte_info.get("biomarker"),
            "value_type": value_type,
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
    df["value_type"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("value_type")
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

    # Convert value to numeric for positive measurements
    df["value_num"] = pd.to_numeric(df["value"], errors="coerce")

    # Calculate summary statistics for each participant and analyte combination
    summary_data = []

    for (participant_id, analyte), group in df.groupby(["participant_id", "analyte"]):
        group = group.sort_values("time_num")

        # Get metadata
        biomarker_val = group["biomarker"].iloc[0]
        specimen_val = group["specimen"].iloc[0]
        reference_event = group["reference_event"].iloc[0]
        value_type_val = group["value_type"].iloc[0]

        # Positive and negative measurements
        positive_df = group[group["is_positive"]]
        negative_df = group[~group["is_positive"]]

        n_positive = len(positive_df)
        n_negative = len(negative_df)
        n_total = len(group)

        # Skip if no positive measurements
        if positive_df.empty:
            summary_data.append(
                {
                    "participant_id": participant_id,
                    "biomarker": biomarker_val,
                    "specimen": specimen_val,
                    "value_type": value_type_val,
                    "reference_event": reference_event,
                    "first_positive_time": np.nan,
                    "last_positive_time": np.nan,
                    "shedding_duration": np.nan,
                    "peak_value": np.nan,
                    "peak_time": np.nan,
                    "n_positive": n_positive,
                    "n_negative": n_negative,
                    "n_total": n_total,
                    "clearance_status": "no_positive",
                    "clearance_time": np.nan,
                }
            )
            continue

        # First and last positive times
        first_positive_time = positive_df["time_num"].min()
        last_positive_time = positive_df["time_num"].max()
        shedding_duration = last_positive_time - first_positive_time

        # Peak value and time (excluding NaN values)
        positive_with_values = positive_df.dropna(subset=["value_num"])
        if not positive_with_values.empty:
            peak_idx = positive_with_values["value_num"].idxmax()
            peak_value = positive_with_values.loc[peak_idx, "value_num"]
            peak_time = positive_with_values.loc[peak_idx, "time_num"]
        else:
            peak_value = np.nan
            peak_time = np.nan

        # Clearance status and time
        last_measurement = group.iloc[-1]
        if last_measurement["is_positive"]:
            # Censored: last measurement is positive
            clearance_status = "censored"
            clearance_time = last_positive_time
        else:
            # Cleared: find first negative after last positive
            clearance_status = "cleared"
            negatives_after = negative_df[negative_df["time_num"] > last_positive_time]
            if not negatives_after.empty:
                clearance_time = negatives_after["time_num"].min()
            else:
                clearance_time = negative_df["time_num"].max()

        summary_data.append(
            {
                "participant_id": participant_id,
                "biomarker": biomarker_val,
                "specimen": specimen_val,
                "value_type": value_type_val,
                "reference_event": reference_event,
                "first_positive_time": first_positive_time,
                "last_positive_time": last_positive_time,
                "shedding_duration": shedding_duration,
                "peak_value": peak_value,
                "peak_time": peak_time,
                "n_positive": n_positive,
                "n_negative": n_negative,
                "n_total": n_total,
                "clearance_status": clearance_status,
                "clearance_time": clearance_time,
            }
        )

    if not summary_data:
        raise ValueError("No valid participant data found")

    return pd.DataFrame(summary_data)


def calc_detection_summary(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    time_bin_size: float = 1.0,
    time_range: tuple[float, float] | None = None,
    min_observations: int = 1,
) -> pd.DataFrame:
    """
    Calculate detection rate statistics by time bin.

    Creates a summary table showing the proportion of positive measurements
    at each time bin, with 95% confidence intervals using Wilson score interval.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, uses all biomarkers.
        specimen: Optional filter for a specific specimen type. If None, uses all specimens.
        time_bin_size: Size of time bins in days. Defaults to 1.0.
        time_range: Optional tuple (min_time, max_time) to limit the time range.
            If None, uses the full range of data.
        min_observations: Minimum number of observations required per time bin.
            Bins with fewer observations are excluded. Defaults to 1.

    Returns:
        pandas.DataFrame with columns:
            - time: Center of the time bin (days)
            - n_tested: Number of measurements in this bin
            - n_positive: Number of positive measurements
            - n_negative: Number of negative measurements
            - proportion: Proportion of positive measurements (0-1)
            - ci_lower: Lower bound of 95% CI (Wilson score interval)
            - ci_upper: Upper bound of 95% CI (Wilson score interval)

    Raises:
        ValueError: If dataset is missing required keys, is empty, or has no valid data.

    Example:
        >>> import shedding_hub as sh
        >>> from shedding_hub.stats import calc_detection_summary
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> detection = calc_detection_summary(data, specimen='sputum', time_bin_size=7)
        >>> print(detection[['time', 'n_tested', 'n_positive', 'proportion']])
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

    # Join with analyte metadata
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "biomarker": analyte_info.get("biomarker"),
        }

    # Add metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
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

    # Round time_min down and time_max up to nearest multiple of bin_size
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

    # Calculate detection statistics per time bin
    summary_data = []

    for time_bin, group in df.groupby("time_bin_num"):
        n_tested = len(group)

        # Skip bins with insufficient observations
        if n_tested < min_observations:
            continue

        n_positive = group["is_positive"].sum()
        n_negative = n_tested - n_positive
        proportion = n_positive / n_tested

        # Wilson score interval for 95% CI
        z = 1.96
        denominator = 1 + z**2 / n_tested
        center = (proportion + z**2 / (2 * n_tested)) / denominator
        margin = (
            z
            * np.sqrt(
                proportion * (1 - proportion) / n_tested + z**2 / (4 * n_tested**2)
            )
            / denominator
        )
        ci_lower = max(0, center - margin)
        ci_upper = min(1, center + margin)

        summary_data.append(
            {
                "time": int(time_bin),
                "n_tested": n_tested,
                "n_positive": n_positive,
                "n_negative": n_negative,
                "proportion": proportion,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
            }
        )

    if not summary_data:
        raise ValueError(
            f"No time bins have at least {min_observations} observations. "
            "Try reducing min_observations or using a larger time_bin_size."
        )

    result_df = pd.DataFrame(summary_data)
    return result_df.sort_values("time").reset_index(drop=True)


def calc_clearance_summary(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    time_points: list[float] | None = None,
) -> Dict[str, Any]:
    """
    Calculate Kaplan-Meier clearance statistics.

    Provides survival analysis statistics including median time to clearance,
    proportion still shedding at specified time points, and number at risk.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, uses all biomarkers.
        specimen: Optional filter for a specific specimen type. If None, uses all specimens.
        time_points: List of time points (days) at which to report proportion still shedding.
            If None, defaults to [7, 14, 21, 28].

    Returns:
        Dictionary containing:
            - 'n_participants': Total number of participants analyzed
            - 'n_cleared': Number who cleared (observed event)
            - 'n_censored': Number censored (still shedding at last observation)
            - 'median_clearance_time': Median time to clearance (None if <50% cleared)
            - 'median_ci_lower': Lower bound of 95% CI for median
            - 'median_ci_upper': Upper bound of 95% CI for median
            - 'survival_table': DataFrame with Kaplan-Meier estimates at each event time
            - 'time_point_summary': DataFrame with proportion still shedding at specified time points

    Raises:
        ValueError: If dataset is missing required keys, is empty, or has no valid data.

    Example:
        >>> import shedding_hub as sh
        >>> from shedding_hub.stats import calc_clearance_summary
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> summary = calc_clearance_summary(data, specimen='sputum')
        >>> print(f"Median clearance: {summary['median_clearance_time']} days")
        >>> print(summary['time_point_summary'])
    """
    if time_points is None:
        time_points = [7, 14, 21, 28]

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

    # Join with analyte metadata
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "biomarker": analyte_info.get("biomarker"),
        }

    # Add metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
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
            negatives_after = negative_df[negative_df["time_num"] > last_positive_time]
            if not negatives_after.empty:
                event_time = negatives_after["time_num"].min()
            else:
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
    clearance_df = clearance_df.sort_values("clearance_time")

    n_participants = len(clearance_df)
    n_censored = clearance_df["censored"].sum()
    n_cleared = n_participants - n_censored

    # Calculate Kaplan-Meier survival estimates
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

    survival_data = []
    n_at_risk = n_participants
    cumulative_var = 0.0
    current_survival = 1.0

    # Add initial point
    survival_data.append(
        {
            "time": 0,
            "n_at_risk": n_at_risk,
            "n_events": 0,
            "n_censored": 0,
            "survival": 1.0,
            "ci_lower": 1.0,
            "ci_upper": 1.0,
        }
    )

    for time, row in event_times.iterrows():
        d = row["events"]  # Number who cleared
        c = row["total"] - d  # Number censored at this time

        if d > 0 and n_at_risk > 0:
            # Kaplan-Meier update
            current_survival = current_survival * (1 - d / n_at_risk)

            # Greenwood's formula for variance
            if n_at_risk > d:
                cumulative_var += d / (n_at_risk * (n_at_risk - d))

            # Calculate CI
            se = current_survival * np.sqrt(cumulative_var) if cumulative_var > 0 else 0
            ci_lower = max(0, current_survival - 1.96 * se)
            ci_upper = min(1, current_survival + 1.96 * se)

            survival_data.append(
                {
                    "time": time,
                    "n_at_risk": n_at_risk,
                    "n_events": d,
                    "n_censored": c,
                    "survival": current_survival,
                    "ci_lower": ci_lower,
                    "ci_upper": ci_upper,
                }
            )

        n_at_risk -= d + c

    survival_table = pd.DataFrame(survival_data)

    # Calculate median clearance time (time when survival <= 0.5)
    median_clearance_time = None
    median_ci_lower = None
    median_ci_upper = None

    if current_survival <= 0.5:
        # Find the first time survival drops to or below 0.5
        below_median = survival_table[survival_table["survival"] <= 0.5]
        if not below_median.empty:
            median_clearance_time = below_median["time"].iloc[0]
            median_ci_lower = below_median["ci_lower"].iloc[0]
            median_ci_upper = below_median["ci_upper"].iloc[0]

    # Calculate survival at specified time points
    time_point_data = []
    for t in time_points:
        # Find survival at time t (last survival estimate <= t)
        at_or_before = survival_table[survival_table["time"] <= t]
        if not at_or_before.empty:
            row = at_or_before.iloc[-1]
            survival_at_t = row["survival"]
            ci_lower_at_t = row["ci_lower"]
            ci_upper_at_t = row["ci_upper"]
            n_at_risk_at_t = row["n_at_risk"]
        else:
            survival_at_t = 1.0
            ci_lower_at_t = 1.0
            ci_upper_at_t = 1.0
            n_at_risk_at_t = n_participants

        time_point_data.append(
            {
                "time": t,
                "proportion_shedding": survival_at_t,
                "proportion_cleared": 1 - survival_at_t,
                "ci_lower": ci_lower_at_t,
                "ci_upper": ci_upper_at_t,
                "n_at_risk": n_at_risk_at_t,
            }
        )

    time_point_summary = pd.DataFrame(time_point_data)

    return {
        "n_participants": n_participants,
        "n_cleared": n_cleared,
        "n_censored": n_censored,
        "median_clearance_time": median_clearance_time,
        "median_ci_lower": median_ci_lower,
        "median_ci_upper": median_ci_upper,
        "survival_table": survival_table,
        "time_point_summary": time_point_summary,
    }


def calc_value_summary(
    dataset: Dict[str, Any],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    value: str | None = None,
    time_bin_size: float = 1.0,
    time_range: tuple[float, float] | None = None,
    min_observations: int = 1,
) -> pd.DataFrame:
    """
    Calculate summary statistics of measurement values by time bin.

    Creates a summary table showing mean, median, standard deviation, and
    interquartile range of measurement values at each time bin. This is
    the tabular version of plot_mean_trajectory.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.
        biomarker: Optional filter for a specific biomarker. If None, uses all biomarkers.
        specimen: Optional filter for a specific specimen type. If None, uses all specimens.
        value: Optional filter for value type. Options are "concentration" or "ct".
            If None, uses all data (may raise error if mixed). Defaults to None.
        time_bin_size: Size of time bins in days. Defaults to 1.0.
        time_range: Optional tuple (min_time, max_time) to limit the time range.
            If None, uses the full range of data.
        min_observations: Minimum number of observations required per time bin.
            Bins with fewer observations are excluded. Defaults to 1.

    Returns:
        pandas.DataFrame with columns:
            - time: Center of the time bin (days)
            - n: Number of observations
            - mean: Mean value
            - std: Standard deviation
            - median: Median value
            - q25: 25th percentile
            - q75: 75th percentile
            - min: Minimum value
            - max: Maximum value

    Raises:
        ValueError: If dataset is missing required keys, is empty, has no valid data,
            or contains mixed CT and concentration values without filtering.

    Example:
        >>> import shedding_hub as sh
        >>> from shedding_hub.stats import calc_value_summary
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> summary = calc_value_summary(data, specimen='sputum', time_bin_size=7)
        >>> print(summary[['time', 'n', 'mean', 'median']])
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

    # Join with analyte metadata
    analyte_metadata = {}
    for analyte_name, analyte_info in dataset["analytes"].items():
        specimen_value = analyte_info.get("specimen")
        if isinstance(specimen_value, list):
            specimen_value = "+".join(specimen_value)

        analyte_metadata[analyte_name] = {
            "specimen": specimen_value,
            "unit": analyte_info.get("unit"),
            "biomarker": analyte_info.get("biomarker"),
        }

    # Add metadata to DataFrame
    df["specimen"] = df["analyte"].map(
        lambda x: analyte_metadata.get(x, {}).get("specimen")
    )
    df["unit"] = df["analyte"].map(lambda x: analyte_metadata.get(x, {}).get("unit"))
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

    # Helper function to check if unit is CT value
    def _is_ct_value(unit: str | None) -> bool:
        if unit is None:
            return False
        unit_lower = str(unit).lower()
        return "ct" in unit_lower or "cycle" in unit_lower

    # Determine if each row is CT value or concentration
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

    # Exclude negative values
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

    # Round time_min down and time_max up to nearest multiple of bin_size
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

    # Calculate summary statistics per time bin
    summary_data = []

    for time_bin, group in df.groupby("time_bin_num"):
        values = group["value_num"]
        n = len(values)

        # Skip bins with insufficient observations
        if n < min_observations:
            continue

        summary_data.append(
            {
                "time": int(time_bin),
                "n": n,
                "mean": values.mean(),
                "std": values.std() if n > 1 else np.nan,
                "median": values.median(),
                "q25": values.quantile(0.25),
                "q75": values.quantile(0.75),
                "min": values.min(),
                "max": values.max(),
            }
        )

    if not summary_data:
        raise ValueError(
            f"No time bins have at least {min_observations} observations. "
            "Try reducing min_observations or using a larger time_bin_size."
        )

    result_df = pd.DataFrame(summary_data)
    return result_df.sort_values("time").reset_index(drop=True)


def calc_dataset_summary(
    dataset: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Summarize basic information about a dataset.

    Provides an overview of the dataset including biomarkers, specimens,
    value types, sample sizes, and time ranges.

    Args:
        dataset: Raw dataset dictionary from load_dataset() containing 'analytes',
            'participants', and 'dataset_id' keys.

    Returns:
        Dictionary containing:
            - 'dataset_id': Dataset identifier
            - 'title': Dataset title (if available)
            - 'doi': Dataset DOI (if available)
            - 'n_participants': Number of participants
            - 'n_measurements': Total number of measurements
            - 'n_analytes': Number of analytes
            - 'biomarkers': List of unique biomarkers
            - 'specimens': List of unique specimens
            - 'value_types': List of value types ('concentration' and/or 'ct')
            - 'reference_events': List of unique reference events
            - 'time_range': Tuple of (min_time, max_time) in days
            - 'analyte_details': DataFrame with details per analyte

    Raises:
        ValueError: If dataset is missing required keys or is empty.

    Example:
        >>> import shedding_hub as sh
        >>> from shedding_hub.stats import calc_dataset_summary
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> summary = calc_dataset_summary(data)
        >>> print(f"Participants: {summary['n_participants']}")
        >>> print(f"Biomarkers: {summary['biomarkers']}")
    """
    # Validate input
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")

    required_keys = ["analytes", "participants", "dataset_id"]
    missing_keys = [key for key in required_keys if key not in dataset]
    if missing_keys:
        raise ValueError(f"Dataset missing required keys: {missing_keys}")

    # Basic info
    dataset_id = dataset.get("dataset_id", "Unknown")
    title = dataset.get("title", None)
    doi = dataset.get("doi", dataset.get("url", None))

    n_participants = len(dataset["participants"]) if dataset["participants"] else 0

    # Helper function to check if unit is CT value
    def _is_ct_value(unit: str | None) -> bool:
        if unit is None:
            return False
        unit_lower = str(unit).lower()
        return "ct" in unit_lower or "cycle" in unit_lower

    # Extract analyte details
    analyte_details = []
    biomarkers = set()
    specimens = set()
    reference_events = set()
    value_types = set()

    for analyte_name, analyte_info in dataset["analytes"].items():
        biomarker = analyte_info.get("biomarker")
        specimen = analyte_info.get("specimen")
        if isinstance(specimen, list):
            specimen = "+".join(specimen)
        unit = analyte_info.get("unit")
        reference_event = analyte_info.get("reference_event")
        lod = analyte_info.get("limit_of_detection")
        loq = analyte_info.get("limit_of_quantification")

        if biomarker:
            biomarkers.add(biomarker)
        if specimen:
            specimens.add(specimen)
        if reference_event:
            reference_events.add(reference_event)

        is_ct = _is_ct_value(unit)
        value_type = "ct" if is_ct else "concentration"
        value_types.add(value_type)

        analyte_details.append(
            {
                "analyte": analyte_name,
                "biomarker": biomarker,
                "specimen": specimen,
                "unit": unit,
                "value_type": value_type,
                "reference_event": reference_event,
                "limit_of_detection": lod,
                "limit_of_quantification": loq,
            }
        )

    analyte_df = pd.DataFrame(analyte_details)

    # Count measurements and extract time range
    n_measurements = 0
    times = []
    n_positive = 0
    n_negative = 0

    for participant in dataset["participants"]:
        measurements = participant.get("measurements", [])
        n_measurements += len(measurements)

        for m in measurements:
            time_val = m.get("time")
            if time_val != "unknown":
                try:
                    times.append(float(time_val))
                except (ValueError, TypeError):
                    pass

            value = m.get("value")
            if value == NEGATIVE_VALUE:
                n_negative += 1
            else:
                n_positive += 1

    # Calculate time range
    if times:
        time_range = (min(times), max(times))
    else:
        time_range = (None, None)

    return {
        "dataset_id": dataset_id,
        "title": title,
        "doi": doi,
        "n_participants": n_participants,
        "n_measurements": n_measurements,
        "n_positive": n_positive,
        "n_negative": n_negative,
        "n_analytes": len(dataset["analytes"]),
        "biomarkers": sorted(list(biomarkers)),
        "specimens": sorted(list(specimens)),
        "value_types": sorted(list(value_types)),
        "reference_events": sorted(list(reference_events)),
        "time_range": time_range,
        "analyte_details": analyte_df,
    }


def compare_datasets(
    datasets: list[Dict[str, Any]],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    value: Literal["concentration", "ct"] | None = None,
) -> pd.DataFrame:
    """
    Compare key statistics across multiple datasets.

    Creates a side-by-side comparison table of key shedding metrics
    across multiple studies.

    Args:
        datasets: List of dataset dictionaries from load_dataset().
        biomarker: Optional filter for a specific biomarker. If None and multiple
            biomarkers exist, a warning is issued.
        specimen: Optional filter for a specific specimen type. If None and multiple
            specimens exist, a warning is issued.
        value: Filter for value type, either 'concentration' or 'ct'. Peak values
            will only be calculated for analytes matching this value type. If None
            and both value types exist, peak statistics are excluded from results.

    Returns:
        pandas.DataFrame with one row per dataset and columns:
            - dataset_id: Dataset identifier
            - n_participants: Number of participants
            - n_measurements: Total measurements
            - pct_positive: Percentage of positive measurements
            - median_shedding_duration: Median duration from first to last positive (days)
            - iqr_shedding_duration: IQR of shedding duration (Q25-Q75)
            - median_peak_value: Median peak value across participants (if value type specified)
            - iqr_peak_value: IQR of peak values (Q25-Q75)
            - median_peak_time: Median time to peak (days)
            - pct_cleared: Percentage of participants who cleared
            - median_clearance_time: Median time to clearance (days, if >50% cleared)

    Raises:
        ValueError: If datasets is empty or contains invalid datasets.

    Example:
        >>> import shedding_hub as sh
        >>> from shedding_hub.stats import compare_datasets
        >>> data1 = sh.load_dataset('woelfel2020virological', local='./data')
        >>> data2 = sh.load_dataset('young2020epidemiologic', local='./data')
        >>> comparison = compare_datasets([data1, data2], specimen='sputum', value='concentration')
        >>> print(comparison)
    """
    if not datasets:
        raise ValueError("datasets list cannot be empty")

    # Helper function to check if unit is CT value
    def _is_ct_value(unit: str | None) -> bool:
        if unit is None:
            return False
        unit_lower = str(unit).lower()
        return "ct" in unit_lower or "cycle" in unit_lower

    # Check for multiple biomarkers/specimens across all datasets
    all_biomarkers = set()
    all_specimens = set()
    all_value_types = set()

    for dataset in datasets:
        if not dataset or not isinstance(dataset, dict):
            continue
        for analyte_info in dataset.get("analytes", {}).values():
            if analyte_info.get("biomarker"):
                all_biomarkers.add(analyte_info["biomarker"])
            specimen_val = analyte_info.get("specimen")
            if specimen_val:
                if isinstance(specimen_val, list):
                    specimen_val = "+".join(specimen_val)
                all_specimens.add(specimen_val)
            unit = analyte_info.get("unit")
            value_type = "ct" if _is_ct_value(unit) else "concentration"
            all_value_types.add(value_type)

    # Warn if multiple biomarkers/specimens exist without filtering
    if len(all_biomarkers) > 1 and biomarker is None:
        warnings.warn(
            f"Multiple biomarkers found ({', '.join(sorted(all_biomarkers))}) "
            "but no biomarker filter specified. Results will merge data from all biomarkers.",
            UserWarning,
            stacklevel=2,
        )

    if len(all_specimens) > 1 and specimen is None:
        warnings.warn(
            f"Multiple specimens found ({', '.join(sorted(all_specimens))}) "
            "but no specimen filter specified. Results will merge data from all specimens.",
            UserWarning,
            stacklevel=2,
        )

    # Determine if we can compute peak values
    compute_peak_values = True
    if len(all_value_types) > 1 and value is None:
        warnings.warn(
            "Both concentration and CT value types found but no value filter specified. "
            "Peak value statistics will be excluded from results.",
            UserWarning,
            stacklevel=2,
        )
        compute_peak_values = False

    comparison_data = []

    for dataset in datasets:
        if not dataset or not isinstance(dataset, dict):
            raise ValueError("Each dataset must be a non-empty dictionary")

        required_keys = ["analytes", "participants", "dataset_id"]
        missing_keys = [key for key in required_keys if key not in dataset]
        if missing_keys:
            raise ValueError(f"Dataset missing required keys: {missing_keys}")

        dataset_id = dataset.get("dataset_id", "Unknown")

        try:
            # Get shedding summary for this dataset (now includes value_type column)
            shedding_summary = calc_shedding_summary(
                dataset, biomarker=biomarker, specimen=specimen
            )

            # Filter to participants with positive measurements
            valid_summary = shedding_summary[
                shedding_summary["clearance_status"] != "no_positive"
            ]

            if valid_summary.empty:
                # No valid data for this dataset with given filters
                comparison_data.append(
                    {
                        "dataset_id": dataset_id,
                        "n_participants": 0,
                        "n_measurements": np.nan,
                        "pct_positive": np.nan,
                        "median_shedding_duration": np.nan,
                        "iqr_shedding_duration": None,
                        "median_peak_value": np.nan,
                        "iqr_peak_value": None,
                        "median_peak_time": np.nan,
                        "pct_cleared": np.nan,
                        "median_clearance_time": np.nan,
                    }
                )
                continue

            n_participants = len(valid_summary)

            # Calculate total measurements
            n_measurements = valid_summary["n_total"].sum()
            n_positive = valid_summary["n_positive"].sum()
            pct_positive = (
                (n_positive / n_measurements * 100) if n_measurements > 0 else np.nan
            )

            # Shedding duration statistics
            durations = valid_summary["shedding_duration"].dropna()
            if not durations.empty:
                median_duration = durations.median()
                q25_duration = durations.quantile(0.25)
                q75_duration = durations.quantile(0.75)
                iqr_duration = f"{q25_duration:.1f}-{q75_duration:.1f}"
            else:
                median_duration = np.nan
                iqr_duration = None

            # Peak value statistics (filtered by value type if specified)
            if compute_peak_values:
                if value is not None:
                    # Filter to specified value type
                    peak_summary = valid_summary[valid_summary["value_type"] == value]
                else:
                    peak_summary = valid_summary

                peak_values = peak_summary["peak_value"].dropna()
                if not peak_values.empty:
                    median_peak = peak_values.median()
                    q25_peak = peak_values.quantile(0.25)
                    q75_peak = peak_values.quantile(0.75)
                    iqr_peak = f"{q25_peak:.2e}-{q75_peak:.2e}"
                else:
                    median_peak = np.nan
                    iqr_peak = None

                # Peak time statistics (also filtered by value type)
                peak_times = peak_summary["peak_time"].dropna()
                median_peak_time = (
                    peak_times.median() if not peak_times.empty else np.nan
                )
            else:
                median_peak = np.nan
                iqr_peak = None
                median_peak_time = np.nan

            # Clearance statistics
            n_cleared = (valid_summary["clearance_status"] == "cleared").sum()
            pct_cleared = (
                (n_cleared / n_participants * 100) if n_participants > 0 else np.nan
            )

            # Median clearance time (from clearance_summary if >50% cleared)
            if pct_cleared >= 50:
                try:
                    clearance_summary = calc_clearance_summary(
                        dataset, biomarker=biomarker, specimen=specimen
                    )
                    median_clearance_time = clearance_summary["median_clearance_time"]
                except (ValueError, KeyError):
                    median_clearance_time = np.nan
            else:
                median_clearance_time = np.nan

            comparison_data.append(
                {
                    "dataset_id": dataset_id,
                    "n_participants": n_participants,
                    "n_measurements": n_measurements,
                    "pct_positive": round(pct_positive, 1),
                    "median_shedding_duration": median_duration,
                    "iqr_shedding_duration": iqr_duration,
                    "median_peak_value": median_peak,
                    "iqr_peak_value": iqr_peak,
                    "median_peak_time": median_peak_time,
                    "pct_cleared": round(pct_cleared, 1),
                    "median_clearance_time": median_clearance_time,
                }
            )

        except ValueError as e:
            # Handle datasets that don't have matching data
            comparison_data.append(
                {
                    "dataset_id": dataset_id,
                    "n_participants": 0,
                    "n_measurements": np.nan,
                    "pct_positive": np.nan,
                    "median_shedding_duration": np.nan,
                    "iqr_shedding_duration": None,
                    "median_peak_value": np.nan,
                    "iqr_peak_value": None,
                    "median_peak_time": np.nan,
                    "pct_cleared": np.nan,
                    "median_clearance_time": np.nan,
                }
            )

    return pd.DataFrame(comparison_data)
