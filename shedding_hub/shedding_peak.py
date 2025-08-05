import shedding_hub as sh
import pandas as pd
from typing import Optional, Dict, Any, Iterable, Tuple, Union
import matplotlib.pyplot as plt
import glob
from pathlib import Path
from matplotlib import cm

def calc_shedding_peak(
    dataset: Dict[str, Any],
    *,
    plotting: bool = False,
    output: str, 
) -> pd.DataFrame:
    """
    Calcualte summary statistics for the shedding peak using a loaded dataset by the `load_dataset` function.

    Args:
        dataset: Loaded dataset using the `load_dataset` function.
        plotting: Create a plot for individual level of shedding peak by specimen type.
        output: Type of dataframe returned. 
            summary: Summary table of shedding peak (min, max, mean) by biomarker and specimen.
            individual: Individual entries (Not summarized)


    Returns:
        Dataframe based on selected output. 

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
    

    # initialize a pandas dataframe
    df_shedding_peak = pd.DataFrame(columns=["dataset_id","participant_ID",#"age","sex",
                                                 "analyte","n_sample","first_sample","last_sample","shedding_peak"])
    df_analyte = pd.DataFrame(columns=["analyte","specimen","biomarker","reference_event"])

    # extract analyte data from the standardized shedding data loaded
    for key in dataset["analytes"]:
        row_new = {"analyte":key,
                   "specimen":dataset["analytes"][key]["specimen"],
                   "biomarker":dataset["analytes"][key]["biomarker"],
                   "reference_event":dataset["analytes"][key]["reference_event"]}
        df_analyte.loc[len(df_analyte)] = row_new
    
    # extract participant and measurement data from the standardized shedding data loaded
    participant_counter = 0
    for item in dataset["participants"]:
        participant_counter += 1
        participant_ID = participant_counter
        #participant_age = item["attributes"]["age"]
        #participant_sex = item["attributes"]["sex"]
        for name, group in pd.DataFrame.from_dict(item["measurements"]).groupby("analyte"):
            # Filter numeric values for shedding peak calculation
            numeric_values = group[group["value"].apply(lambda x: isinstance(x, (int, float)))]
            
            # Check if there are any numeric values
            if len(numeric_values) > 0:
                shedding_peak = group.loc[numeric_values["value"].idxmax(), "time"]
            else:
                # If all values are non-numeric (e.g., "negative"), set to null
                shedding_peak = None
            
            row_new = {
                "dataset_id": dataset["dataset_id"],
                "participant_ID" : participant_ID,
                #"age" : participant_age,
                #"sex" : participant_sex,
                "analyte" : name,
                "n_sample" : group[group["value"].notna()]["time"].size,
                "first_sample": group[
                    (group["value"].notna()) & (group["time"] != "unknown")
                ]["time"].min(),
                "last_sample": group[
                    (group["value"].notna()) & (group["time"] != "unknown")
                ]["time"].max(),
                "shedding_peak": shedding_peak
            }
            df_shedding_peak.loc[len(df_shedding_peak)] = row_new

    # merge analyte information and drop analyte column
    df_shedding_peak = df_shedding_peak.merge(df_analyte, how="left", on="analyte").drop(columns=["analyte"])
    # concatenate list of specimen types to string;
    df_shedding_peak["specimen"] = df_shedding_peak["specimen"].apply(
        lambda x: ", ".join(map(str, x)) if isinstance(x, list) else str(x)
    )
    # calculate individual level shedding peak
    df_shedding_peak["first_sample"] = pd.to_numeric(
        df_shedding_peak["first_sample"], errors="coerce"
    )
    df_shedding_peak["last_sample"] = pd.to_numeric(
        df_shedding_peak["last_sample"], errors="coerce"
    )
    df_shedding_peak["shedding_peak"] = pd.to_numeric(df_shedding_peak["shedding_peak"], errors='coerce')
    
    #if plotting == True:
    #    plt_shedding = plot_shedding_peak(df_shedding_peak, df_ID = df_ID)
    #    plt_shedding.show()

    df_return = df_shedding_peak.groupby(
        ["dataset_id", "biomarker", "specimen", "reference_event"]
    ).agg(
        shedding_peak_min   = ("shedding_peak", "min"),
        shedding_peak_q25   = ("shedding_peak", lambda x: x.quantile(0.25)),
        shedding_peak_median= ("shedding_peak", "median"),
        shedding_peak_q75   = ("shedding_peak", lambda x: x.quantile(0.75)),
        shedding_peak_max   = ("shedding_peak", "max"),
        shedding_peak_mean  = ("shedding_peak", "mean"),
        n_sample            = ("n_sample", "sum"),
        n_participant       = ("participant_ID", "nunique"),
    ).reset_index()
    
    if output == "summary":
        return df_return
    elif output == "individual":
        return df_shedding_peak
    else:
        raise ValueError("`output` must be either 'summary' or 'individual'")
    




def plot_shedding_peak_individual(
    df_shedding_peak: pd.DataFrame,
    max_nparticipant: int = 30,
    random_seed: int = 12345,
):
    """
    Faceted error-bar plot of each participant’s shedding window and peak,
    grouped by specimen.

    Parameters
    ----------
    df_shedding_peak : pd.DataFrame
        Must contain columns:
            ['specimen', 'participant_ID',
             'first_sample', 'last_sample', 'shedding_peak']
    max_nparticipant : int, default 30
        Maximum rows sampled per specimen facet.
    random_seed : int, default 12345
        Seed for deterministic sampling.

    Returns
    -------
    (fig, axes) : tuple
        Matplotlib Figure and list of Axes for further tweaking.
    """
    # ------------------------------------------------------------------
    # Input validation -------------------------------------------------
    # ------------------------------------------------------------------
    required_cols = {
        "specimen",
        "participant_ID",
        "first_sample",
        "last_sample",
        "shedding_peak",
    }

    # 1) Correct type
    if not isinstance(df_shedding_peak, pd.DataFrame):
        raise ValueError(
            "df_shedding_peak must be a pandas DataFrame, "
            f"got {type(df_shedding_peak).__name__}."
        )

    # 2) Required columns present
    missing = required_cols - set(df_shedding_peak.columns)
    if missing:
        raise ValueError(
            "Input DataFrame is missing required column(s): "
            f"{', '.join(sorted(missing))}"
        )

    # 3) Numeric time columns and at least one non-NA row
    time_cols = ["first_sample", "last_sample", "shedding_peak"]
    try:
        df_shedding_peak[time_cols] = df_shedding_peak[time_cols].apply(
            pd.to_numeric, errors="raise"
        )
    except Exception as e:
        raise ValueError(
            "Columns first_sample, last_sample, and shedding_peak "
            "must be numeric (or coercible to numeric)."
        ) from e

    if df_shedding_peak.dropna(subset=time_cols).empty:
        raise ValueError(
            "After dropping NA rows in first_sample, last_sample, "
            "and shedding_peak, no data remain."
        )

    # ------------------------------------------------------------------
    # Plotting logic (unchanged except for moving df name to `df`) -----
    # ------------------------------------------------------------------
    df = (
        df_shedding_peak
        .dropna(subset=["first_sample", "last_sample", "shedding_peak"])
        .copy()
    )

    # limit rows per specimen for legibility
    df = (
        df.groupby("specimen", group_keys=False)
          .apply(
              lambda g: g.sample(
                  n=max_nparticipant, random_state=random_seed
              ) if len(g) > max_nparticipant else g
          )
    )

    # error-bar extents
    df["err_plus"]  = df["last_sample"] - df["first_sample"]
    df["err_minus"] = 0

    specimens = df["specimen"].unique()
    n = len(specimens)

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 6), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, spec in zip(axes, specimens):
        g = df[df["specimen"] == spec].sort_values("participant_ID")
        y = list(range(len(g)))[::-1]

        # horizontal error bars (shedding window)
        ax.errorbar(
            x=g["first_sample"],
            y=y,
            xerr=[g["err_minus"], g["err_plus"]],
            fmt="none",
            ecolor="gray",
            elinewidth=2,
            capsize=0,
            zorder=1,
        )

        # peak markers
        ax.scatter(
            g["shedding_peak"],
            y,
            marker="D",
            s=35,
            color="tab:red",
            zorder=2,
            label="Peak",
        )

        ax.set_yticks(y)
        ax.set_yticklabels("P" + g["participant_ID"].astype(str))
        ax.set_xlabel("Days from reference event")
        ax.set_title(spec)
        ax.grid(axis="x", alpha=0.3)

    axes[0].set_ylabel("Participant")
    axes[-1].legend(loc="upper right")
    fig.suptitle(
        "Shedding window (error bar) + peak, grouped by specimen",
        y=1.02,
        fontsize=14,
    )
    plt.tight_layout()
    return fig, axes





def plot_shedding_peak_summary(
    dataset_ids: Union[str, Iterable[str]],                      # ← NEW sole “loader” input
    *,
    selected_biomarker: str = "SARS-CoV-2",
    selected_reference_event: str = "symptom onset",
    selected_min_nparticipant: int = 5,
    x_axis_upper_limit: int | None = 50,
    ax: plt.Axes | None = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Load multiple datasets, build shedding-peak summaries, and plot them.

    Parameters
    ----------
    dataset_ids : iterable of str
        IDs understood by `sh.load_dataset`.
    selected_biomarker, selected_reference_event, selected_min_nparticipant,
    x_axis_upper_limit, ax : plotting/filter controls (unchanged).

    Returns
    -------
    (fig, ax) : Matplotlib Figure and Axes.
    """

    # ── 0. normalise input -------------------------------------------
    if isinstance(dataset_ids, str):
        dataset_ids = [dataset_ids]           # wrap single ID

    # ---------------- 1. load & summarise ---------------------------------
    summary_frames = []

    for ds_id in dataset_ids:
        try:
            df_raw  = sh.load_dataset(ds_id)
            df_peak = calc_shedding_peak(df_raw, output="summary")
            summary_frames.append(df_peak)
            print(f"✓ {ds_id} processed")
        except Exception as e:
            print(f"✗ {ds_id} failed → {e}")

    if not summary_frames:
        raise ValueError("All dataset IDs failed; nothing to plot.")

    df_summary = pd.concat(summary_frames, ignore_index=True)

    # ---------------- 2. filter summary -----------------------------------
    df_filtered = df_summary.query(
        "biomarker == @selected_biomarker and "
        "reference_event == @selected_reference_event and "
        "n_participant >= @selected_min_nparticipant"
    ).copy()

    if df_filtered.empty:
        raise ValueError(
            "No rows match the chosen biomarker/reference_event "
            "with ≥ selected_min_nparticipant."
        )

    df_filtered = df_filtered.sort_values(["specimen", "dataset_id"])

    # ---------------- 3. build box-plot stats -----------------------------
    bxp_stats, specimen_order = [], []
    for _, r in df_filtered.iterrows():
        bxp_stats.append(
            {
                "label": f"{r['dataset_id']} (n={r['n_participant']})",
                "whislo": r["shedding_peak_min"],
                "q1":     r["shedding_peak_q25"],
                "med":    r["shedding_peak_median"],
                "q3":     r["shedding_peak_q75"],
                "whishi": r["shedding_peak_max"],
            }
        )
        specimen_order.append(r["specimen"])

    unique_specimens = pd.unique(specimen_order)
    cmap = cm.get_cmap("tab10", len(unique_specimens))
    spec_colors = {spec: cmap(i) for i, spec in enumerate(unique_specimens)}

    # ---------------- 4. plotting -----------------------------------------
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, max(4, 0.4 * len(bxp_stats))))
    else:
        fig = ax.figure

    bp = ax.bxp(bxp_stats, vert=False, patch_artist=True, showfliers=False)

    for line in bp["medians"]:
        line.set_color("black")
        line.set_linewidth(1.5)

    for patch, spec in zip(bp["boxes"], specimen_order):
        patch.set_facecolor(spec_colors[spec])
        patch.set_edgecolor("black")
        patch.set_alpha(0.75)

    if x_axis_upper_limit is not None:
        ax.set_xlim(-5, x_axis_upper_limit)

    ax.set_xlabel("Shedding peak (time units)")
    ax.set_ylabel("Study")
    ax.set_title(
        f"{selected_biomarker} shedding peak by study "
        f"(ref = {selected_reference_event})"
    )
    ax.axvline(0, color="gray", linestyle="--", linewidth=1)

    handles = [
        plt.Line2D(
            [0], [0], marker="s", linestyle="",
            markerfacecolor=spec_colors[s], markeredgecolor="black",
            label=s, markersize=12
        )
        for s in unique_specimens
    ]
    ax.legend(
        handles, [h.get_label() for h in handles],
        title="Specimen",
        bbox_to_anchor=(1.02, 1), loc="upper left"
    )

    plt.tight_layout()
    return fig, ax