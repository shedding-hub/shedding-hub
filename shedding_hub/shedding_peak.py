import shedding_hub as sh
import pandas as pd
from typing import Optional
import matplotlib.pyplot as plt
import glob
from pathlib import Path

def calc_shedding_peak(
    df_ID: str,
    *,
    plotting: bool = False,
    ref: Optional[str] = None,
    pr: Optional[int] = None,
    local: Optional[str] = None,
) -> pd.DataFrame:
    """
    Calcualte summary statistics for the shedding peak using a loaded dataset by the `load_dataset` function.

    Args:
        df_ID: Dataset identifier, e.g., :code:`woelfel2020virological`.
        plotting: Create a plot for individual level of shedding peak by specimen type.
        ref: Git reference to load. Defaults to the most recent data on the :code:`main`
            branch of https://github.com/shedding-hub/shedding-hub and is automatically
            fetched if a :code:`pr` number is specified.
        pr: Pull request to fetch data from.
        local: Local directory to load data from.

    Returns:
        Summary table of shedding peak (min, max, mean) by biomarker and specimen.
    """
    # load the dataset;
    df = sh.load_dataset(dataset = df_ID, ref = ref, pr = pr, local = local)
    # initialize a pandas dataframe
    df_shedding_peak = pd.DataFrame(columns=["study_ID","participant_ID",#"age","sex",
                                                 "analyte","n_sample","first_sample","last_sample","shedding_peak"])
    df_analyte = pd.DataFrame(columns=["analyte","specimen","biomarker","reference_event"])

    # extract analyte data from the standardized shedding data loaded
    for key in df["analytes"]:
        row_new = {"analyte":key,
                   "specimen":df["analytes"][key]["specimen"],
                   "biomarker":df["analytes"][key]["biomarker"],
                   "reference_event":df["analytes"][key]["reference_event"]}
        df_analyte.loc[len(df_analyte)] = row_new
    
    # extract participant and measurement data from the standardized shedding data loaded
    participant_counter = 0
    for item in df["participants"]:
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
                "study_ID": df_ID,
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
        ["study_ID", "biomarker", "specimen", "reference_event"]
    ).agg(
        shedding_peak_min=("shedding_peak", "min"),
        shedding_peak_max=("shedding_peak", "max"),
        shedding_peak_mean=("shedding_peak", "mean"),
        n_sample=("n_sample", "sum"),
        n_participant = ("participant_ID", "nunique")
    ).reset_index()
    
    return(df_return)