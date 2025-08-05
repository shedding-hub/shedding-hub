import shedding_hub as sh
import pandas as pd
from typing import Optional, Dict, Any
import matplotlib.pyplot as plt
import glob
from pathlib import Path

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
    df_shedding_peak = pd.DataFrame(columns=["study_ID","participant_ID",#"age","sex",
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
    


def plot_individual_shedding(
        dataset = 