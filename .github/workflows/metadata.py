import argparse
import pandas as pd
import requests
import yaml
import shedding_hub as sh
import glob
from pathlib import Path


def create_meta(
    schema:dict
) -> pd.DataFrame:
    '''
    Create data structure of metadata using schema file.

    Args:
        schema: Schema file loaded from GitHub repository.

    Returns:
        Empty Metadata dataframe with column names.
    '''
    # Initialize column names as a list with ID column
    column_names = ["ID"]
    # Add required fields of analyte into it;
    column_names = column_names + schema["$defs"]["analyte_specification"]["required"]
    # Add summary statistics into it;
    column_names = column_names + ["n_samples", "n_unique_participants", "n_negative", "n_positive", "n_quantified"]
    return pd.DataFrame(columns = column_names)


def safe_compare(a, b, operator):
    def compare(x, y):
        # Allow int and float to compare with each other
        if not (
            isinstance(x, (int, float)) and isinstance(y, (int, float))
        ) and type(x) != type(y):
            return False

        if operator == "==":
            return x == y
        elif operator == "<":
            return x < y
        elif operator == ">":
            return x > y
        elif operator == "<=":
            return x <= y
        elif operator == ">=":
            return x >= y
        else:
            raise ValueError(f"Unsupported operator: {operator}")

    if isinstance(a, list):
        return [compare(item, b) for item in a]
    else:
        return compare(a, b)
    

def append_meta(
    metadata: pd.DataFrame,
    dataset: str
) -> pd.DataFrame:
    '''
    Extract data from yaml file and append it to the metadata.

    Args:
        metadata: Metadata dataframe with column names.
        dataset: Dataset identifier, e.g., :code:`woelfel2020virological`.

    Returns:
        Metadata dataframe with extracted data appended.
    '''
    data = sh.load_dataset(dataset)

    #create measurements dataframe for the new dataset
    measurements = pd.DataFrame(columns=["participant_ID","analyte","time","value"])
    counter = 1
    for participant in data["participants"]:
        new_data = pd.DataFrame.from_dict(participant["measurements"])
        new_data["participant_ID"] = counter
        if len(measurements) == 0:
            measurements = new_data[["participant_ID","analyte","time","value"]]
        else:
            measurements = pd.concat([measurements, new_data], axis=0, join='inner', ignore_index=True)
        counter += 1

    for analyte in data["analytes"]:
        measurements_filtered = measurements.loc[measurements["analyte"] == analyte]
        #calculate summary statistics
        n_samples = len(measurements_filtered)
        n_unique_participants = len(set(measurements_filtered["participant_ID"]))
        n_negative = len(measurements_filtered.loc[
                         [(value == "negative") | 
                          ((data["analytes"][analyte]["unit"] != "cycle threshold") & 
                           (safe_compare(value, data["analytes"][analyte]["limit_of_detection"], "<"))) | 
                          ((data["analytes"][analyte]["unit"] == "cycle threshold") & 
                           (safe_compare(value, data["analytes"][analyte]["limit_of_detection"], ">")))
                         for value in measurements_filtered["value"]]
                        ])
        n_positive = len(measurements_filtered.loc[
                         [(value == "positive") | 
                          ((data["analytes"][analyte]["unit"] != "cycle threshold") &  
                           (safe_compare(value, data["analytes"][analyte]["limit_of_quantification"], "<"))) | 
                          ((data["analytes"][analyte]["unit"] == "cycle threshold") & 
                           (safe_compare(value, data["analytes"][analyte]["limit_of_detection"], "<") | 
                            ((data["analytes"][analyte]["limit_of_detection"] == "unknown") & (safe_compare(value, 0, ">")))))
                         for value in measurements_filtered["value"]]
                        ])
        n_quantified = len(measurements_filtered.loc[
                           [((data["analytes"][analyte]["unit"] != "cycle threshold") & 
                             (data["analytes"][analyte]["limit_of_quantification"] == "unknown") & 
                             (safe_compare(value, 0, ">"))) | 
                            ((data["analytes"][analyte]["unit"] != "cycle threshold") & 
                             (data["analytes"][analyte]["limit_of_quantification"] != "unknown") & 
                             (safe_compare(value, data["analytes"][analyte]["limit_of_quantification"], ">")))  
                         for value in measurements_filtered["value"]]  
                        ])
        new_line = [dataset] + [data["analytes"][analyte][key] for key in metadata.columns if key in list(data["analytes"][analyte].keys())] + [n_samples, n_unique_participants, n_negative, n_positive, n_quantified]
        metadata.loc[len(metadata)] = new_line
    return metadata


def __main__() -> None:
    url = "https://raw.githubusercontent.com/shedding-hub/shedding-hub/refs/heads/main/data/.schema.yaml"
    response = requests.get(url)
    schema = yaml.safe_load(response.text)
    metadata = create_meta(schema)

    # Append all the data into the metadata
    for filename in [Path(file).stem for file in Path("data").glob("*/*.yaml")]:
        print(f"Load the data: {filename}")
        metadata = append_meta(metadata, filename)

    # Remove "\n" from description
    metadata["description"] = [item.replace("\n", "") for item in metadata["description"]]

    # Convert to dict, then dump to YAML
    with open('data/metadata.yaml', 'w') as f:
        yaml.safe_dump(metadata.to_dict(orient="records"), f)
        

if __name__ == "__main__":
    __main__()