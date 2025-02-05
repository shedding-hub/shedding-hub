import os
import pandas as pd
import yaml
from shedding_hub import folded_str


hakki2022 = pd.read_csv("trajectories.csv")



# Define a dictionary containing patient information (ID, Sex, Age) from [Table 1](https://journals.asm.org/doi/10.1128/msphere.00827-20#tab1) in Vetter et al. (2020).
patient_info = {
    '1': {'PatientID': 1, 'Sex': 'F', 'Age': 35},
    '2': {'PatientID': 2, 'Sex': 'M', 'Age': 37},
    '3': {'PatientID': 3, 'Sex': 'M', 'Age': 29},
    '4': {'PatientID': 4, 'Sex': 'M', 'Age': 27},
    '5': {'PatientID': 5, 'Sex': 'F', 'Age': 23},
    '6': {'PatientID': 6, 'Sex': 'F', 'Age': 27},
    '7': {'PatientID': 7, 'Sex': 'F', 'Age': 25},
    '8': {'PatientID': 8, 'Sex': 'F', 'Age': 28},
    '9': {'PatientID': 9, 'Sex': 'M', 'Age': 47},
    '10': {'PatientID': 10, 'Sex': 'F', 'Age': 60},
    '11': {'PatientID': 11, 'Sex': 'F', 'Age': 32},
    '12': {'PatientID': 12, 'Sex': 'M', 'Age': 50},
    '13': {'PatientID': 13, 'Sex': 'M', 'Age': 46},
    '14': {'PatientID': 14, 'Sex': 'M', 'Age': 51},
    '15': {'PatientID': 15, 'Sex': 'F', 'Age': 45},
    '16': {'PatientID': 16, 'Sex': 'F', 'Age': 32},
    '17': {'PatientID': 17, 'Sex': 'F', 'Age': 23},
    '18': {'PatientID': 18, 'Sex': 'F', 'Age': 52},
    '19': {'PatientID': 19, 'Sex': 'M', 'Age': 26},
    '20': {'PatientID': 20, 'Sex': 'M', 'Age': 63},
    '21': {'PatientID': 21, 'Sex': 'F', 'Age': 49},
    '22': {'PatientID': 22, 'Sex': 'F', 'Age': 16},
    '23': {'PatientID': 23, 'Sex': 'F', 'Age': 40},
    '24': {'PatientID': 24, 'Sex': 'F', 'Age': 49},
    '25': {'PatientID': 25, 'Sex': 'F', 'Age': 48},
    '26': {'PatientID': 26, 'Sex': 'M', 'Age': 12},
    '27': {'PatientID': 27, 'Sex': 'M', 'Age': 15},
    '28': {'PatientID': 28, 'Sex': 'M', 'Age': 15},
    '29': {'PatientID': 29, 'Sex': 'F', 'Age': 33},
    '30': {'PatientID': 30, 'Sex': 'M', 'Age': 36},
    '31': {'PatientID': 31, 'Sex': 'M', 'Age': 13},
    '32': {'PatientID': 32, 'Sex': 'F', 'Age': 10},
    '33': {'PatientID': 33, 'Sex': 'F', 'Age': 26},
    '34': {'PatientID': 34, 'Sex': 'M', 'Age': 41},
    '35': {'PatientID': 35, 'Sex': 'F', 'Age': 42},
    '36': {'PatientID': 36, 'Sex': 'M', 'Age': 49},
    '37': {'PatientID': 37, 'Sex': 'F', 'Age': 55},
    '38': {'PatientID': 38, 'Sex': 'F', 'Age': 44},
    '39': {'PatientID': 39, 'Sex': 'F', 'Age': 41},
    '40': {'PatientID': 40, 'Sex': 'F', 'Age': 43},
    '41': {'PatientID': 41, 'Sex': 'M', 'Age': 53},
    '42': {'PatientID': 42, 'Sex': 'F', 'Age': 55},
    '43': {'PatientID': 43, 'Sex': 'M', 'Age': 36},
    '44': {'PatientID': 44, 'Sex': 'F', 'Age': 48},
    '45': {'PatientID': 45, 'Sex': 'M', 'Age': 43},
    '46': {'PatientID': 46, 'Sex': 'M', 'Age': 49},
    '47': {'PatientID': 47, 'Sex': 'F', 'Age': 36},
    '48': {'PatientID': 48, 'Sex': 'F', 'Age': 60},
    '49': {'PatientID': 49, 'Sex': 'F', 'Age': 46},
    '50': {'PatientID': 50, 'Sex': 'F', 'Age': 54},
    '51': {'PatientID': 51, 'Sex': 'M', 'Age': 39},
    '52': {'PatientID': 52, 'Sex': 'F', 'Age': 49},
    '53': {'PatientID': 53, 'Sex': 'M', 'Age': 41},
    '54': {'PatientID': 54, 'Sex': 'M', 'Age': 41},
    '55': {'PatientID': 55, 'Sex': 'M', 'Age': 15},
    '56': {'PatientID': 56, 'Sex': 'F', 'Age': 45},
    '57': {'PatientID': 57, 'Sex': 'F', 'Age': 45},
}


def map_patient_info(df):
    df = df.copy() # Create a copy of the DataFrame to avoid modifying the original data
    df['participant'] = df['participant'].astype(str)
     # Map 'PatientID', 'Sex', and 'Age' based on 'participant' column using the patient_info dictionary
    df.loc[:, 'PatientID'] = df['participant'].map(lambda x: patient_info.get(x, {}).get('PatientID'))
    df.loc[:, 'Sex'] = df['participant'].map(lambda x: patient_info.get(x, {}).get('Sex'))
    df.loc[:, 'Age'] = df['participant'].map(lambda x: patient_info.get(x, {}).get('Age'))

    return df



# Apply the mapping function to df_1 and save the updated DataFrame to a CSV file
hakki2022 = map_patient_info(hakki2022)
hakki2022['participant'] = hakki2022['participant'].astype(int)
hakki2022 = hakki2022.sort_values(by=['PatientID','day'])
columns_to_drop = ["pfu", "LFD", "copy_exist", "pfu_exist", "LFD_exist", "vaccinated", "WGS"]
hakki2022 = hakki2022.drop(columns=columns_to_drop)

# Group by participant and extract measurements
participants = []

for patient_id, group in hakki2022.groupby("participant"):
    # Initialize participant with a default structure
    participant = {"attributes": {}, "measurements": []}

    # Check if 'Age' or 'Sex' is not null and add attributes
    if pd.notnull(group["Age"].iloc[0]) and pd.notnull(group["Sex"].iloc[0]):
        participant["attributes"] = {
            "age": int(group["Age"].iloc[0]),
            "sex": "female" if group["Sex"].iloc[0] == "F" else "male",
        }

    # Process measurements
    for _, row in group.iterrows():
        # Skip rows where 'copy' is NaN
        if pd.isna(row["copy"]):
            continue
        if row["copy"] == 1:
            value = "negative"
        else:
            value = row["copy"]

        measurementN = {
            "analyte": "swab_SARSCoV2_N",
            "time": row["day"],
            "value": value,
        }
        participant["measurements"].append(measurementN)

    participants.append(participant)




hakki2022 = dict(
    title="Onset and window of SARS-CoV-2 infectiousness and temporal correlation with symptom onset: a prospective, longitudinal, community cohort study",
    doi="10.1016/S2213-2600(22)00226-0",
    description=folded_str(
        "Through a prospective, longitudinal, community cohort study that captures the critical growth phase and peak of viral replication, the goal is to characterize the window of SARS-CoV-2 infectiousness and its temporal relationship with symptom onset.\n"
    ),
    analytes=dict(
        swab_upper_res=dict(
            description=folded_str(
                "This analyte represents the detection and quantification of SARS-CoV-2 viral RNA from throat and nose swabs specimens collected from participants. The analysis focuses on measuring the viral load expressed in log10 copies/mL, with the primary reference event being the onset of symptoms.\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen="nose_and_throat_swab",
            biomarker="SARS-CoV-2",
            unit="gc/mL",
            reference_event="symptom onset",
        )
    ),
    participants=participants,
)

with open("hakki2022onset.yaml", "w") as outfile:
    outfile.write("# yaml-language-server:$schema=../.schema.yaml\n")
    yaml.dump(hakki2022, outfile, default_flow_style=False, sort_keys=False)
outfile.close()