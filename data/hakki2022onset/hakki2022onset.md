# Extraction for hakki et al. (2022)

[Hakki et al. (2022)](https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(22)00226-0/fulltext) analyzed the infectious window of SARS-CoV-2 using longitudinal upper respiratory tract (URT) sampling. A cohort of 57 patients underwent daily URT sampling for up to 20 days. SARS-CoV-2 viral RNA levels were quantified from URT swabs, with data extracted from the **trajectories.csv** dataset, which is stored at [GitHub](https://github.com/HPRURespMed/SARS-CoV-2-viral-shedding-dynamics/tree/main). Demographic variables such as age and sex were not included in the dataset.
First, we import `python` modules needed:
```python
import pandas as pd
import yaml
from shedding_hub import folded_str
```
We clean data and add the demographic information in datasets:
```python
hakki2022 = pd.read_csv(
    "trajectories.csv"
)  # The data was obtained from https://github.com/HPRURespMed/SARS-CoV-2-viral-shedding-dynamics/blob/main/Data/trajectories.csv


# Define a dictionary containing patient information (ID, Sex, Age) from [Figure 3](https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(22)00226-0/fulltext) in hakki et al. (2022).
patient_info = {
    "1": {"PatientID": 1, "Sex": "F", "Age": 35},
    "2": {"PatientID": 2, "Sex": "M", "Age": 37},
    "3": {"PatientID": 3, "Sex": "M", "Age": 29},
    "4": {"PatientID": 4, "Sex": "M", "Age": 27},
    "5": {"PatientID": 5, "Sex": "F", "Age": 23},
    "6": {"PatientID": 6, "Sex": "F", "Age": 27},
    "7": {"PatientID": 7, "Sex": "F", "Age": 25},
    "8": {"PatientID": 8, "Sex": "F", "Age": 28},
    "9": {"PatientID": 9, "Sex": "M", "Age": 47},
    "10": {"PatientID": 10, "Sex": "F", "Age": 60},
    "11": {"PatientID": 11, "Sex": "F", "Age": 32},
    "12": {"PatientID": 12, "Sex": "M", "Age": 50},
    "13": {"PatientID": 13, "Sex": "M", "Age": 46},
    "14": {"PatientID": 14, "Sex": "M", "Age": 51},
    "15": {"PatientID": 15, "Sex": "F", "Age": 45},
    "16": {"PatientID": 16, "Sex": "F", "Age": 32},
    "17": {"PatientID": 17, "Sex": "F", "Age": 23},
    "18": {"PatientID": 18, "Sex": "F", "Age": 52},
    "19": {"PatientID": 19, "Sex": "M", "Age": 26},
    "20": {"PatientID": 20, "Sex": "M", "Age": 63},
    "21": {"PatientID": 21, "Sex": "F", "Age": 49},
    "22": {"PatientID": 22, "Sex": "F", "Age": 16},
    "23": {"PatientID": 23, "Sex": "F", "Age": 40},
    "24": {"PatientID": 24, "Sex": "F", "Age": 49},
    "25": {"PatientID": 25, "Sex": "F", "Age": 48},
    "26": {"PatientID": 26, "Sex": "M", "Age": 12},
    "27": {"PatientID": 27, "Sex": "M", "Age": 15},
    "28": {"PatientID": 28, "Sex": "M", "Age": 15},
    "29": {"PatientID": 29, "Sex": "F", "Age": 33},
    "30": {"PatientID": 30, "Sex": "M", "Age": 36},
    "31": {"PatientID": 31, "Sex": "M", "Age": 13},
    "32": {"PatientID": 32, "Sex": "F", "Age": 10},
    "33": {"PatientID": 33, "Sex": "F", "Age": 26},
    "34": {"PatientID": 34, "Sex": "M", "Age": 41},
    "35": {"PatientID": 35, "Sex": "F", "Age": 42},
    "36": {"PatientID": 36, "Sex": "M", "Age": 49},
    "37": {"PatientID": 37, "Sex": "F", "Age": 55},
    "38": {"PatientID": 38, "Sex": "F", "Age": 44},
    "39": {"PatientID": 39, "Sex": "F", "Age": 41},
    "40": {"PatientID": 40, "Sex": "F", "Age": 43},
    "41": {"PatientID": 41, "Sex": "M", "Age": 53},
    "42": {"PatientID": 42, "Sex": "F", "Age": 55},
    "43": {"PatientID": 43, "Sex": "M", "Age": 36},
    "44": {"PatientID": 44, "Sex": "F", "Age": 48},
    "45": {"PatientID": 45, "Sex": "M", "Age": 43},
    "46": {"PatientID": 46, "Sex": "M", "Age": 49},
    "47": {"PatientID": 47, "Sex": "F", "Age": 36},
    "48": {"PatientID": 48, "Sex": "F", "Age": 60},
    "49": {"PatientID": 49, "Sex": "F", "Age": 46},
    "50": {"PatientID": 50, "Sex": "F", "Age": 54},
    "51": {"PatientID": 51, "Sex": "M", "Age": 39},
    "52": {"PatientID": 52, "Sex": "F", "Age": 49},
    "53": {"PatientID": 53, "Sex": "M", "Age": 41},
    "54": {"PatientID": 54, "Sex": "M", "Age": 41},
    "55": {"PatientID": 55, "Sex": "M", "Age": 15},
    "56": {"PatientID": 56, "Sex": "F", "Age": 45},
    "57": {"PatientID": 57, "Sex": "F", "Age": 45},
}


def map_patient_info(df):
    df = (
        df.copy()
    )  # Create a copy of the DataFrame to avoid modifying the original data
    df["participant"] = df["participant"].astype(str)
    # Map 'PatientID', 'Sex', and 'Age' based on 'participant' column using the patient_info dictionary
    df.loc[:, "PatientID"] = df["participant"].map(
        lambda x: patient_info.get(x, {}).get("PatientID")
    )
    df.loc[:, "Sex"] = df["participant"].map(
        lambda x: patient_info.get(x, {}).get("Sex")
    )
    df.loc[:, "Age"] = df["participant"].map(
        lambda x: patient_info.get(x, {}).get("Age")
    )

    return df


# Apply the mapping function to df_1 and save the updated DataFrame to a CSV file
hakki2022 = map_patient_info(hakki2022)
hakki2022["participant"] = hakki2022["participant"].astype(int)
hakki2022 = hakki2022.sort_values(by=["PatientID", "day"])
columns_to_drop = ["LFD", "copy_exist", "pfu_exist", "LFD_exist"]
hakki2022 = hakki2022.drop(columns=columns_to_drop)

patient_ids = [
    1,
    5,
    6,
    7,
    8,
    12,
    16,
    20,
    21,
    22,
    23,
    25,
    30,
    32,
    45,
    47,
    50,
    52,
    56,
    57,
]

# Ensure 'PatientID' exists
if "PatientID" in hakki2022.columns:
    # Patients without symptoms
    nosymptomdataset = hakki2022[hakki2022["PatientID"].isin(patient_ids)]

    # Patients with symptoms (all others not in patient_ids)
    symptomdataset = hakki2022[~hakki2022["PatientID"].isin(patient_ids)]

patient_ids_symptom = [
    2,
    3,
    4,
    9,
    10,
    11,
    13,
    14,
    15,
    17,
    18,
    19,
    24,
    26,
    27,
    28,
    29,
    31,
    33,
    34,
    35,
    36,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    46,
    48,
    49,
    50,
    51,
    53,
    54,
    55,
]
day_values = [
    -1,
    -1,
    1,
    -4,
    -1,
    -1,
    -1,
    -3,
    -1,
    0,
    -2,
    -1,
    -1,
    0,
    -3,
    1,
    -1,
    -4,
    -2,
    -3,
    0,
    -1,
    -2,
    -1,
    -2,
    0,
    -2,
    -1,
    -3,
    -3,
    -4,
    -3,
    -2,
    2,
    -4,
    -2,
    -1,
    -2,
]
# Modify the 'day' column for the symptom group cycling through day_values
if "PatientID" in symptomdataset.columns and "day" in symptomdataset.columns:

    for i, patient in enumerate(patient_ids_symptom):
        if i < len(day_values):
            symptomdataset.loc[
                symptomdataset["PatientID"] == patient, "day"
            ] += day_values[i]

# Group by participant and extract measurements
participants_nasooroph = []

for patient_id, group in nosymptomdataset.groupby("participant"):
    # Initialize participant with a default structure
    participant = {"attributes": {}, "measurements": []}

    # Check if 'Age' or 'Sex' is not null and add attributes
    if pd.notnull(group["Age"].iloc[0]) and pd.notnull(group["Sex"].iloc[0]):
        participant["attributes"] = {
            "age": int(group["Age"].iloc[0]),
            "sex": "female" if group["Sex"].iloc[0] == "F" else "male",
            # "vaccinated": False if group["vaccinated"].iloc[0] == "FALSE" else True,
            "vaccinated": (
                False
                if str(group["vaccinated"].iloc[0]).strip().lower() == "false"
                else True
            ),
            "variant_wave": (
                "Pre-Alpha"
                if group["WGS"].iloc[0] == "Pre-Alpha"
                else (
                    "Alpha"
                    if group["WGS"].iloc[0] == "Alpha"
                    else "Delta" if group["WGS"].iloc[0] == "Delta" else "Unknown"
                )
            ),
        }

    measurements = []
    for _, row in group.iterrows():
        if pd.isna(row["copy"]):
            continue

        value_naor = "negative" if row["copy"] == 1 else row["copy"]

        measurement_entry1 = {
            "analyte": "asymptomatic_PCR",
            "time": row["day"],
            "value": value_naor,
        }

        if pd.isna(row["pfu"]):
            continue

        value_pfu = "negative" if row["pfu"] == 1 else row["pfu"]

        measurement_entry2 = {
            "analyte": "asymptomatic_cultivable",
            "time": row["day"],
            "value": value_pfu,
        }

        measurements.append(measurement_entry1)
        measurements.append(measurement_entry2)

    participant["measurements"].extend(measurements)

    participants_nasooroph.append(participant)


participants_cultivable = []

for patient_id, group in symptomdataset.groupby("participant"):
    # Initialize participant with a default structure
    participant = {"attributes": {}, "measurements": []}

    # Check if 'Age' or 'Sex' is not null and add attributes
    if pd.notnull(group["Age"].iloc[0]) and pd.notnull(group["Sex"].iloc[0]):
        participant["attributes"] = {
            "age": int(group["Age"].iloc[0]),
            "sex": "female" if group["Sex"].iloc[0] == "F" else "male",
            # "vaccinated": False if group["vaccinated"].iloc[0] == "FALSE" else True,
            "vaccinated": (
                False
                if str(group["vaccinated"].iloc[0]).strip().lower() == "false"
                else True
            ),
            "variant_wave": (
                "Pre-Alpha"
                if group["WGS"].iloc[0] == "Pre-Alpha"
                else (
                    "Alpha"
                    if group["WGS"].iloc[0] == "Alpha"
                    else "Delta" if group["WGS"].iloc[0] == "Delta" else "Unknown"
                )
            ),
        }

    measurements = []
    for _, row in group.iterrows():
        if pd.isna(row["copy"]):
            continue

        value_naor = "negative" if row["copy"] == 1 else row["copy"]

        measurement_entry1 = {
            "analyte": "symptomatic_PCR",
            "time": row["day"],
            "value": value_naor,
        }

        if pd.isna(row["pfu"]):
            continue

        value_pfu = "negative" if row["pfu"] == 1 else row["pfu"]

        measurement_entry2 = {
            "analyte": "symptomatic_cultivable",
            "time": row["day"],
            "value": value_pfu,
        }

        measurements.append(measurement_entry1)
        measurements.append(measurement_entry2)

    participant["measurements"].extend(measurements)

    participants_cultivable.append(participant)


participants = []
participants.extend(participants_nasooroph)
participants.extend(participants_cultivable)

```

Finally, the data is formatted and output as a YAML file.
```python
hakki2022 = dict(
    title="Onset and window of SARS-CoV-2 infectiousness and temporal correlation with symptom onset: a prospective, longitudinal, community cohort study",
    doi="10.1016/S2213-2600(22)00226-0",
    description=folded_str(
        "Through a prospective, longitudinal, community cohort study that captures the critical growth phase and peak of viral replication, the goal is to characterize the window of SARS-CoV-2 infectiousness and its temporal relationship with symptom onset.\n"
    ),
    analytes=dict(
        asymptomatic_PCR=dict(
            description=folded_str(
                "This analyte represents the detection and quantification of SARS-CoV-2 viral RNA from throat and nose swabs specimens collected from asymptomatic participants. The analysis focuses on measuring the viral load expressed in log10 copies/mL, with the primary reference event being the enrollment.\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
            biomarker="SARS-CoV-2",
            unit="gc/mL",  # Nose and throat swabs were placed in 3 mL viral transport medium (VTM) of two brands (Copan Diagnostics, Murrieta, CA, USA; or MANTACC, Guangdong, China).
            reference_event="enrollment",
        ),
        asymptomatic_cultivable=dict(
            description=folded_str(
                "The infectious virus analyte consisted of culturable SARS-CoV-2 derived from throat and nasal swabs of asymptomatic participants, measured using a plaque assay in vitro. The unit of measurement was PFU/mL, with enrollment as the primary reference event.\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
            biomarker="SARS-CoV-2",
            unit="pfu/mL",  # Nose and throat swabs were placed in 3 mL viral transport medium (VTM) of two brands (Copan Diagnostics, Murrieta, CA, USA; or MANTACC, Guangdong, China).
            reference_event="enrollment",
        ),
        symptomatic_PCR=dict(
            description=folded_str(
                "This analyte represents the detection and quantification of SARS-CoV-2 viral RNA from throat and nose swabs specimens collected from symptomatic participants. The analysis focuses on measuring the viral load expressed in log10 copies/mL, with the primary reference event being the symptom onset.\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
            biomarker="SARS-CoV-2",
            unit="gc/mL",
            # Nose and throat swabs were placed in 3 mL viral transport medium (VTM) of two brands (Copan Diagnostics, Murrieta, CA, USA; or MANTACC, Guangdong, China).
            reference_event="symptom onset",
        ),
        symptomatic_cultivable=dict(
            description=folded_str(
                "The infectious virus analyte consisted of culturable SARS-CoV-2 derived from throat and nasal swabs of symptomatic participants, measured using a plaque assay in vitro. The unit of measurement was PFU/mL, with symptom onset as the primary reference event.\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
            biomarker="SARS-CoV-2",
            unit="pfu/mL",
            # Nose and throat swabs were placed in 3 mL viral transport medium (VTM) of two brands (Copan Diagnostics, Murrieta, CA, USA; or MANTACC, Guangdong, China).
            reference_event="symptom onset",
        ),
    ),
    participants=participants,
)

with open("hakki2022onset.yaml", "w") as outfile:
    outfile.write("# yaml-language-server:$schema=../.schema.yaml\n")
    yaml.dump(hakki2022, outfile, default_flow_style=False, sort_keys=False)
```
