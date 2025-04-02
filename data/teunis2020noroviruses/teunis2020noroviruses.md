# Extraction for teunis et al.(2020)
[teunis et al. (2020)] analyzed the variation in infectivity and pathogenicity of norovirus. The dosed date was used as the reference for the 2 patients. Norovirus PCR and RNA results were extracted from stool and emesis from the 1-SMA database--012104.xls dataset.
First, we import `python` modules needed:
```python
import pandas as pd
import yaml
from shedding_hub import folded_str
```

We clean data and set new time value:
```python
# Load dataset
df = pd.read_excel("1-SMA Database--012104.xls", sheet_name=["Stool", "Emesis"])
stool_df = df["Stool"]
subset_stool = stool_df[stool_df.iloc[:, 0].astype(str).str.startswith("HV")]
# Define the desired column names for Stool
column_names_stool = [
     "Sample ID", "date", "time", "color", "cons", "amt", "day",
     "PCR", "RNA", "Qiagen", "", "", "", "", "PDU", "samples given to", "notes"
 ]
# Set the column names in your filtered HV stool data
subset_stool.columns = column_names_stool
columns_to_drop = ["color", "cons", "day", "Qiagen", "", "", "", "", "samples given to", "notes"]
subset_stool = subset_stool.drop(columns=columns_to_drop)


subset_stool["date"] = pd.to_datetime(subset_stool["date"], errors='coerce')
subset_stool["time"] = pd.to_datetime(subset_stool["time"], format="%H:%M:%S", errors='coerce')

# new time value for combine day and hour
subset_stool["datetime"] = pd.to_datetime(
    subset_stool["date"].dt.date.astype(str) + " " + subset_stool["time"].dt.time.astype(str)
)




# Extract the reference sample identifier (e.g., HV-AD, HV-TK)
subset_stool['reference'] = subset_stool['Sample ID'].str.extract(r'(^[A-Z]+-[A-Z]+)')

# Calculate the first time for each reference
# Get the first time of the reference sample (e.g., HV-AD, HV-TK)
first_times = subset_stool[subset_stool['Sample ID'].str.match(r'^[A-Z]+-[A-Z]+$')][['reference', 'datetime']]
first_times.columns = ['reference', 'first_time']

# Merge back to calculate the hour difference from the first reference time
subset_stool = pd.merge(subset_stool, first_times, on='reference')
subset_stool['hour_diff_from_first'] = (subset_stool['datetime'] - subset_stool['first_time']).dt.total_seconds() / 3600

# Replace the 'time' column with the hour difference
subset_stool['time'] = subset_stool['hour_diff_from_first']
subset_stool = subset_stool[subset_stool['Sample ID'].str.contains(r'-\d+-\d+$')]




emesis_df = df["Emesis"]
subset_emesis = emesis_df[emesis_df .iloc[:, 0].astype(str).str.startswith("HV")]
# Define the desired column names for Emesis
column_names_emesis = [
     "Sample ID", "date", "time", "amt", "day", "date p/u",
     "PCR", "RNA"
 ]
# Set the column names in your filtered HV Emesis data
subset_emesis.columns = column_names_emesis
columns_to_drop = [ "day", "date p/u"]
subset_emesis = subset_emesis.drop(columns=columns_to_drop)

subset_emesis["date"] = pd.to_datetime(subset_emesis["date"], errors='coerce')
subset_emesis["time"] = pd.to_datetime(subset_emesis["time"], format="%H:%M:%S", errors='coerce')

# new time value for combine day and hour
subset_emesis["datetime"] = pd.to_datetime(
    subset_emesis["date"].dt.date.astype(str) + " " + subset_emesis["time"].dt.time.astype(str)
)

# Extract the reference sample identifier (e.g., HV-AD, HV-TK)
subset_emesis['reference'] = subset_emesis['Sample ID'].str.extract(r'(^[A-Z]+-[A-Z]+)')

# Get the first time of the reference sample (e.g., HV-AD, HV-TK)
first_times = subset_emesis[subset_emesis['Sample ID'].str.match(r'^[A-Z]+-[A-Z]+$')][['reference', 'datetime']]
first_times.columns = ['reference', 'first_time']

# Merge back to calculate the hour difference from the first reference time
subset_emesis = pd.merge(subset_emesis, first_times, on='reference', how='left')
subset_emesis['hour_diff_from_first_e'] = (subset_emesis['datetime'] - subset_emesis['first_time']).dt.total_seconds() / 3600

# Replace the 'time' column with the hour difference
subset_emesis['time'] = subset_emesis['hour_diff_from_first_e']




participants_stool = []

# create a new column, exact the first two part of the patients'ID
subset_stool["Patient_ID"] = subset_stool['Sample ID'].str.split('-').str[:2].str.join('-')


df = subset_stool.groupby("Patient_ID")


for patient_id, group in df:
    participant = {"attributes": {}, "measurements": []}



    measurements = []
    for _, row in group.iterrows():
        if pd.isna(row["PCR"]):
            continue
        if row["PCR"] == '+':
            value_stoolPCR = "positive"
        else:
            value_stoolPCR = "negative"

        measurement_entry1 = {
            "analyte": "stool_PCR",
            "time": row["hour_diff_from_first"],
            "value": value_stoolPCR,
            "amount": row["amt"],
        }

        if pd.isna(row["RNA"]):
            continue
        if row["RNA"] == '+':
            value_stoolRNA = "positive"
        else:
            value_stoolRNA = "negative"

        measurement_entry2 = {
            "analyte": "stool_RNA",
            "time": row["hour_diff_from_first"],
            "value": value_stoolRNA,
            "amount": row["amt"],
        }

        measurements.append(measurement_entry1)
        measurements.append(measurement_entry2)

    participant["measurements"].extend(measurements)

    participants_stool.append(participant)



participants_emesis = []

# create a new column, exact the first two part of the patients'ID
subset_emesis["Patient_ID"] = subset_emesis["Sample ID"].str.split('-').str[:2].str.join('-')
df = subset_emesis.groupby("Patient_ID")
for patient_id, group in df:
    participant = {"attributes": {}, "measurements": []}
    measurements = []
    for _, row in group.iterrows():
        if pd.isna(row["PCR"]):
            continue
        if row["PCR"] == '+':
            value_emesisPCR = "positive"
        else:
            value_emesisPCR = "negative"
        measurement_entry1 = {
            "analyte": "emesis_PCR",
            "time": row["hour_diff_from_first_e"],
            "value": value_emesisPCR,
            "amount": row["amt"],
        }

        if pd.isna(row["RNA"]):
            continue
        if row["PCR"] == '+':
            value_emesisRNA = "positive"
        else:
            value_emesisRNA = "negative"

        measurement_entry2 = {
            "analyte": "emesis_RNA",
            "time": row["hour_diff_from_first_e"],
            "value": value_emesisRNA,
            "amount": row["amt"],
        }

        measurements.append(measurement_entry1)
        measurements.append(measurement_entry2)

    participant["measurements"].extend(measurements)

    participants_emesis.append(participant)

participants = []
participants.extend(participants_stool)
participants.extend(participants_emesis)
```
Finally, the data is formatted and output as a YAML file.

```python
teunis2020 = dict(
    title="Noroviruses are highly infectious but there is strong variation in host susceptibility and virus pathogenicity",
    doi="10.1016/j.epidem.2020.100401",
    description=folded_str(
        "By analyzing the variation in infectivity and pathogenicity of norovirus, the study compares these differences and describes how virus infectivity and host susceptibility vary. The results indicate that secretor-positive individuals experience significantly higher infection and illness rates compared to secretor-negative individuals, who exhibit protective effects.\n"
    ),
    analytes=dict(
        stool_PCR=dict(
            description=folded_str(
                "This analyte represents the PCR value in the stool\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen="stool",
            biomarker="norovirus",
            unit="NA",
            reference_event="dosed",
            sample_amount_unit="gram"

        ),
        stool_RNA=dict(
            description=folded_str(
                "The analyte represents the RNA value in the stool\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen="stool",
            biomarker="norovirus",
            unit="NA",
            reference_event="dosed",
            sample_amount_unit="gram",
        ),
        emesis_PCR=dict(
            description=folded_str(
                "This analyte represents the PCR in emesis\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen="emesis",
            biomarker="norovirus",
            unit="NA",
            reference_event="dosed",
            sample_amount_unit="gram",
        ),
        emesis_RNA=dict(
            description=folded_str(
                "The analyte represents the RNA in emesis\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection="unknown",
            specimen="emesis",
            biomarker="norovirus",
            unit="NA",
            reference_event="dosed",
            sample_amount_unit="gram",
        ),
    ),
    participants=participants,
)

with open("teunis2020noroviruses.yaml", "w") as outfile:
    outfile.write("# yaml-language-server:$schema=../.schema.yaml\n")
    yaml.dump(teunis2020, outfile, default_flow_style=False, sort_keys=False)
```
