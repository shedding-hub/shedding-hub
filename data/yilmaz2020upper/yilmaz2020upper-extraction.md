# Extraction for ylimaz2020 et al. (2020)

[ylimaz2020 et al. (2020)](https://academic.oup.com/jid/article/223/1/15/5918189?login=false#221942476) The authors reports longitudinal viral RNA loads from the nasopharynx/throat in patients with mild and severe/critical coronavirus disease 2019 (COVID-19). He also investigated whether the duration of symptoms correlated with the duration of viral RNA shedding. A total of 56 patients were included. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub). 

First, we `import` python modules needed:
```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```
```python
# Step 1: Read in data from Excel files
data = pd.read_excel("CombinedDataset.xlsx") 
data = data[data['StudyNum'] == 14] 

# Step 2: Replace values in the Series
data = data.replace({
    "M": "male", 
    "F": "female", 
    "Moderate": "moderate", 
    "Mild": "mild", 
    "Severe": "severe", 
    "14-1": "1", 
    "14-2": "2",
    "14-3": "3", "14-4": "4", "14-5": "5", "14-6": "6",
        "14-7": "7", "14-8": "8", "14-9": "9", "14-10": "10", "14-11": "11", "14-12": "12",
        "14-13": "13", "14-14": "14", "14-15": "15", "14-16": "16", "14-17": "17",
        "14-18": "18", "14-19": "19", "14-20": "20", "14-21": "21", "14-22": "22",
        "14-23": "23", "14-24": "24", "14-25": "25", "14-26": "26", "14-27": "27",
        "14-28": "28", "14-29": "29", "14-30": "30", "14-31": "31", "14-32": "32",
        "14-33": "33", "14-34": "34", "14-35": "35", "14-36": "36", "14-37": "37",
        "14-38": "38", "14-39": "39", "14-40": "40", "14-41": "41", "14-42": "42",
        "14-43": "43", "14-44": "44", "14-45": "45", "14-46": "46", "14-47": "47",
        "14-48": "48", "14-49": "49", "14-50": "50", "14-51": "51", "14-52": "52",
        "14-53": "53", "14-54": "54"
    
})


```
```python
df = pd.DataFrame(data) if not isinstance(data, pd.DataFrame) else data

participants = []

for patient_id, patient_data in df.groupby("PatientID"):
#patient_id in df["PatientID"].unique():
#     patient_data = df[df["PatientID"] == patient_id]
    participant = {
        "attributes": {
            # "day": int(patient_data["Day"].iloc[0]),
            "age": float(patient_data["Age"].iloc[0]),
            "sex": str(patient_data["Sex"].iloc[0]),
            # "estimated": int(patient_data["Estimated"].iloc[0]),
            # "sevmax": float(patient_data["SevMax"].iloc[0]),
        },
        # "sev1st": float(patient_data["Sev1st"].iloc[0]),
        #"death": int(patient_data["Died"].iloc[0]),
        "measurements": []
    }
    # print(participant)
    for _, row in patient_data.iterrows():
        measurement = {
            "analyte": "throatswab_SARSCoV2",
            "time": int(row["Day"]),
            "value": float(row["value"])
        }
        participant["measurements"].append(measurement)


    participants.append(participant)

    # Check column names
print(df.columns)

# Strip spaces from column names to avoid issues with whitespace
df.columns = df.columns.str.strip()

# Verify the column exists after cleaning
if 'PatientID' not in df.columns:
    print("Column 'PatientID' not found!")
else:
    participants = []

    for patient_id, patient_data in df.groupby("PatientID"):
        participant = {
            "attributes": {
                "age": float(patient_data["Age"].iloc[0]),
                "sex": str(patient_data["Sex"].iloc[0]),
            },
            "measurements": []
        }

        for _, row in patient_data.iterrows():
            measurement = {
                "analyte": "throatswab_SARSCoV2",
                "time": int(row["Day"]),
                "value": float(row["value"])
            }
            participant["measurements"].append(measurement)

        participants.append(participant)

    print(participants)
```

Finally, the data is formatted and output as a YAML file.
```python
output_data = {
    "title": "Upper Respiratory Tract Levels of Severe Acute Respiratory Syndrome Coronavirus 2 RNA and Duration of Viral RNA Shedding Do Not Differ Between Patients With Mild and Severe/Critical Coronavirus Disease 2019",
    "doi": "10.1093/infdis/jiaa632",
    "description": folded_str(
        "The authors reports longitudinal viral RNA loads from the nasopharynx/throat in patients with mild and severe/critical coronavirus disease 2019 (COVID-19). He also investigated whether the duration of symptoms correlated with the duration of viral RNA shedding. A total of 56 patients were included.\n"
    ),
    "analytes": {
        "throatswab_SARSCoV2": {
            "description": folded_str(
                "The author collected serial upper respiratory tract samples (1 nasopharyngeal swab and 1 throat swab put in a single collection tube with 1 mL of trans- port medium) for real-time PCR of SARS-CoV-2 RNA for all patients.\n"
            ),
            "specimen": "throat_swab",
            "biomarker": "SARS-CoV-2",
            "gene_target": "RdRp",
            "limit_of_quantification": "unknown",
            "limit_of_detection": "unknown",
            "unit": "gc/swab",
            "reference_event": "symptom onset"
        }
    },
    "participants": participants
}
with open("yilmaz2020.yaml","w") as outfile:
    yaml.dump(output_data, outfile, default_flow_style=False, allow_unicode=True, sort_keys=False)
```
