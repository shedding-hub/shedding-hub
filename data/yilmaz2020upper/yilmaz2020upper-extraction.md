# Extraction for ylimaz2020 et al. (2020)

[ylimaz2020 et al. (2020)](https://academic.oup.com/jid/article/223/1/15/5918189?login=false#221942476) The authors reports longitudinal viral RNA loads from the nasopharynx/throat in patients with mild and severe/critical coronavirus disease 2019 (COVID-19). He also investigated whether the duration of symptoms correlated with the duration of viral RNA shedding. A total of 56 patients were included. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub). 

First, we `import` python modules needed:
```python
import yaml
import pandas as pd
import numpy as np
class folded_str(str): pass
class literal_str(str): pass

def folded_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='>')
def literal_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(folded_str, folded_str_representer)
yaml.add_representer(literal_str, literal_str_representer)
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

# Step 3: Create the list of participants
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




```
Finally, the data is formatted and output as a YAML file.
```python
output_data = {
    "title": "Viral kinetics of SARS-CoV-2 in asymptomatic carriers and presymptomatic patients",
    "doi": "10.1016/j.ijid.2020.04.083",
    "description": folded_str(
        "The authors measured SARS-CoV-2 in longitudinal throat swab samples collected from 71 COVID-19 patients between February 4 and April 7, 2020. Abundances were quantified using real-time reverse transcription polymerase chain reaction (RT-PCR). Specimens were collected from all patients at least 2 days after hospitalization and physicians checked their symptoms and signs daily. Patients who had asymptomatic carrier and incubation period were analyzed.\n"
    ),
    "analytes": {
        "throatswab_SARSCoV2": {
            "description": folded_str(
                "SARS-CoV-2 RNA genome copy concentration in throat swab samples. Specimens were collected from all patients at least 2 days after hospitalization and physicians checked their symptoms and signs daily.\n"
            ),
            "specimen": "throat_swab",
            "biomarker": "SARS-CoV-2",
            "gene_target": "RdRp",
            "limit_of_quantification": "unknown",
            "limit_of_detection": 40,
            "unit": "value",
            "reference_event": "symptom onset"
        }
    },
    "participants": participants
}
with open("kimse2020.yaml","w") as outfile:
    yaml.dump(output_data, outfile, default_flow_style=False, allow_unicode=True, sort_keys=False)

```