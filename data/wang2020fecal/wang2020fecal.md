```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```


```python
df_a = pd.read_excel("A data.xlsx")
df_b = pd.read_excel("B data.xlsx")
df_a.columns = ["PatientID", "Day", "Ctvalue"]
df_b.columns = ["PatientID", "Day", "Ctvalue"]
df_a["Day"] = df_a["Day"].astype(int)
df_a["Ctvalue"] = pd.to_numeric(df_a["Ctvalue"], errors="coerce")
df_b["Day"] = df_b["Day"].astype(int)
df_b["Ctvalue"] = pd.to_numeric(df_b["Ctvalue"], errors="coerce")
```


```python
participants_ORF1ab = []    
for patient_id, patient_data in df_a.groupby("PatientID"):
        participant = {
            "attributes": {},  
            "measurements": []
        }
        for _, row in patient_data.iterrows():
            if pd.notna(row['Ctvalue']):
                measurement_ct = {
                    "analyte": "stool_SARSCoV2",
                    "time": row["Day"],
                    "value": float(row["Ctvalue"]) if float(row["Ctvalue"]) < 40 else "negative"
                }
                participant['measurements'].append(measurement_ct)

participants_ORF1ab.append(participant)
```


```python
participants_N = []
for patient_id, patient_data in df_b.groupby("PatientID"):
    participant = {
        "attributes": {},  
        "measurements": []
    }
    for _, row in patient_data.iterrows():
        measurement = {
            "analyte": "stool_SARSCoV2",
            "time": round(float(row["Day"])),
            "value": float(row["Ctvalue"]) if float(row["Ctvalue"]) < 40 else "negative"
        }
        participant["measurements"].append(measurement)

    participants_N.append(participant)
```


```python
participants = []
participants.extend(participants_ORF1ab)
participants.extend(participants_N)
```


```python
output_data = {
    "title": "Fecal viral shedding in COVID-19 patients: Clinical significance, viral load dynamics and survival analysis",
    "doi": "10.1016/j.virusres.2020.198147",
    "description": "This study investigates the fecal shedding of SARS-CoV-2 in COVID-19 patients, analyzing viral load dynamics, clinical significance, and survival analysis.",
    "analytes": {
        "stool_SARSCoV2_ORF1ab": {
            "description": "qPCR analysis of SARS-CoV-2 RNA in stool samples targeting ORF1ab gene. Ct = 35 is the cut-off for a positive result, and Ct = 40 is a negative sample.",
            "specimen": "stool",
            "biomarker": "SARS-CoV-2",
            "gene_target": "ORF1ab",
            "limit_of_detection": 40,
            "unit": "cycle threshold",
            "reference_event": "symptom onset"
        },
        "stool_SARSCoV2_N": {
            "description": "qPCR analysis of SARS-CoV-2 RNA in stool samples targeting N gene. Ct = 35 is the cut-off for a positive result, and Ct = 40 is a negative sample.",
            "specimen": "stool",
            "biomarker": "SARS-CoV-2",
            "gene_target": "N",
            "limit_of_detection": 40,
            "unit": "cycle threshold",
            "reference_event": "hospital admission"
        }
    },  
    "participants": participants
}
with open("wang2020fecal.yaml","w") as outfile:
    outfile.write("# yaml-language-server:$schema=../.schema.yaml\n")
    yaml.dump(output_data, outfile, default_flow_style=False, sort_keys=False)
```
