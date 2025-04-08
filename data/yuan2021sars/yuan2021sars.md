```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```


```python
df = pd.read_csv("pone.0247367.s001.csv")
df_case = df[df['case'].notna()].copy()
df_case = df_case[[
    "id", "age", "gender1", "ORF1ab", "N",
    "date_collection", "date_onset", "date_detection", "type" 
]].copy()

df_case.columns = [
    "patient_id", "age", "sex", "ctvalue_ORF1ab", "ctvalue_N",
    "date_collection", "date_onset", "date_detection", "symptom_status"  
]

df_case["ctvalue_ORF1ab"] = pd.to_numeric(df_case["ctvalue_ORF1ab"], errors="coerce")
df_case["ctvalue_N"] = pd.to_numeric(df_case["ctvalue_N"], errors="coerce")
df_case["date_collection"] = pd.to_datetime(df_case["date_collection"], errors="coerce")
df_case["date_onset"] = pd.to_datetime(df_case["date_onset"], errors="coerce")
df_case["date_detection"] = pd.to_datetime(df_case["date_detection"], errors="coerce")

df_case["sex"] = df_case["sex"].replace({"M": "male", "F": "female"})
df_case["symptom_status"] = df_case["symptom_status"].str.lower().str.strip()
df_case["symptom_status"] = df_case["symptom_status"].map({
    "symptomatic": "symptomatic",
    "asymptomatic": "asymptomatic"
}).fillna("unknown")  

def compute_time(row):
    if row["symptom_status"] == "symptomatic" and pd.notna(row["date_onset"]):
        return (row["date_collection"] - row["date_onset"]).days
    elif row["symptom_status"] == "asymptomatic" and pd.notna(row["date_detection"]):
        return (row["date_collection"] - row["date_detection"]).days
    else:
        return None  

df_case["time"] = df_case.apply(compute_time, axis=1)

positive_patients = df_case[
    (df_case["ctvalue_ORF1ab"] < 40) | (df_case["ctvalue_N"] < 40)
]["patient_id"].unique()
df_filtered = df_case[df_case["patient_id"].isin(positive_patients)].copy()

```


```python
participants = [] 

for patient_id, group in df_filtered.groupby("patient_id"):
    participant = {
        "attributes": {
            "id": patient_id,
            "age": float(group["age"].iloc[0]) if pd.notna(group["age"].iloc[0]) else "unknown",
            "sex": group["sex"].iloc[0] if pd.notna(group["sex"].iloc[0]) else "unknown",
            "symptom_status": group["symptom_status"].iloc[0] if pd.notna(group["symptom_status"].iloc[0]) else "unknown"
        },
        "measurements": []
    }

    for _, row in group.iterrows():
        time = int(row["time"]) if pd.notna(row["time"]) and row["time"] >= 0 else "unknown"

        if pd.notna(row["ctvalue_ORF1ab"]):
            participant["measurements"].append({
                "analyte": "stool_SARSCoV2_ORF1ab",
                "time": time,
                "value": float(row["ctvalue_ORF1ab"]) if row["ctvalue_ORF1ab"] < 40 else "negative"
            })

        if pd.notna(row["ctvalue_N"]):
            participant["measurements"].append({
                "analyte": "stool_SARSCoV2_N",
                "time": time,
                "value": float(row["ctvalue_N"]) if row["ctvalue_N"] < 40 else "negative"
            })

    if participant["measurements"]:
        participants.append(participant)

```


```python
output_data = {
    "title": "Fecal viral shedding in COVID-19 patients: Clinical significance, viral load dynamics and survival analysis",
    "doi": "10.1016/j.virusres.2020.198147",
    "description": folded_str(
        "This study investigates the fecal shedding of SARS-CoV-2 in COVID-19 patients, "
        "analyzing viral load dynamics, time course of shedding, and clinical significance. "
        "It reports detection of viral RNA in stool over extended time following symptom onset."
    ),
    "analytes": {
        "stool_SARSCoV2_ORF1ab": {
            "description": folded_str(
                "qPCR analysis of SARS-CoV-2 RNA in stool samples targeting ORF1ab gene. "
                "Ct less than 40 is considered positive in this analysis."
            ),
            "specimen": "stool",
            "biomarker": "SARS-CoV-2",
            "gene_target": "ORF1ab",
            "limit_of_quantification": "unknown",
            "limit_of_detection": 40,
            "unit": "cycle threshold",
            "reference_event": "symptom onset"
        },
        "stool_SARSCoV2_N": {
            "description": folded_str(
                "qPCR analysis of SARS-CoV-2 RNA in stool samples targeting N gene. "
                "Ct less than 40 is considered positive in this analysis."
            ),
            "specimen": "stool",
            "biomarker": "SARS-CoV-2",
            "gene_target": "N",
            "limit_of_quantification": "unknown",
            "limit_of_detection": 40,
            "unit": "cycle threshold",
            "reference_event": "symptom onset"
        }
    },
    "participants": participants 
}
with open("yuan2021sars.yaml","w") as outfile:
    outfile.write("# yaml-language-server:$schema=../.schema.yaml\n")
    yaml.dump(output_data, outfile, default_flow_style=False, sort_keys=False)
```

```python

```
