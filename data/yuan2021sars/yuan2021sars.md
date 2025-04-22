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
    "date_collection", "date_onset", "date_detection", "type", "specimen"
]].copy()

df_case.columns = [
    "patient_id", "age", "sex", "ctvalue_ORF1ab", "ctvalue_N",
    "date_collection", "date_onset", "date_detection", "symptom_status", "specimen"
]

df_case["ctvalue_ORF1ab"] = pd.to_numeric(df_case["ctvalue_ORF1ab"], errors="coerce")
df_case["ctvalue_N"] = pd.to_numeric(df_case["ctvalue_N"], errors="coerce")
df_case["date_collection"] = pd.to_datetime(df_case["date_collection"], errors="coerce")
df_case["date_onset"] = pd.to_datetime(df_case["date_onset"], errors="coerce")
df_case["date_detection"] = pd.to_datetime(df_case["date_detection"], errors="coerce")
df_case["sex"] = df_case["sex"].replace({"M": "male", "F": "female"})
df_case["symptom_status"] = df_case["symptom_status"].str.lower().str.strip() 
df_case["specimen"] = df_case["specimen"].str.lower().str.strip()
df_case["specimen"] = df_case["specimen"].replace({
    "faeces": "stool",
    "respiratory secretions": "nasopharyngeal_swab"
})

asymptomatic_positive = df_case[
    (df_case["symptom_status"] == "asymptomatic") &
    ((df_case["ctvalue_ORF1ab"] < 40) | (df_case["ctvalue_N"] < 40))
]

asymptomatic_min_dates = asymptomatic_positive.groupby("patient_id")["date_collection"].min()

df_case["confirmation_date"] = df_case["patient_id"].map(asymptomatic_min_dates)

df_case["first_collection_date"] = df_case["patient_id"].map(asymptomatic_min_dates)

def compute_time(row):
    if row["symptom_status"] == "symptomatic" and pd.notna(row["date_onset"]):
        return (row["date_collection"] - row["date_onset"]).days
    elif row["symptom_status"] == "asymptomatic" and pd.notna(row["first_collection_date"]):
        return (row["date_collection"] - row["first_collection_date"]).days
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
            "age": float(group["age"].iloc[0]) if pd.notna(group["age"].iloc[0]) else "unknown",
            "sex": group["sex"].iloc[0] if pd.notna(group["sex"].iloc[0]) else "unknown",
            "symptom_status": group["symptom_status"].iloc[0] if pd.notna(group["symptom_status"].iloc[0]) else "unknown"
        },
        "measurements": []
    }

    for _, row in group.iterrows():
        time = int(row["time"]) if pd.notna(row["time"]) else "unknown"

        if row["symptom_status"] not in ["symptomatic", "asymptomatic"]:
            continue

        status_suffix = (
            "_symptomatic" if row["symptom_status"] == "symptomatic"
            else "_asymptomatic"
        )

        if pd.notna(row["ctvalue_ORF1ab"]):
            analyte = f"{row['specimen'].lower().replace(' ', '_')}_SARSCoV2_ORF1ab{status_suffix}"
            participant["measurements"].append({
                "analyte": analyte,
                "time": time,
                "value": float(row["ctvalue_ORF1ab"]) if row["ctvalue_ORF1ab"] < 40 else "negative"
            })

        if pd.notna(row["ctvalue_N"]):
            analyte = f"{row['specimen'].lower().replace(' ', '_')}_SARSCoV2_N{status_suffix}"
            participant["measurements"].append({
                "analyte": analyte,
                "time": time,
                "value": float(row["ctvalue_N"]) if row["ctvalue_N"] < 40 else "negative"
            })

    if participant["measurements"]:
        participants.append(participant)
```


```python
analytes = {}
for specimen in df_filtered["specimen"].dropna().unique():
    for gene in ["ORF1ab", "N"]:
        for status in ["symptomatic", "asymptomatic"]:
            analyte_key = f"{specimen.lower().replace(' ', '_')}_SARSCoV2_{gene}_{status}"
            reference_event = "symptom onset" if status == "symptomatic" else "confirmation date"
            
            description_text = (
                f"This entry describes RT-qPCR detection of SARS-CoV-2 RNA "
                f"(Ct < 40 considered positive) targeting the {gene} gene in "
                f"{specimen.lower()} from {status} individuals.\n"
            )

            analytes[analyte_key] = {
                "description": folded_str(description_text),  
                "specimen": specimen.lower(),
                "biomarker": "SARS-CoV-2",
                "gene_target": gene,
                "limit_of_quantification": "unknown",
                "limit_of_detection": 40,
                "unit": "cycle threshold",
                "reference_event": reference_event
            }


output_data = {
    "title": "SARS-CoV-2 viral shedding characteristics and potential evidence for the priority for faecal specimen testing in diagnosis",
    "doi": "10.1016/j.virusres.2020.198147",
    "description": folded_str(
        "This study investigates the shedding of SARS-CoV-2 RNA across multiple specimen types—including stool, respiratory secretions, urine, and serum—in both symptomatic and asymptomatic COVID-19 patients. It evaluates viral load dynamics and time to clearance across specimen types.\n"
    ),
    "analytes": analytes,
    "participants": participants
}

with open("yuan2021sars.yaml", "w") as outfile:
    outfile.write("# yaml-language-server:$schema=../schema.yaml\n")
    yaml.dump(output_data, outfile, default_flow_style=False, sort_keys=False)
```

```python

```
