# Extraction for obara2008single et al. (2008)

[obara et al. (2008)](https://journals.asm.org/doi/10.1128/jcm.01932-07) reported longitudinal norovirus RNA shedding in two food handlers—one symptomatic and one asymptomatic—during a foodborne outbreak investigation. The study used real-time reverse transcription PCR (RT-PCR) to quantify norovirus RNA in stool samples collected at multiple time points over a period of more than two months. One subject (Employee A) reported gastrointestinal symptoms and exhibited prolonged shedding, while the other (Employee B) remained asymptomatic but also showed detectable viral RNA levels. The authors reported both qualitative PCR results (positive/negative) and quantitative viral loads (log₁₀ copies per gram of feces), allowing detailed tracking of shedding dynamics over time. Data were extracted directly from Tables and Figures (notably Figure 1) of the original publication and structured using standardized analyte and participant formats. Sample type was recorded as stool in the standardized dataset. The data is stored at[Shedding Hub]

First, we `import` python modules needed:
```python
import yaml
import pandas as pd
from shedding_hub import folded_str
from datetime import datetime
```
```python
df_log = pd.read_excel("log10_copy_number_data.xlsx")
df_pcr = pd.read_excel("employee_pcr_data.xlsx")

df_log.columns = ["PatientID", "Group", "Date", "Log10CopyNumber"]
df_pcr.columns = ["Employee", "PatientID", "Date", "PCR_Result_Raw", "PCR_Value", "PCR_Qualitative"]
df_log["PersonID"] = df_log["PatientID"].str.extract(r"(a|b)")
df_pcr["PersonID"] = df_pcr["PatientID"].str.extract(r"(a|b)")

merged_df = pd.merge(
    df_log[["PatientID", "PersonID", "Date", "Log10CopyNumber"]],
    df_pcr[["PatientID", "Date", "PCR_Qualitative"]],
    on=["PatientID", "Date"],
    how="left"
)
merged_df["Date"] = pd.to_datetime(merged_df["Date"], errors="coerce")
```

```python
participants = []

for person_id, person_data in merged_df.groupby("PersonID"):
    group_type = "symptomatic" if person_id == "a" else "asymptomatic"
    participant = {
        "attributes": {"group": group_type},
        "measurements": []
    }

    person_data = person_data.sort_values("Date")
    reference_date = person_data["Date"].min()

    for _, row in person_data.iterrows():
        if pd.isna(row["Date"]):
            continue 

        day_offset = (row["Date"] - reference_date).days

        if pd.notna(row["Log10CopyNumber"]):
            participant["measurements"].append({
                "analyte": "stool_Norovirus_log10copies",
                "time": day_offset,
                "value": float(row["Log10CopyNumber"])
            })

        if pd.notna(row["PCR_Qualitative"]) and row["PCR_Qualitative"] in ["positive", "negative"]:
            participant["measurements"].append({
                "analyte": "norovirus_presence_qualitative",
                "time": day_offset,
                "value": row["PCR_Qualitative"]
            })

    if participant["measurements"]:
        participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.
```python
output_data = {
    "title": "Duration of Norovirus Shedding in Symptomatic and Asymptomatic Food Handlers",
    "doi": "10.1128/jcm.01932-07",
    "description": folded_str(
        "This study investigates the duration of norovirus RNA shedding in two food handlers—"
        "one symptomatic and one asymptomatic—using repeated real-time reverse transcription-PCR "
        "tests on fecal specimens. Norovirus RNA was detected and quantified over time, and the duration "
        "of viral shedding was compared between symptomatic and asymptomatic individuals.\n"
    ),
    "analytes": {
        "stool_Norovirus_log10copies": {
            "description": folded_str(
                "Quantitative measurement of norovirus RNA in stool samples using real-time RT-PCR. "
                "Results are expressed as log10 copies per gram of feces.\n"
            ),
            "specimen": "stool",
            "biomarker": "sapovirus",
            "gene_target": "GII",
            "limit_of_quantification": "unknown",
            "limit_of_detection": "unknown",
            "unit": "gc/wet gram",
            "reference_event": "confirmation date"
        },
        "norovirus_presence_qualitative": {
            "description": folded_str(
                "Qualitative detection of norovirus RNA in stool samples using real-time RT-PCR. "
                "Results are recorded as 'positive' or 'negative'.\n"
            ),
            "specimen": "stool",
            "biomarker": "sapovirus",
            "gene_target": "GII",
            "limit_of_quantification": "unknown",
            "limit_of_detection": "unknown",
            "unit": "gc/wet gram",
            "reference_event": "confirmation date"
        }
    },
    "participants": participants
}

with open("obara2008single.yaml", "w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(output_data, outfile, default_flow_style=False, allow_unicode=True, sort_keys=False)
```
