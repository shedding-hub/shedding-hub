# Extraction for obara2008single et al. (2008)

[Obara et al. (2008)](https://journals.asm.org/doi/10.1128/jcm.01932-07) reported longitudinal norovirus RNA shedding in two food handlers—one symptomatic and one asymptomatic—during a foodborne outbreak investigation. The study used real-time reverse transcription PCR (RT-PCR) to quantify norovirus RNA in stool samples collected at multiple time points over a period of more than two months. One subject (Employee A) reported gastrointestinal symptoms and exhibited prolonged shedding, while the other (Employee B) remained asymptomatic but also showed detectable viral RNA levels. The authors reported both qualitative PCR results (positive/negative) and quantitative viral loads (log₁₀ copies per gram of feces), allowing detailed tracking of shedding dynamics over time. Data were extracted directly from Tables and Figures (notably Figure 1) of the original publication and structured using standardized analyte and participant formats. Sample type was recorded as stool in the standardized dataset. The data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/obara2008single).

First, we `import` python modules needed:
```python
import yaml
import pandas as pd
import numpy as np
from shedding_hub import folded_str
from datetime import datetime
```
```python
# Load log-transformed norovirus copy number data
df_log = pd.read_excel("log10_copy_number_data.xlsx")
# Load corresponding PCR result data
df_pcr = pd.read_excel("employee_pcr_data.xlsx")

# Rename columns for clarity and consistency
df_log.columns = ["PatientID", "Group", "Date", "Log10CopyNumber"]
df_pcr.columns = ["Employee", "PatientID", "Date", "PCR_Result_Raw", "PCR_Value", "PCR_Qualitative"]

# Extract shared participant ID (a or b) from sample-level PatientID (e.g., a-1, b-3)
df_log["PersonID"] = df_log["PatientID"].str.extract(r"(a|b)")
df_pcr["PersonID"] = df_pcr["PatientID"].str.extract(r"(a|b)")

# Merge the log data and PCR data on PatientID and Date, preserving PersonID from df_log
merged_df = pd.merge(
    df_log[["PatientID", "PersonID", "Date", "Log10CopyNumber"]],
    df_pcr[["PatientID", "Date", "PCR_Qualitative"]],
    on=["PatientID", "Date"],
    how="left"
)

# Convert Date to datetime format for subsequent time calculations
merged_df["Date"] = pd.to_datetime(merged_df["Date"], errors="coerce")
```

```python
participants = []

# Group merged sample-level data by participant ID ('a' or 'b')
for person_id, person_data in merged_df.groupby("PersonID"):
    # Determine group type based on symptomatic (a) or asymptomatic (b)
    group_type = True if person_id == "a" else False
    # Initialize participant record with group attribute and empty measurement list
    participant = {
        "attributes": {"symptomatic": group_type},
        "measurements": []
    }

    # Sort samples by date and establish earliest date as Day 0 for this participant
    person_data = person_data.sort_values("Date")
    reference_date = person_data["Date"].min()

    # Iterate over all samples from this participant
    for _, row in person_data.iterrows():
        if pd.isna(row["Date"]):
            continue 
 
        # Calculate time offset (in days) from reference_date
        day_offset = (row["Date"] - reference_date).days

        # If quantitative viral load data is available, record it
        if pd.notna(row["Log10CopyNumber"]):
           val = float(row["Log10CopyNumber"])

   
           if val > 100000:  
               val = np.log10(val)

           participant["measurements"].append({
               "analyte": "stool_Norovirus_copies",
               "time": day_offset,
               "value": 10**(val)
           })

        # Append qualitative PCR result (binary: positive/negative), not subject to log10 transformation
        if pd.notna(row["PCR_Qualitative"]) and row["PCR_Qualitative"] in ["positive", "negative"]:
            participant["measurements"].append({
                "analyte": "norovirus_presence_qualitative",
                "time": day_offset,
                "value": row["PCR_Qualitative"]
            })

    # Only include this participant if at least one measurement is recorded
    if participant["measurements"]:
        participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.
```python
output_data = {
    "title": "Single Base Substitutions in the Capsid Region of the Norovirus Genome during Viral Shedding in Cases of Infection in Areas Where Norovirus Infection Is Endemic",
    "doi": "10.1128/jcm.01932-07",
    "description": folded_str(
        "This study investigates the duration of norovirus RNA shedding in two food handlers-"
        "one symptomatic and one asymptomatic-using repeated real-time reverse transcription-PCR "
        "tests on fecal specimens. Norovirus RNA was detected and quantified over time, and the duration "
        "of viral shedding was compared between symptomatic and asymptomatic individuals.\n"
    ),
   "analytes": {
        "stool_Norovirus_copies": {
            "description": folded_str(
                "Quantitative measurement of norovirus RNA in stool samples using real-time RT-PCR. "
                "Results are expressed as copies per gram of feces.\n"
            ),
            "specimen": "stool",
            "biomarker": "norovirus",
            "genotype": "GII.4",
            "limit_of_quantification": "unknown",
            "limit_of_detection": 10000,
            "unit": "gc/wet gram",
            "reference_event": "confirmation date"
        },
        "norovirus_presence_qualitative": {
            "description": folded_str(
                "Detection of norovirus RNA in stool samples using real-time RT-PCR. "
                "Results are recorded as 'positive' or 'negative'.\n"
            ),
            "specimen": "stool",
            "biomarker": "norovirus",
            "genotype": "GII.4",
            "limit_of_quantification": "unknown",
            "limit_of_detection": "unknown",
            "unit": None,
            "reference_event": "confirmation date"
        },
        "norovirus_ct_value": {
            "description": folded_str(
                "Cycle threshold (Ct) values from real-time RT-PCR detection of norovirus RNA in stool samples. "
                "Lower Ct values indicate higher viral load.\n"
            ),
            "specimen": "stool",
            "biomarker": "norovirus",
            "genotype": "GII.4",
            "limit_of_quantification": "unknown",
            "limit_of_detection": "unknown",
            "unit": "cycle threshold",
            "reference_event": "confirmation date"
        }
    },
    "participants": participants
}

with open("obara2008single.yaml", "w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(output_data, outfile, default_flow_style=False, sort_keys=False)
```
