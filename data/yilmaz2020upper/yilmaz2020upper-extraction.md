# Extraction for ylimaz2020 et al. (2020)

[Ylimaz et al. (2020)](https://doi.org/10.1093/infdis/jiaa632) reported longitudinal viral RNA loads from the nasopharynx/throat in patients with mild and severe/critical coronavirus disease 2019 (COVID-19). The authors also investigated whether the duration of symptoms correlated with the duration of viral RNA shedding. A total of 56 patients were included. The raw data was obtained from [Challenger et al. (2020)](https://doi.org/10.1186/s12916-021-02220-0) and validated with Figure 1 in [Ylimaz et al. (2020)](https://doi.org/10.1093/infdis/jiaa632). The data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/yilmaz2020upper). 

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
    **{f"14-{i}": str(i) for i in range(1, 55)}  
})

```

```python
# Ensure 'data' is a DataFrame
df = pd.DataFrame(data) if not isinstance(data, pd.DataFrame) else data

# Strip spaces from column names to avoid issues with whitespace
df.columns = df.columns.str.strip()

# Verify the 'PatientID' column exists
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

        
    # Iterate through each patient's data
for _, row in patient_data.iterrows():
    # Add viral load (VL) measurement
    if pd.notna(row['value']):
        measurement_vl = {
            "analyte": "throatswab_SARSCoV2",
            "time": int(row["Day"]),
            "value": "negative" if row['value'] == 1.0 else row['value']  # 1.0 -> negative, 其他保留
        }
        participant['measurements'].append(measurement_vl)

    # Add Ct value measurement
    if pd.notna(row['Ctvalue']):
        measurement_ct = {
            "analyte": "throatswab_SARSCoV2",
            "time": int(row["Day"]),
            "value": "negative" if row['Ctvalue'] == 40.0 else row['Ctvalue']  # 40.0 -> negative, 其他保留
        }
        participant['measurements'].append(measurement_ct)


# Append participant data after processing each patient
participants.append(participant)


```

Finally, the data is formatted and output as a YAML file.
```python
output_data = {
    "title": "Upper Respiratory Tract Levels of Severe Acute Respiratory Syndrome Coronavirus 2 RNA and Duration of Viral RNA Shedding Do Not Differ Between Patients With Mild and Severe/Critical Coronavirus Disease 2019",
    "doi": "10.1093/infdis/jiaa632",
    "description": folded_str(
        "The authors reported longitudinal viral RNA loads from the nasopharynx/throat in patients with mild and severe/critical coronavirus disease 2019 (COVID-19). They also investigated whether the duration of symptoms correlated with the duration of viral RNA shedding. A total of 56 patients were included.\n"
    ),
    "analytes": {
        "throatswab_SARSCoV2": {
            "description": folded_str(
                "The author collected serial upper respiratory tract samples (one nasopharyngeal swab and one throat swab were put in a single collection tube with 1 mL of trans- port medium) for real-time PCR of SARS-CoV-2 RNA for all patients.\n"
            ),
            "specimen": ["nasopharyngeal_swab", "throat_swab"],
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
with open("yilmaz2020upper.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(output_data, outfile, default_flow_style=False, allow_unicode=True, sort_keys=False)
```
