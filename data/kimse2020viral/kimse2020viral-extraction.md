# Extraction for kimse et al. (2020)

[kimse et al. (2020)](https://www.ijidonline.com/article/S1201-9712(20)30299-X/fulltext) The authors measured SARS-CoV-2 in longitudinal throat swab samples collected from 71 COVID-19 patients between February 4 and April 7, 2020. Abundances were quantified using real-time reverse transcription polymerase chain reaction (RT-PCR). Specimens were collected from all patients at least 2 days after hospitalization and physicians checked their symptoms and signs daily. Patients who had asymptomatic carrier and incubation period were analyzed. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub). 

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
import numpy as np
from shedding_hub import folded_str
```
```python
combineddataset_df
```

```python
asymptomatic_df
```

```python
merged_df
```

```python
df
```

```python
merged_df
```

```python
combineddataset_df = pd.read_excel("CombinedDataset.xlsx") 
combineddataset_df = combineddataset_df[combineddataset_df['StudyNum'] == 4]
asymptomatic_df = pd.read_excel("asymptomatic.xlsx") 
combineddataset_df = combineddataset_df.replace({"Sex": {"M": "male", "F": "female"},
                           "PatientID":{"4-1":"1", "4-2":"2","4-3":"3"}})
# Step 2: Rename patients in asymptomatic data
asymptomatic_df['patient'] = asymptomatic_df['patient'].str.strip()
patient_mapping = {patient: f"4-{i+1}" for i, patient in enumerate(asymptomatic_df['patient'].unique())}
asymptomatic_df['patient'] = asymptomatic_df['patient'].map(patient_mapping)

# Step 4: Merge datasets on Patient and PatientID
merged_df = pd.merge(asymptomatic_df, combineddataset_df, left_on='patient', right_on='PatientID', how='outer')

# Step 5: Fill missing values with 0
#merged_df.fillna(0, inplace=True)
#df = pd.DataFrame(combineddataset_df) if not isinstance(combineddataset_df, pd.DataFrame) else combineddataset_df

participants = []

for patient_id, patient_data in merged_df("PatientID"):
#patient_id in df["PatientID"].unique():
#     patient_data = df[df["PatientID"] == patient_id]
    participant = {
        "attributes": {
            # "day": int(patient_data["Day"].iloc[0]),
            "age": float(patient_data["Age"].iloc[0]),
            "sex": str(patient_data["Sex"].iloc[0]),
            #"death": int(patient_data["Died"].iloc[0]),
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
            "value": float(row["value"]) if float(row["value"]) !=1.0 else 'negative'
        }
        participant["measurements"].append(measurement)

#不要merge，分别做成两个study，先loop前三个人用第一个analyte，再loop后十个人
    #先做，再把participants merge
    participants.append(participant)

```
Finally, the data is formatted and output as a YAML file.
```python
output_data = {
    "title": "Viral kinetics of SARS-CoV-2 in asymptomatic carriers and presymptomatic patients",
    "doi": "10.1016/j.ijid.2020.04.083",
    "description": folded_str(
        "The authors measured SARS-CoV-2 in longitudinal throat swab samples collected from 71 COVID-19 patients between February 4 and April 7, 2020..."
    ),
    "analytes": {
        "throatswab_SARSCoV2_sy": {
            "description": folded_str("SARS-CoV-2 RNA genome copy concentration in throat swab samples..."),
            "specimen": "throat_swab",
            "biomarker": "SARS-CoV-2",
            "gene_target": "RdRp",
            "limit_of_quantification": "unknown",
            "limit_of_detection": 40,
            "unit": "value",
            "reference_event": "symptom onset"
        },
         "throatswab_SARSCoV2_asy": {
            "description": folded_str("SARS-CoV-2 RNA genome copy concentration in throat swab samples..."),
            "specimen": "throat_swab",
            "biomarker": "SARS-CoV-2",
            "gene_target": "RdRp",
            "limit_of_quantification": "unknown",
            "limit_of_detection": 40,
            "unit": "value",
            "reference_event": "confirmation date"
        }
    },  
    "participants": participants #放merge之后的participants
}
with open("kimse2020viral.yaml","w") as outfile:
    outfile.write("# yaml-language-server:$schema=../.schema.yaml\n")
    yaml.dump(output_data, outfile, default_flow_style=False, sort_keys=False)
```

```python

```
