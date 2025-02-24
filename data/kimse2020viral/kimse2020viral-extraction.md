# Extraction for kimse et al. (2020)

[kim et al. (2020)](https://www.ijidonline.com/article/S1201-9712(20)30299-X/fulltext) The authors measured SARS-CoV-2 in longitudinal throat swab samples collected from 71 COVID-19 patients between February 4 and April 7, 2020. Abundances were quantified using real-time reverse transcription polymerase chain reaction (RT-PCR). Specimens were collected from all patients at least two days after hospitalization and physicians checked their symptoms and signs daily. Patients who had asymptomatic carrier and incubation period were analyzed. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub). 

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```
```python
combineddataset_df = pd.read_excel("CombinedDataset.xlsx") 
combineddataset_df = combineddataset_df[combineddataset_df['StudyNum'] == 4]
asymptomatic_df = pd.read_excel("asymptomatic.xlsx") 
combineddataset_df = combineddataset_df.replace({"Sex": {"M": "male", "F": "female"},
                           "PatientID":{"4-1":"1", "4-2":"2","4-3":"3"}})
```

```python
participants_symptomatic = []

for patient_id, patient_data in combineddataset_df.groupby("PatientID"):
    participant = {
        "attributes": {
            "age": float(patient_data["Age"].iloc[0]),
            "sex": patient_data["Sex"].iloc[0].replace("F", "female").replace("M", "male"),
        },
        "measurements": []
    }

    for _, row in patient_data.iterrows():
        measurement = {
            "analyte": "throatswab_SARSCoV2",
            "time": round(float(row["Day"])),
            # 直接在此处理负数为0
            "value": float(row["value"]) if float(row["value"]) > 0 else 0
        }
        participant["measurements"].append(measurement)

    participants_symptomatic.append(participant)

```

```python
participants_asymptomatic = []
for patient_id, patient_data in asymptomatic_df.groupby("PatientID"):
    participant = {
        "attributes": {
            "age": float(patient_data["Age"].iloc[0]),
            "sex": patient_data["Sex"].iloc[0].replace("F", "female").replace("M", "male"),
        },
        "measurements": []
    }
    for _, row in patient_data.iterrows():
        measurement = {
            "analyte": "throatswab_SARSCoV2",
            "time": round(float(row["Day"])), 
            "value": value if value > 0 else 0
        }
        participant["measurements"].append(measurement)
    
    participants_asymptomatic.append(participant)
```

```python
participants = []
participants.extend(participants_symptomatic)
participants.extend(participants_asymptomatic)
```

Finally, the data is formatted and output as a YAML file.
```python
output_data = {
    "title": "Viral kinetics of SARS-CoV-2 in asymptomatic carriers and presymptomatic patients",
    "doi": "10.1016/j.ijid.2020.04.083",
    "description": folded_str(
        "The authors measured SARS-CoV-2 in longitudinal throat swab samples collected from 71 COVID-19 patients between February 4 and April 7, 2020.\n"
    ),
    "analytes": {
        "throatswab_SARSCoV2_symptomatic": {
            "description": folded_str("SARS-CoV-2 RNA genome copy concentration in throat swab samples. Ct = 35 is the cut-off for a positive result and Ct = 40 is a negative sample; Ct = 40 was the limit of detection.\n"),
            "specimen": "throat_swab",
            "biomarker": "SARS-CoV-2",
            "gene_target": "RdRp",
            "limit_of_quantification": "unknown",
            "limit_of_detection": 40,
            "unit": "cycle threshold",
            "reference_event": "symptom onset"
        },
         "throatswab_SARSCoV2_asymptomatic": {
            "description": folded_str("SARS-CoV-2 RNA genome copy concentration in throat swab samples.Ct = 35 is the cut-off for a positive result and Ct = 40 is a negative sample; Ct = 40 was the limit of detection.\n"),
            "specimen": "throat_swab",
            "biomarker": "SARS-CoV-2",
            "gene_target": "RdRp",
            "limit_of_quantification": "unknown",
            "limit_of_detection": 40,
            "unit": "cycle threshold",
            "reference_event": "confirmation date"
        }
    },  
    "participants": participants 
}
with open("kimse2020viral.yaml","w") as outfile:
    outfile.write("# yaml-language-server:$schema=../.schema.yaml\n")
    yaml.dump(output_data, outfile, default_flow_style=False, sort_keys=False)
```

```python
#make sure branch
#git pull
#result
#push
```
