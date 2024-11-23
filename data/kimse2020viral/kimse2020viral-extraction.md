# Extraction for Kim et al. (2020)

[Kim et al. (2020)](https://www.ijidonline.com/article/S1201-9712(20)30299-X/fulltext) The authors measured SARS-CoV-2 in longitudinal throat swab samples collected from 71 COVID-19 patients between February 4 and April 7, 2020. Abundances were quantified using real-time reverse transcription polymerase chain reaction (RT-PCR). Specimens were collected from all patients at least 2 days after hospitalization and physicians checked their symptoms and signs daily. Patients who had asymptomatic carrier and incubation period were analyzed. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/kimse2020viral).

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```
```python
#load the data
data = pd.read_excel("CombinedDataset.xlsx")
data = data[data['StudyNum'] == 4]
data = data.replace({"Sex": {"M": "male", "F": "female"},
                    "PatientID":{"4-1":"1", "4-2":"2","4-3":"3"}})

#some data cleaning to match the schema;
df = pd.DataFrame(data) if not isinstance(data, pd.DataFrame) else data

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
            "value": float(row["value"]) if float(row["value"]) !=1.0 else 'negative'
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
with open("kimse2020viral.yaml","w") as outfile:
    yaml.dump(output_data, outfile, default_flow_style=False, allow_unicode=True, sort_keys=False)

```
