# Extraction for Peiris et al. (2003)

[Peiris et al. (2003)](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(03)13412-5/fulltext) measured SARS detected by real-time reverse transcriptase PCR in nasopharyngeal samples from 14 patients who admitted to the United Christian Hospital from the Amoy Gardens housing estate who fulfilled the modified WHO definition of SARS, up to days 5, 10 and 15 after symptom onset. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/peiris2003clinical). 

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```

Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/peiris2003clinical), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
Peiris2003 = pd.read_csv("data.csv")
# Map the 'PatientID' column to the new 'ID' column
patient_id_mapping = {pid: i+1 for i, pid in enumerate(sorted(Peiris2003 ['PatientID'].unique()))}
Peiris2003['ID'] = Peiris2003['PatientID'].map(patient_id_mapping)
Peiris2003["type"] = "NPS"


participant_list = []

for i in pd.unique(Peiris2003["ID"]):
    patient_data = Peiris2003[Peiris2003["ID"] == i]


    measurements = []
    for _, row in patient_data.iterrows():
        try:
            # Handle missing values safely
            if pd.isna(row['value']):
                value = 'negative'
            else:
                value = 10 ** float(row['value']) 

            # Special condition for 1.0
            if value == 1.0:
                value = 'negative'

        except ValueError:
            value = 'negative'

        # Append measurements separately for each sample type
        if row['type'] == 'NPS':
            measurements.append({
                "analyte": "NPS_SARS", 
                "time": int(row['Day']), 
                "value": value})


    participant_dict = {"measurements": measurements}
    participant_list.append(participant_dict)


```

Finally, the data is formatted and output as a YAML file.

```python
peiris2003 = dict(title="Clinical progression and viral load in a community outbreak of coronavirus-associated SARS pneumonia: a prospective study",
            doi="10.1016/S0140-6736(03)13412-5",
            description=folded_str("This study measured SARS detected by real-time reverse transcriptase PCR in nasopharyngeal samples from 14 patients who admitted to the United Christian Hospital from the Amoy Gardens housing estate who fulfilled the modified WHO definition of SARS, up to days 5, 10 and 15 after symptom onset.\n"),
            analytes=dict(NPS_SARS=dict(description=folded_str("This analyte indicates the detection of SARS RNA in nasopharyngeal, stool and faecal swabs, but only the nasopharyngeal swab values are presented.\n"),
                          specimen="nasopharyngeal_swab",
                          biomarker="SARS",
                          gene_target="unkonw",
                          limit_of_quantification="unknown",
                          limit_of_detection="unknown",
                          unit="gc/mL",
                          reference_event="symptom onset",
        )
    ),
    participants=participant_list,
)


with open("peiris2003clinical.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")
    yaml.dump(peiris2003, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close() 

```