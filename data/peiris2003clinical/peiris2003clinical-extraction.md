# Extraction for Peiris et al. (2003)

[Peiris et al. (2003)](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(03)13412-5/fulltext)measured severe acute respiratory syndrome (SARS) detected by real-time reverse transcriptase PCR in nasopharynx, stool, and urine samples from day 10 to 21 after onset of symptoms. Currently, the dataset includes viral shedding data from the nasopharynx on days 5, 10, and 15 after symptom onset of 14 patients. Demographic information (e.g., age, sex, etc.) is not included. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/peiris2003clinical). 



First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
from shedding_hub import folded_str
```

We extracted the raw data using [automeris.io](https://automeris.io/) from Figure 4 (see below) in [Peiris et al. (2003)](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(03)13412-5/fulltext). The extracted data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/peiris2003clinical), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

[!image](patient_data.png)

```python
Peiris2003 = pd.read_csv("data.csv")

patient_id_mapping = {pid: i+1 for i, pid in enumerate(sorted(Peiris2003 ['PatientID'].unique()))}
# Map the 'PatientID' column to the new 'ID' column
Peiris2003['ID'] = Peiris2003['PatientID'].map(patient_id_mapping)

Peiris2003["type"] = "NPS"


participant_list = []

for i in pd.unique(Peiris2003["ID"]):
    patient_data = Peiris2003[Peiris2003["ID"] == i]

    measurements = []
    for _, row in patient_data.iterrows():
        value = 10 ** float(row['value'])

        # Append measurements separately for each sample type
        if row['type'] == 'NPS':
            measurements.append({
                "analyte": "NPS_SARS", 
                "time": round(row['Day']), 
                "value": value})


    participant_dict = {"measurements": measurements}
    participant_list.append(participant_dict)

```
Finally, the data is formatted and output as a YAML file.

```python
peiris2003 = dict(title="Clinical progression and viral load in a community outbreak of coronavirus-associated SARS pneumonia: a prospective study",
            doi="10.1016/S0140-6736(03)13412-5",
            description=folded_str("This study measured SARS detected by real-time reverse transcriptase PCR in nasopharyngeal samples from 14 patients who admitted to the United Christian Hospital from the Amoy Gardens housing estate who fulfilled the modified WHO definition of SARS, on days 5, 10, and 15 after symptom onset.\n"),
            analytes=dict(NPS_SARS=dict(description=folded_str("This analyte indicates the detection of SARS RNA in nasopharyngeal aspirates.\n"),
                          specimen="nasopharyngeal_swab", 
                          biomarker="SARS",
                          gene_target="unkonw",
                          limit_of_quantification="unknown",
                          limit_of_detection="unknown",  # Essentially, 1 swab sample was transformed to 3 mL of liquid sample, from the original paper (https://pmc.ncbi.nlm.nih.gov/articles/PMC7108127/#R4).
                          unit="gc/mL",
                          reference_event="symptom onset",
        )
    ),
    participants=participant_list,
)


with open("peiris2003clinical.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(peiris2003, outfile, default_style=None, default_flow_style=False, sort_keys=False)
```