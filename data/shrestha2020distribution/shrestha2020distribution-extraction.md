# Extraction for Shrestha et al. (2020)

[Shrestha et al. (2020)](https://doi.org/10.1093/cid/ciaa886) evaluated the transmission potential of COVID-19 by examining viral load over time. Over six weeks, 230 healthcare personnel underwent 528 tests at the Cleveland Clinic. Cycle threshold (Ct) values were obtained using RT-PCR targeting the N gene, and viral loads were calculated. Data were obtained from the combined dataset in the supplementary materials of Challenger et al. BMC Medicine (2022) 20:25 (https://doi.org/10.1186/s12916-021-02220-0).

```python
#import modules;
import json
import pandas as pd 
import numpy as np
import jsonschema
import yaml

# Functions to add folded blocks and literal blocks;
class folded_str(str): pass
class literal_str(str): pass

def folded_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='>')
def literal_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(folded_str, folded_str_representer)
yaml.add_representer(literal_str, literal_str_representer)

# Load dataset
df = pd.read_excel('CombinedDataset.xlsx', sheet_name='Viral_Load')
shrestha2020 = df[df['StudyNum'] == 12].copy()
columns_to_drop = ['Estimated', 'SevMax', 'Sev1st', 'Died', 'Ctvalue', 'SevMax3']
shrestha2020 = shrestha2020.drop(columns=columns_to_drop)
```

```python
# Group by participant and extract measurements
participants = []

for patient_id, group in shrestha2020.groupby('PatientID'):
    participant = {
        'attributes': {
            'age': int(group['Age'].iloc[0]),
            'sex': 'female' if group['Sex'].iloc[0] == 'F' else 'male'
        },
        'measurements': []
    }

    for _, row in group.iterrows():
        if row['value'] == 1:
            value = "negative"
        else:
            value = row['value']
        measurementN = {
            'analyte': 'swab_SARSCoV2_N',
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurementN)
        
    participants.append(participant)

```

The data is formatted and output as a YAML file.

```python
shrestha2020distribution = dict(title="Distribution of Transmission Potential During Nonsevere COVID-19 Illness",
        doi="10.1093/cid/ciaa886",
        description=folded_str('This study evaluated the transmission potential of COVID-19 by examining viral load over time. Over six weeks, 230 healthcare personnel underwent 528 tests at the Cleveland Clinic. Cycle threshold (Ct) values were obtained using RT-PCR targeting the N gene, and viral loads were calculated. Data were obtained from the combined dataset in the supplementary materials of Challenger et al. BMC Medicine (2022) 20:25 (https://doi.org/10.1186/s12916-021-02220-0).\n'),
        analytes=dict(swab_SARSCoV2_N=dict(description=folded_str("Cycle threshold (Ct) values were quantified using RT-PCR targeting the N gene in nasopharyngeal swab samples. Ct values were then converted into copies per mL, with a detection limit of 20 copies/mL. Viral load calculations were based on the minimum detectable viral load (MDVL) and Ct values.\n"),
                                    limit_of_quantification="unknown",
                                    limit_of_detection=20,
                                    specimen="nasopharyngeal_swab", 
                                    biomarker="SARS-CoV-2", 
                                    gene_target="N", 
                                    unit="gc/mL",
                                    reference_event="symptom onset",)),
            participants = participants
                                    )

with open("shrestha2020distribution.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(shrestha2020distribution, outfile, default_flow_style=False, sort_keys=False)
outfile.close()

```
