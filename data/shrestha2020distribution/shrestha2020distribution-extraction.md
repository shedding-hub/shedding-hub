# Extraction for Shrestha et al. (2020)

[Shrestha et al. (2020)](https://doi.org/10.1093/cid/ciaa886) evaluated the transmission potential of COVID-19 by examining viral load over time. Over six weeks, 230 healthcare personnel underwent 528 tests at the Cleveland Clinic. Cycle threshold (Ct) values were obtained using RT-PCR targeting the N gene, and viral loads were calculated. Data were obtained from the combined dataset in the supplementary materials of Challenger et al. BMC Medicine (2022) 20:25 (https://doi.org/10.1186/s12916-021-02220-0).

```python
import pandas as pd
import yaml
from shedding_hub import folded_str

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
        measurement_VL = {
            'analyte': 'swab_SARSCoV2_N_VL',
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurement_VL)
    
    for _, row in group.iterrows():
        value = row['Ctvalue']
        measurement_Ct = {
            'analyte': 'swab_SARSCoV2_N_Ct',
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurement_Ct)
        
    participants.append(participant)
```

The data is formatted and output as a YAML file.

```python
shrestha2020distribution = dict(title="Distribution of Transmission Potential During Nonsevere COVID-19 Illness",
               doi="10.1093/cid/ciaa886",
               description=folded_str('This study evaluated the transmission potential of COVID-19 by examining viral load over time. Over six weeks, 230 healthcare personnel underwent 528 tests at the Cleveland Clinic. Cycle threshold (Ct) values were obtained using RT-PCR targeting the N gene, and viral loads were calculated. Data were obtained from the combined dataset in the supplementary materials of Challenger et al. BMC Medicine (2022) 20:25 (https://doi.org/10.1186/s12916-021-02220-0).\n'),
               analytes=dict(swab_SARSCoV2_N_VL=dict(description=folded_str("Cycle threshold (Ct) values were quantified using RT-PCR targeting the N gene in nasopharyngeal swab samples. Ct values were then converted into copies per mL. Viral load calculations were based on an average standard curve and Ct values obtained from the combined dataset in the supplementary materials of Challenger et al. (2022). Using commercially available plasmids that contained the nucleocapsid (N) gene, the LOD was found to be 20 copies per uL for upper respiratory specimens. We cannot determine the LOD for raw samples.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection="unknown",
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",),
                            swab_SARSCoV2_N_Ct=dict(description=folded_str("Cycle threshold (Ct) values were quantified using RT-PCR targeting the N gene in nasopharyngeal swab samples.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection="unknown",
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N", 
                                        unit="cycle threshold",
                                        reference_event="symptom onset",)
                                        ),
                participants = participants
                                        )

with open("shrestha2020distribution.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(shrestha2020distribution, outfile, default_flow_style=False, sort_keys=False)
outfile.close()
```
