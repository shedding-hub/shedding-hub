# Extraction for The COVID-19 Investigation Team (2020)

[The COVID-19 Investigation Team (2020)](https://www.nature.com/articles/s41591-020-0877-5) studied clinical and virologic characteristics of the first 12 patients with coronavirus disease 2019 (COVID-19) in the United States. Respiratory, stool, serum, and urine specimens were submitted for SARS-CoV-2 real-time reverse-transcription polymerase chain reaction (rRT-PCR) testing, viral culture, and whole genome sequencing. Only nasopharyngeal samples were included in this dataset. Data for the nasopharyngeal swab results were obtained from the combined dataset in the supplementary materials of [Challenger et al. (2022)](https://doi.org/10.1186/s12916-021-02220-0). Attributes of hospitalized patients were obtained from the original paper. 

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
from shedding_hub import folded_str
```

Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/team2020clinical), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml). 

```python
#load the data;
df = pd.read_excel('CombinedDataset.xlsx', sheet_name='Viral_Load')
team2020 = df[df['StudyNum'] == 7].copy()
columns_to_drop = ['Estimated', 'SevMax', 'Sev1st', 'Died', 'SevMax3', 'Age', 'Sex']
team2020 = team2020.drop(columns=columns_to_drop)

# Information from the original paper.
patient_info = {
    'PatientID': ['7-6', '7-7', '7-8', '7-9', '7-10', '7-11', '7-12'],
    'Age': ['30-39', '60-69', '60-69', '30-39', '50-59', '50-59', '50-59'],
    'Sex': ['male', 'female', 'male', 'male', 'male', 'male', 'female']
}
patient_df = pd.DataFrame(patient_info)
team2020 = team2020.merge(patient_df, on='PatientID', how='left')

participants = []

# Group by participant and extract measurements

for patient_id, group in team2020.groupby('PatientID'):
    # Check if both 'Age' and 'Sex' are not NaN
    if pd.notna(group['Age'].iloc[0]) and pd.notna(group['Sex'].iloc[0]):
        participant = {
            'attributes': {
                'age_group': group['Age'].iloc[0],
                'sex': group['Sex'].iloc[0]
            },
            'measurements': []
        }
    else:
        # If age or sex is NaN, exclude 'attributes' from the participant dictionary
        participant = {
            'measurements': []
        }

    for _, row in group.iterrows():
        if row['value'] == 1:
            value = "negative"
        else:
            value = row['value']
        measurement_VL = {
            'analyte': 'naso_swab_SARSCoV2_N_VL',
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurement_VL)
    
    for _, row in group.iterrows():
        value = row['Ctvalue']
        measurement_Ct = {
            'analyte': 'naso_swab_SARSCoV2_N_Ct',
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurement_Ct)
        
    participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.

```python
team2020clinical = dict(title="Clinical and virologic characteristics of the first 12 patients with coronavirus disease 2019 (COVID-19) in the United States",
               doi="10.1038/s41591-020-0877-5",
               description=folded_str('Respiratory, stool, serum, and urine specimens were submitted for SARS-CoV-2 real-time reverse-transcription polymerase chain reaction (rRT-PCR) testing, viral culture, and whole-genome sequencing. Only nasopharyngeal samples were included in this dataset. Data for the nasopharyngeal swab results were obtained from the combined dataset in the supplementary materials of Challenger et al. (2022), while attributes of hospitalized patients were obtained from the original paper.\n'),
               analytes=dict(naso_swab_SARSCoV2_N_VL=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentrations in nasopharyngeal swab samples were obtained from the combined dataset in the supplementary materials of Challenger et al. (2022). The viral load values were estimated using an average standard curve.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection="unknown",
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",),
                            naso_swab_SARSCoV2_N_Ct=dict(description=folded_str("Cycle threshold (Ct) values were quantified using rRT-PCR targeting the N gene in nasopharyngeal swab samples.\n"),
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

with open("team2020clinical.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(team2020clinical, outfile, default_flow_style=False, sort_keys=False)
```
