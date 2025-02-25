# Extraction for Gautret et al. (2020)

[Gautret et al. (2020)](https://www.sciencedirect.com/science/article/pii/S0924857920300996?via%3Dihub) evaluated the effect of hydroxychloroquine on respiratory viral loads. Six patients were asymptomatic, 22 exhibited upper respiratory tract infection symptoms, and eight showed lower respiratory tract infection symptoms. Only 19 patients with a known onset of symptoms were included in the analysis below. Data for the nasopharyngeal wab results were obtained from the combined dataset in the supplementary materials of [Challenger et al. (2022)](https://doi.org/10.1186/s12916-021-02220-0). Attributes of drug treatments were sourced from the supplementary materials of the original paper.

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
from shedding_hub import folded_str
```

Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/gautret2020hydroxychloroquine), will be loaded and cleaned to match the [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml). 

```python
df = pd.read_excel('CombinedDataset.xlsx', sheet_name='Viral_Load')
gautret2020 = df[df['StudyNum'] == 5].copy()
columns_to_drop = ['Estimated', 'SevMax', 'Sev1st', 'Died', 'SevMax3']
gautret2020 = gautret2020.drop(columns=columns_to_drop)

# Information from the original paper.
treatment_info = {
    'Hydroxychloroquine': ['no', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes','yes'],
    'Azithromycin': ['no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'yes', 'yes', 'yes', 'yes', 'yes','yes'],
    'PatientID': ['5-1', '5-11', '5-12', '5-13', '5-14', '5-15', '5-16', '5-17', '5-18', '5-19', '5-20', '5-21', '5-22', '5-23', '5-24', '5-25', '5-26', '5-27', '5-28']
}
treatment_df = pd.DataFrame(treatment_info)
merged_df = pd.merge(gautret2020, treatment_df, on='PatientID', how='left')

participants = []

# Group by participant and extract measurements
for patient_id, group in merged_df.groupby('PatientID'):
    participant = {
        'attributes': {
            'age': int(group['Age'].iloc[0]),
            'sex': 'female' if group['Sex'].iloc[0] == 'F' else 'male',
            'Hydroxychloroquine_treatment': True if group['Hydroxychloroquine'].iloc[0] == 'yes' else False,
            'Azithromycin_treatment':True if group['Azithromycin'].iloc[0] == 'yes' else False
        },
        'measurements': []
    }
    
    for _, row in group.iterrows():
        if row['value'] == 1:
            value = "negative"
        else:
            value = row['Ctvalue']
        measurementN = {
            'analyte': 'nasopharyngeal_swab_SARSCoV2_Ct',
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurementN)
        
    participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.

```python
gautret2020hydroxychloroquine = dict(title="Hydroxychloroquine and azithromycin as a treatment of COVID-19:results of an open-label non-randomized clinical trial",
               doi="10.1016/j.ijantimicag.2020.105949",
               description=folded_str('This study evaluated the effect of hydroxychloroquine on respiratory viral loads. Six patients were asymptomatic, 22 exhibited upper respiratory tract infection symptoms, and eight showed lower respiratory tract infection symptoms. Only 19 patients with a known onset of symptoms were included in the analysis below. Data for the nasopharyngeal swab results were obtained from the combined dataset in the supplementary materials of Challenger et al. (2022). Attributes of drug treatments were sourced from the supplementary materials of the original paper.\n'),
               analytes=dict(nasopharyngeal_swab_SARSCoV2_Ct=dict(description=folded_str("Cycle threshold (Ct) values targeting the E gene in nasopharyngeal swab samples.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection="unknown",
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        unit="cycle threshold",
                                        reference_event="symptom onset",)
                            ),
                participants = participants
                                        )


with open("gautret2020hydroxychloroquine.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(gautret2020hydroxychloroquine, outfile, default_flow_style=False, sort_keys=False)
```
