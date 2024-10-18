# Extraction of Young et al. (2020) using the data published by Challenger et al. (2022) and Goyal et al. (2020)

[Young et al. 2020]((https://jamanetwork.com/journals/jama/fullarticle/2762688)) collect clinical, laboratory, and radiologic data, including PCR cycle threshold values from nasopharyngeal swabs and viral shedding in blood, urine, and stool. The clinical course was summarized, including the requirement for supplemental oxygen, intensive care, and the use of empirical treatment with lopinavir-ritonavir. The numbers of positive stool, blood, and urine samples were small. Data for the nasopharyngeal swab results were obtained from the combined dataset in the supplementary materials of Challenger et al. BMC Medicine (2022) 20:25 (https://doi.org/10.1186/s12916-021-02220-0). The standard curve was calculated based on the concentration provided by Goyal et al. (2020)(https://www.science.org/doi/10.1126/sciadv.abc7112). 

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
import numpy as np
from shedding_hub import folded_str, literal_str
```

Raw data ([CombinedDataset](https://github.com/shedding-hub/shedding-hub/blob/main/data/young2020epidemiologic/CombinedDataset.xlsx)), which is stored on [Shedding Hub](https://github.com/shedding-hub), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
df = pd.read_excel(os.path.join(current_directory, 'Challenger2022/CombinedDataset.xlsx'), sheet_name='Viral_Load')
young2020 = df[df['StudyNum'] == 6].copy()
columns_to_drop = ['Estimated', 'SevMax', 'Sev1st', 'Died', 'SevMax3']
young2020 = young2020.drop(columns=columns_to_drop)
# Transform Ct value to gc/swab using obtained standard curve
young2020['concentration'] = np.where(young2020['Ctvalue'] == 38, 1, 10 ** (14.1147 - 0.3231 * young2020['Ctvalue']))
participants = []
# Group by participant and extract measurements
for patient_id, group in young2020.groupby('PatientID'):
    participant = {
        'measurements': []
    }

    for _, row in group.iterrows():
        if row['concentration'] == 1:
            value = "negative"
        else:
            value = row['concentration']
        measurementN = {
            'analyte': 'swab_SARSCoV2_N',
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurementN)
        
    participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.

```python
young2020epidemiologic = dict(title="Epidemiologic Features and Clinical Course of Patients Infected With SARS-CoV-2 in Singapore",
               doi="10.1001/jama.2020.3204",
               description=folded_str('Clinical, laboratory, and radiologic data were collected, including PCR cycle threshold values from nasopharyngeal swabs and viral shedding in blood, urine, and stool. The clinical course was summarized, including the requirement for supplemental oxygen, intensive care, and the use of empirical treatment with lopinavir-ritonavir. The numbers of positive stool, blood, and urine samples were small. Data for the nasopharyngeal swab results were obtained from the combined dataset in the supplementary materials of Challenger et al. BMC Medicine (2022) 20:25 (https://doi.org/10.1186/s12916-021-02220-0). The standard curve was calculated based on the concentration provided by Goyal et al. (2020).\n'),
               analytes=dict(swab_SARSCoV2_N=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration in nasopharyngeal swab samples. The unit of concentration were converted to gc/swab from Ct values based on standard curve calculated from the concentration provided by Goyal et al. (2020).\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection="unknown",
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N, E, and ORF1lab", 
                                        unit="gc/swab",
                                        reference_event="symptom onset",)),
                participants = participants
                                        )

with open("young2020epidemiologic.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(young2020epidemiologic, outfile, default_flow_style=False, sort_keys=False)
outfile.close() 
```
