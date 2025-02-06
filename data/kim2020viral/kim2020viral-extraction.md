# Extraction for Kim et al. (2020)

[Kim et al. (2020)](https://doi.org/10.3346/jkms.2020.35.e86) presents viral load kinetics for the first two confirmed COVID-19 patients with mild to moderate illness in Korea. Swabs, sputum, serum, plasma, urine, and stool samples were collected throughout the course of the illness. Cycle threshold (Ct) values were quantified using rRT-PCR targeting the RdRp and E genes. Data were obtained from Supplementary Table 3, which contains the estimated number of viral copies in respiratory specimens. The swab RdRp estimated viral load was cross-validated with the combined dataset in the supplementary materials of Challenger et al., BMC Medicine (2022) 20:25 (https://doi.org/10.1186/s12916-021-02220-0). The sputum information from Supplementary Table 1, 2, 3 have been formatted into kim2020viral.xlsx file with an additional column to label the specimen type.

First, we `import` python modules needed:

```python
import pandas as pd
import yaml
from shedding_hub import folded_str
```

```python
KimJY = pd.read_excel('kim2020viral.xlsx')
swab_RdRp = KimJY[KimJY['Type '] == 'swab RdRp estimated viral copy/mL']
sputum_RdRp = KimJY[KimJY['Type '] == 'sputum RdRp']
```


```python
participants = []

# Group by participant and extract measurements
for patient_id, group in KimJY.groupby('PatientID'):
    participant = {
        'attributes': {
            'age': int(group['Age'].iloc[0]),
            'sex': 'female' if group['Sex'].iloc[0] == 'F' else 'male'
        },
        'measurements': []
    }

    for _, row in group.iterrows():
        if row['Ctvalue'] == 'ud':
            value = "negative"
        else:
            value = row['Value'] if pd.notna(row['Value']) else row['Ctvalue']
        measurementN = {
            'analyte': row['Type '],
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurementN)

    participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.

```python
KimJY2020 = dict(title="Viral Load Kinetics of SARS-CoV-2 Infection in First Two Patients in Korea",
               doi="10.3346/jkms.2020.35.e86",
               description=folded_str('The authors present viral load kinetics for the first two confirmed COVID-19 patients with mild to moderate illness in Korea. Swabs, sputum, serum, plasma, urine, and stool samples were collected throughout the illness. Cycle threshold (Ct) values were quantified using rRT-PCR targeting the RdRp and E genes. Data were obtained from the supplementary material.\n'),
               analytes=dict(sputum_SARSCoV2_RdRp_VL=dict(description=folded_str("Cycle threshold (Ct) values were quantified using rRT-PCR targeting the RdRp gene in sputum samples. Ct values were then converted into estimated viral copies per mL, with a detection limit of 2,690 copies/mL. The RNA copy number was calculated by the authors using a standard curve derived from Ct values of plasmid DNA, as noted in the supplementary material.\n"),
                                        limit_of_quantification=2690,
                                        limit_of_detection="unknown",
                                        specimen="sputum",
                                        biomarker="SARS-CoV-2",
                                        gene_target="RdRp",
                                        unit="gc/mL",
                                        reference_event="symptom onset",),
                            swab_SARSCoV2_RdRp_VL=dict(description=folded_str("Cycle threshold (Ct) values were quantified using rRT-PCR targeting the RdRp gene in swab samples. Ct values were then converted into estimated viral copies per mL, with a detection limit of 2,690 copies/mL. The RNA copy number was calculated by the authors using a standard curve derived from Ct values of plasmid DNA, as noted in the supplementary material.\n"),
                                        limit_of_quantification=2690,
                                        limit_of_detection="unknown",
                                        specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
                                        biomarker="SARS-CoV-2",
                                        gene_target="RdRp",
                                        unit="gc/mL",
                                        reference_event="symptom onset",),
                            swab_SARSCoV2_E_Ct=dict(description=folded_str("Cycle threshold (Ct) values were quantified using rRT-PCR targeting the E gene in swab samples.\n"),
                                        limit_of_quantification=35,
                                        limit_of_detection="unknown",
                                        specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
                                        biomarker="SARS-CoV-2",
                                        gene_target="E",
                                        unit="cycle threshold",
                                        reference_event="symptom onset",),
                            sputum_SARSCoV2_E_Ct=dict(description=folded_str("Cycle threshold (Ct) values were quantified using rRT-PCR targeting the E gene in sputum samples.\n"),
                                        limit_of_quantification=35,
                                        limit_of_detection="unknown",
                                        specimen="sputum",
                                        biomarker="SARS-CoV-2",
                                        gene_target="E",
                                        unit="cycle threshold",
                                        reference_event="symptom onset",),
                            stool_SARSCoV2_E_Ct=dict(description=folded_str("Cycle threshold (Ct) values were quantified using rRT-PCR targeting the E gene in stool samples.\n"),
                                        limit_of_quantification=35,
                                        limit_of_detection="unknown",
                                        specimen="stool",
                                        biomarker="SARS-CoV-2",
                                        gene_target="E",
                                        unit="cycle threshold",
                                        reference_event="symptom onset",),
                                        ),
                            stool_SARSCoV2_RdRp_Ct=dict(description=folded_str("Cycle threshold (Ct) values were quantified using rRT-PCR targeting the RdRp gene in stool samples.\n"),
                                        limit_of_quantification=35,
                                        limit_of_detection="unknown",
                                        specimen="stool",
                                        biomarker="SARS-CoV-2",
                                        gene_target="RdRp",
                                        unit="cycle threshold",
                                        reference_event="symptom onset",),
                participants = participants
                                        )

with open("kim2020viral.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(KimJY2020, outfile, default_flow_style=False, sort_keys=False)
outfile.close()
```
