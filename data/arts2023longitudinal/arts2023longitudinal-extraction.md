# Extraction for Arts et al. (2023)

[Arts et al. (2023)](https://journals.asm.org/doi/10.1128/msphere.00132-23) present longitudinal, quantitative fecal shedding data for SARS-CoV-2 RNA, pepper mild mottle virus (PMMoV) RNA, and crAss-like phage (crAssphage) DNA from 48 COVID-19 patients. Abundances were quantified using (RT)-ddPCR assays targeting the N and ORF1a genes. The data were obtained from [supplementary material](https://journals.asm.org/doi/suppl/10.1128/msphere.00132-23/suppl_file/msphere.00132-23-s0002.xlsx).

First, we `import` python modules needed:

```python
import pandas as pd
import yaml
from shedding_hub import folded_str
```

```python
#load the data;
df = pd.read_excel('msphere.00132-23-s0002.xlsx', sheet_name='Shedding Data')
Demographics = pd.read_excel('msphere.00132-23-s0002.xlsx', sheet_name='Demographics')
merged_data = pd.merge(df, Demographics, on='ID')
```


Finally, the data is formatted and output as a YAML file.

```python
participants = []

# Group by participant and extract measurements
for patient_id, group in merged_data.groupby('ID'):
    participant = {
        'attributes': {
            'age': int(group['Age'].iloc[0]),
            'sex': 'female' if group['Sex'].iloc[0] == 'Female' else 'male',
            'vaccinated': True if group['Vaccination Status'].iloc[0] == 'Yes' else False
        },
        'measurements': []
    }

    for _, row in group.iterrows():
        if pd.isna(row['N_det']):
            continue
        elif row['N_det']:
            value = row['N_conc (gc/mg-dw)']
        else:
            value = "negative"
        measurementN = {
            'analyte': 'stool_SARSCoV2_N',
            'time': row['Day'],
            'value': value
        }
        if not row['N_det']:
            measurementN['limit_of_blank'] = row['N_conc (gc/mg-dw)']
        participant['measurements'].append(measurementN)


    for _, row in group.iterrows():
        if pd.isna(row['ORF1a_det']):
            continue
        elif row['ORF1a_det']:
            value = row['ORF1a_conc (gc/gm-dw)']
        else:
            value = "negative"
        measurementP = {
            'analyte': 'stool_SARSCoV2_ORF1a',
            'time': row['Day'],
            'value': value
        }
        if not row['ORF1a_det']:
            measurementP['limit_of_blank'] = row['ORF1a_conc (gc/gm-dw)']
        participant['measurements'].append(measurementP)

    for _, row in group.iterrows():
        if pd.isna(row['PMMoV_det']):
            continue
        elif row['PMMoV_det']:
            value = row['PMMoV_conc (gc/mg-dw)']
        else:
            value = "negative"
        measurementP = {
            'analyte': 'stool_PMMoV',
            'time': row['Day'],
            'value': value
        }
        if not row['PMMoV_det']:
            if not pd.isna(row['PMMoV_conc (gc/mg-dw)']):
                measurementP['limit_of_blank'] = row['PMMoV_conc (gc/mg-dw)']
        participant['measurements'].append(measurementP)

    for _, row in group.iterrows():
        if pd.isna(row['crAss_quant']):
            continue

        if row['crAss_quant']:
            value = row['crAss_conc']
        elif row['crAss_det']:
            value = "positive"
        else:
            value = "negative"

        measurementA = {
            'analyte': 'stool_crAssphage',
            'time': row['Day'],
            'value': value
        }
        if not row['crAss_quant']:
            measurementA['limit_of_quantification'] = row['crAss_conc']
        participant['measurements'].append(measurementA)

    participants.append(participant)

```


```python

Arts2023 = dict(title="Longitudinal and quantitative fecal shedding dynamics of SARS-CoV-2, pepper mild mottle virus, and crAssphage",
               doi="10.1128/msphere.00132-23",
               description=folded_str('The authors present longitudinal, quantitative fecal shedding data for SARS-CoV-2 RNA, pepper mild mottle virus (PMMoV) RNA, and crAss-like phage (crAssphage) DNA from 48 COVID-19 patients. Abundances were quantified using (RT)-ddPCR assays targeting the N and ORF1a genes. The data were obtained from supplementary material.\n'),
               analytes=dict(stool_SARSCoV2_N=dict(description=folded_str("Concentration of RNA of the N gene quantified using (RT)-ddPCR in stool samples. The concentration was quantified in gene copies per dry weight of stool. The limit of blank (LOB), determined as the upper 95% confidence limit of the negative extraction control, ranged from 11.2 to 1,550 gc/mg-dry weight. The reported number is either the measured concentration of SARS-CoV-2 N or the LOB if the concentration wasn't detectable.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection="unknown",
                                        limit_of_blank=1550,
                                        specimen="stool",
                                        biomarker="SARS-CoV-2",
                                        gene_target="N",
                                        unit="gc/dry gram",
                                        reference_event="symptom onset",),
                            stool_SARSCoV2_ORF1a=dict(description=folded_str("Concentration of RNA of the ORF1a gene quantified using (RT)-ddPCR in stool samples. The concentration was quantified in gene copies per dry weight of stool. The limit of blank (LOB), determined as the upper 95% confidence limit of the negative extraction control, ranged from 11.2 to 1,550 gc/mg-dry weight. The reported number is either the measured concentration of SARS-CoV-2 ORF1a or the LOB if the concentration wasn't detectable.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection="unknown",
                                        limit_of_blank=1550,
                                        specimen="stool",
                                        biomarker="SARS-CoV-2",
                                        gene_target="ORF1a",
                                        unit="gc/dry gram",
                                        reference_event="symptom onset",),
                            stool_PMMoV=dict(description=folded_str("Concentration of PMMoV RNA was measured by ddPCR on the same day as the SARS-CoV-2. The concentration was quantified in gene copies per dry weight of stool. The reported number is either the measured concentration of PMMoV or the LOB if the concentration wasn't detectable.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection="unknown",
                                        limit_of_blank=400,
                                        specimen="stool",
                                        biomarker="PMMoV",
                                        unit="gc/dry gram",
                                        reference_event="symptom onset",),
                            stool_crAssphage=dict(description=folded_str("Concentration of crAssphage DNA was measured using qPCR. The reported number is either the measured concentration of crAssphage or the Limit of Quantification (LOQ) if the concentration wasn't quantifiable.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection=25,
                                        specimen="stool",
                                        biomarker="crAssphage",
                                        unit="gc/dry gram",
                                        reference_event="symptom onset",)),
                participants = participants)

with open("arts2023longitudinal.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(Arts2023, outfile, default_flow_style=False, sort_keys=False)
```
