---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.4
  kernelspec:
    display_name: python_clean
    language: python
    name: python3
---

```python
#import modules;
import yaml
import pandas as pd
import numpy as np
from shedding_hub import folded_str
```

```python
#load the data;
df = pd.read_excel('CombinedDataset.xlsx', sheet_name='Viral_Load')
df = df[df['StudyNum']==17]
```

```python
#Drop unnecessary columns
df = df.loc[:,'Day':'Sex'].join(df.loc[:,['LOD','PatientID']])
```

```python
#Modify the values to fit the schema
df['Sex'] = df['Sex'].replace({'M':'male','F':'female'})
df['value'] = df['value'].astype('object')
df['Age'] = df['Age'].astype('object')
df['Day'] = df['Day'].astype('object')
df.loc[df['value']==1,'value'] = 'negative'
#Rename column
df.rename(columns={'PatientID': 'id'}, inplace=True)
```

```python
#Output Yaml file
#Enter the yaml writing stage
participant_list = [dict(attributes=dict(age=df.loc[df.loc[df["id"]==i].index[0],'Age'],
                            sex=df.loc[df.loc[df["id"]==i].index[0],'Sex']),
                            measurements=[dict(analyte='N1_and_N2_probes',time=j,
                                            value=df.loc[(df['Day'] == j) & (df['id'] == i),"value"].item()) for j in np.unique(df.loc[df['id'] == i,'Day'])]) for i in np.unique(df['id'])]

```

The authors stated: "a limit of detection in the range of 1000–3162 viral RNA copies per mL". We used the maximum "3162" to be conservative.

Results were considered positive if signals were detected (Ct < 40) for the RP, N1, and N2 genes. Results were classified as “not detected” if RP was detected but no signal was observed (Ct ≥ 40) from either N1 or N2. Results were classified as inconclusive if RP was detected (Ct < 40) and either N1 or N2 was detected (but not both). Results were classified as invalid if no RP was detected in the sample. Any specimens for which results were inconclusive were retested; specimens which produced inconclusive results after retesting were excluded from the analysis (n = 17 specimens from 6 participants). Ct values for amplification of both viral targets (N1 and N2 probes).

The author of the original paper stated: "For this analysis, we focused on values of the N1 probe." Thus we think the viral load in data is from N1 probe.

```python
phillip2020 = dict(title="Epidemiological Correlates of Polymerase Chain Reaction Cycle Threshold Values in the Detection of Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2)",
               doi="10.1093/cid/ciaa1469",
               description=folded_str('This study was conducted in Utah and Wisconsin between 23 March and 13 May 2020, with testing data collected during a prospective household transmission investigation of outpatient and mild coronavirus disease 2019 cases.\n'),
               analytes=dict(N1_and_N2_probes=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration calculated from evaluation of both N1 and N2 probe.\nThe authors stated: 'a limit of detection in the range of 1000-3162 viral RNA copies per mL'. We used the maximum '3162' to be conservative.\nResults were considered positive if signals were detected (Ct < 40) for the RP, N1, and N2 genes. Results were classified as 'not detected' if RP was detected but no signal was observed (Ct >= 40) from either N1 or N2. Results were classified as inconclusive if RP was detected (Ct < 40) and either N1 or N2 was detected (but not both). Results were classified as invalid if no RP was detected in the sample. Any specimens for which results were inconclusive were retested; specimens which produced inconclusive results after retesting were excluded from the analysis (n = 17 specimens from 6 participants). Ct values for amplification of both viral targets (N1 and N2 probes).\n"),
                                                    specimen="nasopharyngeal_swab",
                                                    biomarker="SARS-CoV-2",
                                                    gene_target="N1 probe",
                                                    limit_of_quantification='unknown',
                                                    limit_of_detection=3162, #1000~3162 in the article.
                                                    unit="gc/mL",
                                                    reference_event="symptom onset")),
               participants=participant_list)
```

```python
with open("salvatore2020epidemiological.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(phillip2020, outfile, default_style=None, default_flow_style=False, sort_keys=False)
```
