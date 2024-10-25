# Extraction for Ke et al. (2022)

[Ke et al. (2022)](https://www.nature.com/articles/s41564-022-01105-z) studied the dynamics of infectious virus and viral RNA shedding for SARS-CoV-2 during acute infection through daily longitudinal sampling of 60 individuals for up to 14 days. Nasal swab and saliva samples were collected daily and tested for University of Illinois at Urbana-Champaign faculty, staff, and students (during the fall of 2020 and spring of 2021) who reported a negative RT–qPCR test result in the past 7 days and were either (1) within 24 h of a positive RT–qPCR result or (2) within 5 days of exposure to someone with a confirmed positive RT–qPCR result. Nasal swab samples were also cultured but the results were not included in our dataset. The raw data was obtained from [here](https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-022-01105-z/MediaObjects/41564_2022_1105_MOESM4_ESM.xlsx).

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
import numpy as np

#functions to add folded blocks and literal blocks;
class folded_str(str): pass
class literal_str(str): pass

def folded_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='>')
def literal_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(folded_str, folded_str_representer)
yaml.add_representer(literal_str, literal_str_representer)
```

Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/ke2022daily), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml). The concentration of SARS-CoV-2 RNA in nasal swab samples and saliva samples were calculated based on the calibration curves provided in [Ke et al. (2022)](https://www.nature.com/articles/s41564-022-01105-z).

```python
#load the data;
Ke2022 = pd.read_excel('41564_2022_1105_MOESM4_ESM.xlsx',sheet_name='data_samples')

#sort by ID and date;
Ke2022 = Ke2022.sort_values(by=['Ind','Time'])

#calculate the concentrations based on nasal CN values using the formula, log10(V)=11.35-0.25CN, in the subsection 
#"viral genome load calibration: nasal samples" in Ke et al. (2022).
def cal_conc_CN(CN_value):
    if CN_value!=48:
        return 10**(11.35-0.25*CN_value)
    elif CN_value==48:
        return 'negative'

Ke2022['nasal_conc']=Ke2022['Nasal_CN'].map(cal_conc_CN)

#calculate the concentrations based on saliva Ct values using the formula, log10(V)=14.24-0.28Ct in the subsection 
#"viral genome load calibration: saliva samples" in Ke et al. (2022).
def cal_conc_Ct(Ct_value):
    if Ct_value!=47:
        return 10**(14.24-0.28*Ct_value)
    elif Ct_value==47:
        return 'negative'

Ke2022['saliva_conc']=Ke2022['Saliva_Ct'].map(cal_conc_Ct)
```

We adjusted the time as the day after the confirmation (first antigen test positive).

```python
#identify the first antigen positive day as the confirmation day for each subject;
Ke2022['cum_cnt'] = Ke2022.groupby(['Ind','Antigen']).cumcount()
Ke2022ConfirmTime = Ke2022.loc[(Ke2022['Antigen']=='Pos') & (Ke2022['cum_cnt']==0),['Ind','Time']]
Ke2022ConfirmTime = Ke2022ConfirmTime.rename(columns={'Time': 'Confirm_Time'})
#two subjects do not have any antigen positive. we used the first saliva positive day as confirmation day for those two subjects: 449614,'451146 *'.
Ke2022ConfirmTime = pd.concat([Ke2022ConfirmTime,pd.DataFrame(dict(Ind=[449614,'451146 *'], Confirm_Time=[-2, -4]))])
#change the 'Time' as the day after confirmation;
Ke2022=Ke2022.merge(Ke2022ConfirmTime, how='left', on='Ind')
Ke2022['Time']=Ke2022['Time']-Ke2022['Confirm_Time']
```

Finally, the data is formatted and output as a YAML file.

```python
participant_list = [dict(attributes=dict(age=int(Ke2022.loc[Ke2022.loc[Ke2022["Ind"]==i].index[0],"Age"]),
                                         lineage=Ke2022.loc[Ke2022.loc[Ke2022["Ind"]==i].index[0],"Lineage"]),
                         measurements=[dict(analyte="nasal_SARSCoV2",
                                             time=int(Ke2022.loc[j,"Time"].item()),
                                             value=Ke2022.loc[j,"nasal_conc"]) for j in Ke2022.loc[(Ke2022["Ind"]==i) & (pd.notna(Ke2022['nasal_conc']))].index] +
                                       [dict(analyte="saliva_SARSCoV2",
                                             time=int(Ke2022.loc[j,"Time"].item()),
                                             value=Ke2022.loc[j,"saliva_conc"]) for j in Ke2022.loc[(Ke2022["Ind"]==i) & (pd.notna(Ke2022['saliva_conc']))].index]) for i in pd.unique(Ke2022["Ind"])]

ke2022 = dict(title="Daily Longitudinal Sampling of SARS-CoV-2 Infection Reveals Substantial Heterogeneity in Infectiousness",
               doi="10.1038/s41564-022-01105-z",
               description=folded_str("The authors studied the dynamics of infectious virus and viral RNA shedding for SARS-CoV-2 during acute infection through daily longitudinal sampling of 60 individuals for up to 14 days. Nasal swab and saliva samples were collected daily and tested for University of Illinois at Urbana-Champaign faculty, staff, and students (during the fall of 2020 and spring of 2021) who reported a negative RT-qPCR test result in the past 7 days and were either within 24 h of a positive RT-qPCR result or within 5 days of exposure to someone with a confirmed positive RT-qPCR result.\n"),
               analytes=dict(nasal_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration in mid-turbinate nasal swab (nasopharyngeal swab) samples. Note that the unit of these measurements is per mL: this is because nasal swab samples were each collected in 3 mL of VTM. The calibration curve for nasal samples was in the Supplemental Table S10 in Ke et al. (2022).\n"),
                                                    specimen="nasopharyngeal_swab",
                                                    biomarker="SARS-CoV-2",
                                                    gene_target="N1 and N2",
                                                    limit_of_quantification="unknown",
                                                    limit_of_detection="unknown",
                                                    unit="gc/mL",
                                                    reference_event="confirmation date"),
                             saliva_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration in saliva samples. The study was not able to measure the calibration curve using saliva samples taken from participants. Instead, the authors used data from calibration experiments in which saliva samples obtained from healthy donors were spiked with SARS-CoV-2 genomic RNA. The calibration curve for saliva samples was in the Supplemental Table S11 in Ke et al. (2022).\n"),
                                              specimen="saliva",
                                              biomarker="SARS-CoV-2",
                                              gene_target="N",
                                              limit_of_quantification="unknown",
                                              limit_of_detection="unknown",
                                              unit="gc/mL",
                                              reference_event="confirmation date")),
               participants=participant_list)

with open("ke2022daily.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(ke2022, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close()
```

```python
Ke2022[0:20]
```
