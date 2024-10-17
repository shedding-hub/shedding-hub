# Extraction for Liu et al. (2024)

[Liu et al. (2024)](https://www.medrxiv.org/content/10.1101/2024.04.22.24305845v1) measured SARS-CoV-2, pepper mild mottle virus (PMMoV), and human mitochondrial DNA (mtDNA) concentrations in longitudinal stool samples collected from 42 COVID-19 patients for up to 42 days after the first sample collection date. Demographic information (e.g., age, sex, race, ethnicity, etc.) is also included in the data. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub). Currently, the raw data does not include any symptom data, such as e.g., fever, cough, short of breath, diarrhea, headache, loss of smell, loss of taste, etc.

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
import numpy as np
```

Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
#load the data;
Liu2024 = pd.read_csv("Liu2024.csv") #Need to change the directory to load the data;

#sort by ID and date;
Liu2024 = Liu2024.sort_values(by=['subject','day_actual'])

#some data cleaning to match the schema;
Liu2024 = Liu2024.replace({"Gender": {"M": "male", "F": "female"},
                           "Race": {"White": "white", "Black or African American": "black", "Asian": "asian"},
                           "Hispanic.or.Latin.Origin": {"Yes": "hispanic", "No": "not hispanic"},
                           "inpatient": {"Yes": "inpatient", "No": "outpatient"}})
Liu2024.loc[np.isnan(Liu2024["Age"]),"Age"]="unknown"
Liu2024.loc[Liu2024["Cohort"]=="Breakthrough","Cohort"]=True
Liu2024.loc[Liu2024["Cohort"]=="Unvaccinated","Cohort"]=False
Liu2024.loc[Liu2024["dpcr_result_class_N1"]=="Negative","gc_dryg_N1"]="negative"
Liu2024.loc[Liu2024["dpcr_result_class_PMMoV"]=="Negative","gc_dryg_PMMoV"]="negative"
Liu2024.loc[Liu2024["dpcr_result_class_mtDNA"]=="Negative","gc_dryg_mtDNA"]="negative"
```

Finally, the data is formatted and output as a YAML file.

```python
participant_list = [dict(attributes=dict(age=Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"Age"],
                                         sex=Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"Gender"],
                                         race=Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"Race"],
                                         ethnicity=Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"Hispanic.or.Latin.Origin"],
                                         inpatient=Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"inpatient"],
                                         vaccinated=Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"Cohort"]),
                         measurements=[dict(analyte="stool_SARSCoV2_N1",
                                             time=int(Liu2024.loc[j,"day_actual"].item()),
                                             value=Liu2024.loc[j,"gc_dryg_N1"]) for j in Liu2024.loc[Liu2024["subject"]==i].index] +
                                       [dict(analyte="stool_PMMoV",
                                             time=int(Liu2024.loc[j,"day_actual"].item()),
                                             value=Liu2024.loc[j,"gc_dryg_PMMoV"]) for j in Liu2024.loc[Liu2024["subject"]==i].index] +
                                       [dict(analyte="stool_mtDNA",
                                             time=int(Liu2024.loc[j,"day_actual"].item()),
                                             value=Liu2024.loc[j,"gc_dryg_mtDNA"]) for j in Liu2024.loc[Liu2024["subject"]==i].index]) for i in pd.unique(Liu2024["subject"])]

liu2024 = dict(title="Longitudinal Fecal Shedding of SARS-CoV-2, Pepper Mild Mottle Virus, and Human Mitochondrial DNA in COVID-19 Patients",
               doi="10.1101/2024.04.22.24305845",
               description=folded_str('The authors measured SARS-CoV-2, pepper mild mottle virus (PMMoV), and human mitochondrial DNA (mtDNA) in longitudinal stool samples collected from 42 COVID-19 patients for up to 42 days after the first sample collection date. Abundances were quantified using Digital PCR assays targeting the N1 genes. The symptom data (e.g., fever, cough, short of breath, diarrhea, headache, loss of smell, loss of taste, etc.) is currently not included in this data.\n'),
               analytes=dict(stool_SARSCoV2_N1=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration in stool samples. The concentration were quantified in genome copies per dry weight of stool.\n"),
                                                    specimen="stool",
                                                    biomarker="SARS-CoV-2",
                                                    gene_target="N1",
                                                    limit_of_quantification=1000,
                                                    limit_of_detection="unknown",
                                                    unit="gc/dry gram",
                                                    reference_event="confirmation date"),
                             stool_PMMoV=dict(description=folded_str("PMMoV genome copy concentration in stool samples. The concentration were quantified in genome copies per dry weight of stool.\n"),
                                              specimen="stool",
                                              biomarker="PMMoV",
                                              limit_of_quantification="unknown",
                                              limit_of_detection="unknown",
                                              unit="gc/dry gram",
                                              reference_event="confirmation date"),
                             stool_mtDNA=dict(description=folded_str("mtDNA genome copy concentration in stool samples. The concentration were quantified in genome copies per dry weight of stool.\n"),
                                              specimen="stool",
                                              biomarker="mtDNA",
                                              limit_of_quantification="unknown",
                                              limit_of_detection="unknown",
                                              unit="gc/dry gram",
                                              reference_event="confirmation date")),
               participants=participant_list)

with open("liu2024longitudinal.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(liu2024, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close()
```
