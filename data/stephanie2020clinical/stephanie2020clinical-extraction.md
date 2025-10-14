# Extraction for Stephanie et al. (2020)

[The COVID-19 Investigation Team (2020)](https://www.nature.com/articles/s41591-020-0877-5) measured SARS-CoV-2 detected by real-time reverse transcriptase PCR in stool, urine, serum, sputum, oropharyngeal and nasopharyngeal swabs from twelve COVID-19 patients in United States, up to 2 to 3 weeks after symptom onset. Demographic information (age and sex) and hospitalization status are also included in the data. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/covid2020clinical). 

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
import numpy as np
from shedding_hub import folded_str
```

Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/covid2020clinical), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml). Patient demographic information (age and sex) was obtained from the [Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0877-5/MediaObjects/41591_2020_877_MOESM1_ESM.pdf). For patients 1–5, who were non-hospitalized, demographic data were not available and are therefore marked as unknown in our list. When age was reported as a range, we used the median value as the representative age.

```python
Stephanie2020 = pd.read_csv("data.csv")       
Stephanie2020["Day"] = Stephanie2020["Day"].round().astype(int)


patient_info = {
    1: {'Sex': 'unknown', 'Age': 'unknown'},
    2: {'Sex': 'unknown', 'Age': 'unknown'},
    3: {'Sex': 'unknown', 'Age': 'unknown'},
    4: {'Sex': 'unknown', 'Age': 'unknown'},
    5: {'Sex': 'unknown', 'Age': 'unknown'},
    6: {'Sex': 'male', 'Age': 35},
    7: {'Sex': 'female', 'Age': 65},
    8: {'Sex': 'male', 'Age': 65},
    9: {'Sex': 'male', 'Age': 35},
    10: {'Sex': 'male', 'Age': 55},
    11: {'Sex': 'male', 'Age': 55},
    12: {'Sex': 'female', 'Age': 55}
}
info_df = pd.DataFrame.from_dict(patient_info, orient="index").reset_index()
info_df = info_df.rename(columns={'index': 'ID'})
Stephanie2020 = Stephanie2020.merge(info_df, on="ID", how="left")

participant_list = []

for i in pd.unique(Stephanie2020["ID"]):
    patient_data = Stephanie2020[Stephanie2020["ID"] == i]
    age = patient_data['Age'].iloc[0] 
    sex = patient_data['Sex'].iloc[0]
   
    measurements = []
    for _, row in patient_data.iterrows():
        try:
            value = float(row['Value'])  
            if value >= 40:
                value = 'negative'
        except ValueError:
            value = 'negative'  
      
        # Append only for the specific sample type
        if row['Specimen'] == 'Stool':
            measurements.append({
                "analyte": "Stool_SARSCoV2",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Specimen'] == 'Serum':
            measurements.append({
                "analyte": "Serum_SARSCoV2",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Specimen'] == 'Sputum':
            measurements.append({
                "analyte": "Sputum_SARSCoV2",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Specimen'] == 'OP_swab':
            measurements.append({
                "analyte": "OPS_SARSCoV2",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Specimen'] == 'NP_swab':
            measurements.append({
                "analyte": "NPS_SARSCoV2",
                "time": int(row['Day']),
                "value": value
            })

    participant_dict = {
        "attributes": {
            "age": age,
            "sex": sex
        },
        "measurements": measurements
    }
    participant_list.append(participant_dict)


```
Finally, the data is formatted and output as a YAML file.

```python
stephanie2020 = dict(title="Clinical and virologic characteristics of the first 12 patients with coronavirus disease 2019 (COVID-19) in the United States",
               doi="10.1038/s41591-020-0877-5",
               description=folded_str('This study describes the first 12 COVID-19 patients identified in the United States, tracking their clinical progression and virological characteristics. SARS-CoV-2 was detected by real-time reverse transcriptase PCR in stool, urine, serum, sputum, oropharyngeal and nasopharyngeal swabs for 2 to 3 weeks after symptom onset. Results were reported in cycle threshold (Ct) values.\n'),
               analytes=dict(Stool_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in stool samples. The concentration were quantified in gene copies per milliliter.\n"),
                                                    specimen="stool",
                                                    biomarker="SARS-CoV-2",
                                                    gene_target="N",  
                                                    limit_of_quantification="unknown",
                                                    limit_of_detection=40,
                                                    unit="gc/mL",
                                                    reference_event="symptom onset"),
                             Serum_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in serum samples. The results were reported in cycle threshold numbers.\n"),
                                              specimen="serum",
                                              biomarker="SARS-CoV-2",
                                              gene_target="N", 
                                              limit_of_quantification="unknown",
                                              limit_of_detection=40,
                                              unit="gc/mL",
                                              reference_event="symptom onset"),
                             Sputum_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in sputum samples. The concentration were quantified in gene copies per milliliter.\n"),
                                              specimen="sputum",
                                              biomarker="SARS-CoV-2",
                                              gene_target="N", 
                                              limit_of_quantification="unknown",
                                              limit_of_detection=40,
                                              unit="gc/mL",
                                              reference_event="symptom onset"), 
                             OPS_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in oropharyngeal samples. The concentration were quantified in gene copies per milliliter.\n"),
                                              specimen="oropharyngeal_swab",
                                              biomarker="SARS-CoV-2",
                                              gene_target="N", 
                                              limit_of_quantification="unknown",
                                              limit_of_detection=40,
                                              unit="gc/mL",
                                              reference_event="symptom onset"),
                             NPS_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in nasopharyngeal samples. The concentration were quantified in gene copies per milliliter.\n"),
                                              specimen="nasopharyngeal_swab",
                                              biomarker="SARS-CoV-2",
                                              gene_target="N", 
                                              limit_of_quantification="unknown",
                                              limit_of_detection=40,
                                              unit="gc/mL",
                                              reference_event="symptom onset")),               
                             
               participants=participant_list)

with open("stephanie2020clinical.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")
    yaml.dump(stephanie2020, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close() 
```

