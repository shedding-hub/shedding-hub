# Extraction for Han et al. (2020)

[Han et al. (2020)](https://doi.org/10.1093/cid/ciaa447) reports SARS-CoV-2 viral loads in different specimen types for a neonate and her mother diagnosed on 2020-03-20. The study includes nasopharyngeal, oropharyngeal, stool, plasma, saliva, and urine samples. . The viral load in the respiratory specimens gradually decreased with time and was undetectable after 17 days from the onset of symptoms. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/han2020sequential).

First, we `import` python modules needed:

```python
#import modules;
import os
import yaml
import math
import pandas as pd
import numpy as np
from shedding_hub import folded_str
```
Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/han2020sequential), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
# Read in the CSV file containing data and store it in Han2020
Han2020 = pd.read_csv("combinedataset.csv")


# Define a dictionary containing patient information (ID, Sex, Age) from information provided in Clinical Analysis in Han et al. (2020).
patient_info = {
    1: {'Sex': 'female', 'Age': 0.074},
    2: {'Sex': 'female', 'Age': 'unknown'},
}

# Define a function to map patient information into the DataFrame
def map_patient_info(df):
    df = df.copy()
    
   
    # Using .loc to ensure that the changes are done correctly on the DataFrame
    df['Sex'] = df['Patient'].map(lambda x: patient_info.get(x, {}).get('Sex'))
    df['Age'] = df['Patient'].map(lambda x: patient_info.get(x, {}).get('Age'))
    
    
    return df

Han2020 = map_patient_info(Han2020)
#Data below the dashline is undetectable, assign to 0.0. The dashed line indicates the detection limit (5.7 × 10 ** 3 copies/mL). 
value = 5.7 * 10 ** 3
Han2020.loc[(Han2020['Patient'] == 1) & (Han2020['LogValue'] < math.log10(value)), 'LogValue'] = 0.0
Han2020.loc[(Han2020['Patient'] == 2) & (Han2020['LogValue'] < math.log10(value)), 'LogValue'] = 0.0

# Initialize an empty list to store participant information
participant_list = []
# Loop through each unique 'Patient' in the DataFrame
for i in pd.unique(Han2020["Patient"]):
    patient_data = Han2020[Han2020["Patient"] == i]
    try:
        age = float(patient_data['Age'].iloc[0])  # Convert to float first and then to int to handle numeric strings
    except ValueError:
        age = 'unknown'  # If the conversion fails, keep 'unknown' 
    sex = str(patient_data['Sex'].iloc[0]) 
    
    measurements = []
    for _, row in patient_data.iterrows():
        try:
            
            value = 10 ** float(row['LogValue'])  # Use float for scientific notation
            
            if value == 1.0:
                value = 'negative'
                
        except ValueError:
            value = 'negative'  # or assign 0
        # Append only for the specific sample type
        
        if row['Type'] == 'nasopharynx':
            measurements.append({
                "analyte": "nasopharynx_E",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Type'] == 'oropharynx':
            measurements.append({
                "analyte": "oropharynx_E",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Type'] == 'plasma':
            measurements.append({
                "analyte": "plasma_E",
                "time": int(row['Day']),
                "value": value
            })    
        elif row['Type'] == 'saliva':
            measurements.append({
                "analyte": "saliva_E",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Type'] == 'urine':
            measurements.append({
                "analyte": "urine_E",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Type'] == 'sputum':
            measurements.append({
                "analyte": "sputum_E",
                "time": int(row['Day']),
                "value": value
            })

        elif row['Type'] == 'nasopharynx+oropharynx':
            measurements.append({
                "analyte": "NPSOPS_E",
                "time": int(row['Day']),
                "value": value
            })
        elif row['Type'] == 'stool':
            measurements.append({
                "analyte": "stool_E",
                "time": int(row['Day']),
                "value": value
            })

        elif row['Type'] == 'stool':
            measurements.append({
                "analyte": "stool_E",
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
han2020 = dict(
    title="Sequential Analysis of Viral Load in a Neonate and Her Mother Infected With Severe Acute Respiratory Syndrome Coronavirus 2",
    doi="10.1093/cid/ciaa447",
    description=folded_str('The study reports SARS-CoV-2 viral loads in different specimen types for a neonate and her mother diagnosed on 2020-03-20. The study includes nasopharyngeal, oropharyngeal, stool, plasma, saliva, and urine samples. Viral loads were extracted manually from Figure 1 using [WebPlotDigitizer](https://automeris.io).\n'),
    analytes=dict(nasopharynx_E=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in nasopharynx samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen="nasopharyngeal_swab",
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
            unit="gc/mL",
            reference_event="symptom onset"
        ),
        oropharynx_E=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in oropharynx samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen="oropharyngeal_swab",
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
            unit="gc/mL",
            reference_event="symptom onset"
        ),
        NPSOPS_E=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in oropharynx and nasopharynx samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
            unit="gc/mL",
            reference_event="symptom onset"
        ),
        plasma_E=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in plasma samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen="plasma",
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
            unit="gc/mL",
            reference_event="symptom onset"
        ),
        saliva_E=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in saliva samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen="saliva",
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
            unit="gc/mL",
            reference_event="symptom onset"
        ),
        urine_E=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in urine samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen="urine",
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
            unit="gc/mL",
            reference_event="symptom onset"
        ),
        sputum_E=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in sputum samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen="sputum",
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
            unit="gc/mL",
            reference_event="symptom onset"
        ),
<<<<<<< HEAD
        stool_E=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in stool samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen="stool",
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
=======
        stool_E=dict(description=folded_str("From manuscript: \"Viral RNA was detected using the PowerChek 2019-nCoV real-time polymerase chain reaction kit (Kogene Biotech, Seoul, Korea) for amplification of the E gene and the RNA-dependent RNA polymerase region of the ORF1b gene, and quantified with a standard curve that was constructed using in vitro transcribed RNA provided from the European Virus Archive.\" However, data reported in Figure 1 explicitly refer to the E gene target.\n"),
            specimen="sputum",
            biomarker="SARS-CoV-2",
            limit_of_quantification="unknown",
            limit_of_detection=5700,
            gene_target= "E",
>>>>>>> b7cf71eb5b0ba2704560bc84f07d3813018ff40f
            unit="gc/mL",
            reference_event="symptom onset"
        )
    ),
    participants=participant_list
)



with open("han2020sequential.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(han2020, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close() 
