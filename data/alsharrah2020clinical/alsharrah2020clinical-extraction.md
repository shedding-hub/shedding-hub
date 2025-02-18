# Extraction for Alsharrah et al. (2020)

[Alsharrah et al. (2020)](https://onlinelibrary.wiley.com/doi/10.1002/jmv.26684) measured SARS-CoV-2 detected by real-time reverse transcriptase PCR in paired oropharyngeal and nasopharyngeal swabs from thirty three COVID-19 patients in in Jaber Alahmad Hospital (JAH) after symptom onset. Demographic information (e.g., age, sex, etc.) is also included in the data. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/alsharrah2020clinical). 
First, we `import` python modules needed:

```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```
Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/alsharrah2020clinical), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
# Read in the CSV file containing data and store it in Alsharrah2020
Alsharrah2020 = pd.read_csv("data.csv")
Alsharrah2020  = Alsharrah2020 .replace({"Sex": {"M": "male", "F": "female"}})
#Adds a new column "type" with the value "NPS_OPS" for all rows, indicating nasopharyngeal and oropharyngeal swabs.
Alsharrah2020 ["type"] = "NPS"
#Extracts numerical IDs from the PatientID column
Alsharrah2020["ID"] = Alsharrah2020["PatientID"].str.split("-").str[1].astype(int)


participant_list = []
# Loops through each unique patient ID in the dataset
for i in pd.unique(Alsharrah2020["ID"]):
    patient_data = Alsharrah2020[Alsharrah2020["ID"] == i]
    age = int(patient_data['Age'].iloc[0])  # Convert to Python int
    sex = str(patient_data['Sex'].iloc[0])  # Convert to Python str

    measurements = []
    for _, row in patient_data.iterrows():
        try:
            # Handle missing values safely
            if pd.isna(row['value']):
                value = 'negative'
            else:
                value = float(row['value'])

            # Special condition for 1.0
            if value == 1.0:
                value = 'negative'

        except ValueError:
            value = 'negative'

        # Append measurements separately for each sample type
        if row['type'] == 'NPS':
            measurements.append({
                "analyte": "NPS_SARSCoV2", 
                "time": int(row['Day']), 
                "value": value})


    participant_dict = {"attributes": {"age": age, "sex": sex}, "measurements": measurements}
    participant_list.append(participant_dict)

```
Finally, the data is formatted and output as a YAML file.

```python
alsharrah2020 = dict(title="Clinical characteristics of pediatric SARS-CoV-2 infection and coronavirus disease 2019 (COVID-19) in Kuwait",
            doi="10.1002/jmv.26684",
            description=folded_str("This study measured SARS-CoV-2 detected by real-time reverse transcriptase PCR in paired oropharyngeal and nasopharyngeal samples from 33 COVID-19 patients in Jaber Alahmad Hospital (JAH).\n"),
            analytes=dict(NPS_SARSCoV2=dict(description=folded_str("This analyte indicates the detection of SARS-CoV-2 RNA in both nasopharyngeal and oropharyngeal swabs, but only the nasopharyngeal swab values are presented.\n"),
                          specimen="nasopharyngeal_swab",
                          biomarker="SARS-CoV-2",
                          gene_target="E, RdRP",
                          limit_of_quantification="unknown",
                          limit_of_detection="unknown",
                          unit="gc/mL",
                          reference_event="symptom onset",
        )
    ),
    participants=participant_list,
)


with open("alsharrah2020clinical.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")
    yaml.dump(alsharrah2020, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close() 
```