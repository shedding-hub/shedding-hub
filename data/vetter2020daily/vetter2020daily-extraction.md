# Extraction for Vetter et al. (2020)

[Vetter et al. (2020)](https://journals.asm.org/doi/10.1128/msphere.00827-20) measured SARS-CoV-2 detected by real-time reverse transcriptase PCR in both oropharyngeal (OPS) and nasopharyngeal (NPS) swabs from five COVID-19 patients in Geneva, Switzerland, up to days 7 and 19 after symptom onset. Demographic information (e.g., age, sex, etc.) is also included in the data. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/vetter2020daily). Currently, this dataset does not include cellular and humoral SARS-CoV-2-specific adaptive responses and serology data.

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
import numpy as np
from shedding_hub import folded_str
```

Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/vetter2020daily), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
# Read in the CSV file containing data and store it in df_1
Vetter2020 = pd.read_csv("msphere.00827-20-st002.csv")

# Define a dictionary containing patient information (ID, Sex, Age) from [Table 1](https://journals.asm.org/doi/10.1128/msphere.00827-20#tab1) in Vetter et al. (2020).
patient_info = {
    'P1': {'PatientID': 1, 'Sex': 'M', 'Age': 28},
    'P2': {'PatientID': 2, 'Sex': 'M', 'Age': 30},
    'P3': {'PatientID': 3, 'Sex': 'M', 'Age': 55},
    'P4': {'PatientID': 4, 'Sex': 'M', 'Age': 24},
    'P5': {'PatientID': 5, 'Sex': 'M', 'Age': 66}
}

# Define a function to map patient information into the DataFrame
def map_patient_info(df):
    df = df.copy() # Create a copy of the DataFrame to avoid modifying the original data

     # Map 'PatientID', 'Sex', and 'Age' based on 'patient number' column using the patient_info dictionary
    df.loc[:, 'PatientID'] = df['patient number'].map(lambda x: patient_info.get(x, {}).get('PatientID'))
    df.loc[:, 'Sex'] = df['patient number'].map(lambda x: patient_info.get(x, {}).get('Sex'))
    df.loc[:, 'Age'] = df['patient number'].map(lambda x: patient_info.get(x, {}).get('Age'))

    return df
# Apply the mapping function to df_1 and save the updated DataFrame to a CSV file
Vetter2020 = map_patient_info(Vetter2020)
# Sort the DataFrame by 'PatientID' and 'dpo' (days post onset)
Vetter2020 = Vetter2020.sort_values(by=['PatientID','dpo'])
# Replace values to match the schema
Vetter2020 = Vetter2020.replace({"Sex": {"M": "male", "F": "female"},
                                 "Virus isolation": {"yes": "virus isolated", "No": "virus unisolated", "n.d.": "not available"},
                                 "Log copies/ml": {"neg": "negative"}})

# Initialize an empty list to store participant information
participant_list = []
# Loop through each unique 'PatientID' in the DataFrame
for i in pd.unique(Vetter2020["PatientID"]):
    patient_data = Vetter2020[Vetter2020["PatientID"] == i]
    age = int(patient_data['Age'].iloc[0])  # Convert to Python int
    sex = str(patient_data['Sex'].iloc[0])  # Convert to Python str

    measurements = []
    # Iterate over each row of the patient's data
    for _, row in patient_data.iterrows():
        try:
            # Convert 'Log copies/ml' to float
            value = float(row['Log copies/ml'])  # Use float for scientific notation

            # Add the new condition here
            if value == 0.0:
                value = 'negative'

        except ValueError:
            # Handle non-numeric values (like 'negative')
            value = 'negative'

        # Append only for the specific sample type
        if row['Sample'] == 'FNP':
            measurements.append({
                "analyte": "NPS_SARSCoV2",
                "time": int(row['dpo']),
                "value": value
            })
        elif row['Sample'] == 'OPS':
            measurements.append({
                "analyte": "OPS_SARSCoV2",
                "time": int(row['dpo']),
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
vetter2020 = dict(title="Daily Viral Kinetics and Innate and Adaptive Immune Response Assessment in COVID-19: a Case Series",
               doi="10.1128/msphere.00827-20",
               description=folded_str('The author measured SARS-CoV-2 detected by real-time reverse transcriptase PCR in both oropharyngeal (OPS) and nasopharyngeal (NPS) swabs from five COVID-19 patients in Geneva, Switzerland, up to days 7 and 19 after symptom onset. Cellular and humoral SARS-CoV-2-specific adaptive responses and serology data are currently not included in this dataset.\n'),
               analytes=dict(NPS_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentrations in nasopharyngeal samples. The concentrations were quantified in gene copies per mL.\n"),
                                                    specimen="nasopharyngeal_swab",
                                                    biomarker="SARS-CoV-2",
                                                    gene_target="unknown",
                                                    limit_of_quantification="unknown",
                                                    limit_of_detection="unknown",
                                                    unit="gc/mL",
                                                    reference_event="symptom onset"),
                             OPS_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentrations in oropharyngeal samples. The concentrations were quantified in gene copies per mL.\n"),
                                              specimen="oropharyngeal_swab",
                                              biomarker="SARS-CoV-2",
                                              gene_target="unknown",
                                              limit_of_quantification="unknown",
                                              limit_of_detection="unknown",
                                              unit="gc/mL",
                                             reference_event="symptom onset")),

               participants=participant_list)

with open("vetter2020daily.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")
    yaml.dump(vetter2020, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close()
```
