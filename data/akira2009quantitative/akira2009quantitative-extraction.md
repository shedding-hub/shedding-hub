# Extraction for Akira et al. (2009)

[Iwakiri et al. (2009)](https://doi.org/10.1007/s00705-009-0358-0) quantified Sapovirus (SaV) RNA shedding in stool samples from cases of two outbreaks using real-time RT-PCR, revealing that SaV excretion generally declined within two weeks but can persist at high concentrations for up to four weeks in some individuals. The study also identified nucleotide substitutions in the VP1 gene during prolonged excretion, suggesting potential viral evolution.
               

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```

Raw data, extracting from Table 1: Fecal specimens collected from two outbreaks (https://link.springer.com/article/10.1007/s00705-009-0358-0/tables/1), is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/akira2009quantitative), and will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
# Read in the CSV file containing data and stored as Akira2009
Akira2009 = pd.read_csv("data.csv")
Akira2009 = Akira2009.replace({"Sex": {"M": "male", "F": "female"}})
Akira2009["type"] = "stool"

#First create a 'ID' columne to index each participant based on their information in column 'Specimens'
# Convert "1-10" -> (1,10) instead of lexicographic sorting
def sort_key(specimen):
    x, y = map(int, specimen.split("-"))  
    return (x, y)
# Sort unique specimens in numerical order
sorted_specimens = sorted(Akira2009['Specimens'].unique(), key=sort_key)
# Create a mapping dictionary for sorted specimens
specimen_ids = {specimen: i+1 for i, specimen in enumerate(sorted_specimens)}
# Map specimens to their numerical ID
Akira2009['ID'] = Akira2009['Specimens'].map(specimen_ids)


participant_list = []

for i in pd.unique(Akira2009["ID"]):
    patient_data = Akira2009[Akira2009["ID"] == i]
    age = int(patient_data['Age'].iloc[0])  # Convert to Python int
    sex = str(patient_data['Sex'].iloc[0])  # Convert to Python str

    measurements = []
    # Iterate over each row of the patient's data
    for _, row in patient_data.iterrows():
        try:
            # Convert to Python float
            value = float(row['value'])  # Use float for scientific notation

            # Add the new condition here
            if value == 0.0:
                value = 'negative'

        except ValueError:
            # Handle non-numeric values (like 'negative')
            value = 'negative'

        # Append only for the specific sample type
        if row['type'] == 'stool':
            measurements.append({
                "analyte": "stool_SaV",
                "time": int(row['days_after_onset']),
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
akira2009 = dict(title="Quantitative analysis of fecal sapovirus shedding: identification of nucleotide substitutions in the capsid protein during prolonged excretion",
               doi="10.1007/s00705-009-0358-0",
               description=folded_str("This study quantifies Sapovirus (SaV) RNA shedding in stool from two outbreak cases using real-time RT-PCR, revealing that SaV excretion generally declines within two weeks but can persist at high concentrations for up to four weeks in some individuals. The study also identifies nucleotide substitutions in the VP1 gene during prolonged excretion, suggesting potential viral evolution.\n"),
               analytes=dict(stool_SaV=dict(description=folded_str("Sapovirus RNA gene copy concentration in stool samples. The concentration were quantified in cDNA copies per gram.\n"),
                                                    specimen="stool",
                                                    biomarker="sapovirus",
                                                    gene_target="VP1",  
                                                    limit_of_quantification="unknown",
                                                    limit_of_detection=129000,
                                                    unit="gc/wet gram",
                                                    reference_event="symptom onset")),
                             
               participants=participant_list)

with open("akira2009quantitative.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")
    yaml.dump(akira2009, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close() 
```