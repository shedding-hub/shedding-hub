# Extraction for Kissler et al. (2021)

[Kissler et al. (2021)](https://doi.org/10.1101/2021.02.16.21251535) analyzed densely sampled longitudinal RT-qPCR data from 65 individuals infected with SARS-CoV-2, including 7 infected with the B.1.1.7 (Alpha) variant, using each person's lowest Ct value as the reference event time point. Inclusion criteria required each individual to have at least 5 positive PCR tests (Ct < 40), including at least one with Ct < 35. Ct value were collected using the Roche cobas SARS-CoV-2 assay and converted to estimated RNA viral concentrations via a standard curve and log-linear transformation. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/kissler2021densely).

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```

The raw data, available on the [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/kissler2021densely), originates from [skissler](https://github.com/skissler/CtTrajectories_B117/tree/main/data). It will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).


```python
# Load Ct trajectory dataset
Kissler2021 = pd.read_csv("ct_dat_refined.csv")

# Assign sample type as nasopharyngeal swab for all records
Kissler2021["type"] = "NPS"

# Round test date index to nearest integer for consistent daily time units
Kissler2021["RoundedDay"] = Kissler2021["TestDateIndex"].round().astype(int)


# Function to assign reference day per individual (first day with Ct < 40)
def assign_reference_day(group):
    # Reference day = first rounded day where CtT1 < 40 (i.e., first confirmed positive)
    ref_day = group[group["CtT1"] < 40]["RoundedDay"].min()
    # Shift all days so that reference day becomes Day 1
    group["DayFromReference"] = group["RoundedDay"] - ref_day + 1
    return group

# Apply reference day assignment to each participant
Kissler2021 = Kissler2021.groupby("PersonIDClean", group_keys=False).apply(assign_reference_day)


# Construct participant-level measurement records
participant_list = []

# Iterate over unique participant IDs
for i in pd.unique(Kissler2021["PersonIDClean"]):
    patient_data = Kissler2021[Kissler2021["PersonIDClean"] == i]

    measurements = []
    # Loop through each row of this participant's data
    for _, row in patient_data.iterrows():

        # Skip data before Day 1 (i.e., prior to first confirmed positive)
        if row["DayFromReference"] < 1:
            continue

        try:
            # Handle missing Ct values
            if pd.isna(row['CtT1']):
                value = 'negative'
            else:
                value = float(row['CtT1']) 

            # Treat Ct = 40 as negative (at limit of detection)
            if value == 40:
                value = 'negative'

        except ValueError:
            # Catch parsing errors, treat as negative
            value = 'negative'

        # Only include nasopharyngeal swab measurements
        if row['type'] == 'NPS':
            measurements.append({
                "analyte": "NPS_SARSCoV2",             # Name of biomarker
                "time": round(row['DayFromReference']), # Time since first positive (Day 1 = reference)
                "value": value                          # Ct value or 'negative'
            })

    # Store all measurements for this participant
    participant_dict = {"measurements": measurements}
    participant_list.append(participant_dict)
```


Finally, the data is formatted and output as a YAML file.


```python

# YAML dictionary structure describing the dataset
kissler2021 = dict(
    title="Densely sampled viral trajectories suggest longer duration of acute infection with B.1.1.7 variant relative to non-B.1.1.7 SARS-CoV-2",
    doi="10.1056/nejmc2102507",
    description=folded_str("This study analyzed densely sampled longitudinal RT-qPCR data from 65 individuals infected with SARS-CoV-2, including 7 infected with the B.1.1.7 (Alpha) variant, using each person's lowest Ct value as the reference event time point. Inclusion criteria required each individual to have at least 5 positive PCR tests (Ct < 40), including at least one with Ct < 35. Ct value were collected using the Roche cobas SARS-CoV-2 assay and converted to estimated RNA viral concentrations via a standard curve and log-linear transformation.\n"),
    
    # Define the biomarker being measured
    analytes=dict(
        NPS_SARSCoV2=dict(
            description=folded_str("SARS-CoV-2 RNA gene copy concentration in nasopharynx samples. The concentration was quantified in gene copies per milliliter.\n"),
            specimen="nasopharyngeal_swab",   
            biomarker="SARS-CoV-2",           
            gene_target="ORF1ab",             # Public documentation of Roche cobas shows that Target 1 = ORF1ab and Target 2 = E gene
            limit_of_quantification="unknown",
            limit_of_detection="unknown",     # Not specified numerically; Ct=40 treated as cutoff
            unit="gc/mL",                     
            reference_event="confirmation date"  # Day 1 = first PCR-positive test (Ct < 40)
        )
    ),

    # Add participant data
    participants=participant_list,
)

# Output YAML to file
with open("kissler2021densely.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")  # Optional schema annotation
    yaml.dump(kissler2021, outfile, default_style=None, default_flow_style=False, sort_keys=False)

outfile.close()
```