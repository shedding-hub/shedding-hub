# Extraction for Kissler et al. (2021)

[Kissler et al. (2021)]( https://doi.org/10.1371/journal.pbio.3001333) analyzed densely sampled longitudinal RT-qPCR data from 65 individuals infected with SARS-CoV-2, including 7 infected with the B.1.1.7 (Alpha) variant, using each person's lowest Ct value as the reference event time point. Inclusion criteria required each individual to have at least 5 positive PCR tests (Ct < 40), including at least one with Ct < 35. Ct value were collected using the Roche cobas SARS-CoV-2 assay and converted to estimated RNA viral concentrations via a standard curve and log-linear transformation. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/kissler2021densely).

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```

The raw data, available on the [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/kissler2021densely), originates from [skissler](https://github.com/skissler/CtTrajectories_B117/tree/main/data). It will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).


```python
Kissler2021 = pd.read_csv("ct_dat_refined.csv")
Kissler2021["type"] = "AN+OPS"    

# Round test date index to nearest integer for consistent daily time units
Kissler2021["RoundedDay"] = Kissler2021["TestDateIndex"].round().astype(int)

# Function to assign reference day per individual (first day with Ct < 40)

def assign_reference_day(group):
    # first positive (Ct < 40) time using raw TestDateIndex
    pos_times = group.loc[group["CtT1"] < 40, "TestDateIndex"]
    if pos_times.empty:
        group["DayFromReference"] = pd.NA
    else:
        ref_time = pos_times.min()
        group["DayFromReference"] = group["TestDateIndex"] - ref_time + 1
    return group

Kissler2021 = Kissler2021.groupby("PersonIDClean", group_keys=False).apply(assign_reference_day)
Kissler2021["DayFromReferenceRounded"] = (
    Kissler2021["DayFromReference"].round().astype("Int64")  
)


def extract_lineage(patient_data):
    if "B117Status" in patient_data.columns:
        status = patient_data["B117Status"].astype(str)
        if not status.empty:
            if status.iloc[0] == "Yes":
                return "B.1.1.7"
            elif status.iloc[0] == "No":
                return "non-B.1.1.7"
    return "unknown"


participant_list = []
# Iterate over unique participant IDs
for i in pd.unique(Kissler2021["PersonIDClean"]):
    patient_data = Kissler2021[Kissler2021["PersonIDClean"] == i]
    
    lineage_val = extract_lineage(patient_data)

    measurements = []
    for _, row in patient_data.iterrows():
        if pd.isna(row["DayFromReferenceRounded"]):
            continue
        time_idx = int(row["DayFromReferenceRounded"])
        # Treat Ct = 40 as negative (at limit of detection)
        ct_raw = row.get("CtT1", pd.NA)
        if pd.isna(ct_raw) or ct_raw == 40:
            value = "negative"
        else:
            value = float(ct_raw)

        if row["type"] == "AN+OPS":
            measurements.append({
                "analyte": "AN_OPS_SARSCoV2",
                "time": time_idx,
                "value": value
            })
    

    participant_dict = {
        "attributes": {
            "lineage":lineage_val
            },
        "measurements": measurements}
    
    participant_list.append(participant_dict)

```


Finally, the data is formatted and output as a YAML file.


```python

# YAML dictionary structure describing the dataset
kissler2021 = dict(title="Densely sampled viral trajectories suggest longer duration of acute infection with B.1.1.7 variant relative to non-B.1.1.7 SARS-CoV-2",
            doi="10.1056/nejmc2102507",
            description=folded_str("This study analyzed densely sampled longitudinal RT-qPCR data from 65 individuals infected with SARS-CoV-2, including 7 infected with the B.1.1.7 (Alpha) variant, using each person's lowest Ct value as the reference event time point. Inclusion criteria required each individual to have at least 5 positive PCR tests (Ct < 40), including at least one with Ct < 35. Ct value were collected using the Roche cobas SARS-CoV-2 assay and converted to estimated RNA viral concentrations via a standard curve and log-linear transformation.\n"),
            analytes=dict(AN_OPS_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in combined anterior nares and oropharyngeal swabs. The concentration was quantified in gene copies per milliliter.\n"),
                          specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
                          biomarker="AN_OPS_SARSCoV2",
                          gene_target=["N1, N2, RdRp"],
                          limit_of_quantification="unknown",
                          limit_of_detection=40,
                          unit="gc/mL",
                          reference_event="confirmation date",
        )
    ),
    participants=participant_list,
)


with open("kissler2021viral.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")
    yaml.dump(kissler2021, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close()
```