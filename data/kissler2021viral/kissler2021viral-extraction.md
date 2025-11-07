# Extraction for Kissler et al. (2021)

[Kissler et al. (2021)](https://doi.org/10.1371/journal.pbio.3001333) analyzed the viral dynamics of acute SARS-CoV-2 infection using repeated PCR testing of 68 individuals (90% male) during the NBA 2019–2020 season, with testing occurring almost daily for most individuals. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/kissler2021viral).

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
from shedding_hub import folded_str
```


The raw data, available on the [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/kissler2021viral), originates from [gradlab](https://github.com/gradlab/CtTrajectories/blob/main/data/ct_dat_clean.csv). It will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
Kissler2021 = pd.read_csv("ct_dat_clean.csv")
Kissler2021["type"] = "AN+OPS"          


def assign_reference_day(group):
    # first positive (Ct < 40) time using raw Date.Index
    pos_times = group.loc[group["CT.Mean"] < 40, "Date.Index"]
    if pos_times.empty:
        group["DayFromReference"] = pd.NA
    else:
        ref_time = pos_times.min()
        group["DayFromReference"] = group["Date.Index"] - ref_time 
    return group

Kissler2021 = Kissler2021.groupby("Person.ID", group_keys=False).apply(assign_reference_day)

participant_list = []

for i in pd.unique(Kissler2021["Person.ID"]):
    patient_data = Kissler2021[Kissler2021["Person.ID"] == i]
    person_id = int(patient_data["Person.ID"].iloc[0])


    measurements = []
    for _, row in patient_data.iterrows():
        
        ct_raw = row.get("CT.Mean", pd.NA)
        if pd.isna(ct_raw) or ct_raw >= 40:
            ct_value = "negative"
        else:
            ct_value = float(ct_raw)


        measurements.append({
                "analyte": "AN_OPS_SARSCoV2_ct",
                "time": int(row['DayFromReference']),
                "value": ct_value
            })
    
        viral_load_log = row.get("log10_GEperML", pd.NA)
        if pd.isna(viral_load_log) or viral_load_log==2.65760968215156:
            viral_load_value = "negative"
        else:
            viral_load_value = float(10 ** viral_load_log)
        
        measurements.append({
                "analyte": "AN_OPS_SARSCoV2_viral",
                "time": int(row["DayFromReference"]),
                "value": viral_load_value
            })

    participant_dict = {
        "attributes": {
            "person_id":person_id
            },
        "measurements": measurements}
    
    participant_list.append(participant_dict)



```

Finally, the data is formatted and output as a YAML file.


```python
kissler2021 = dict(title="Viral dynamics of acute SARS-CoV-2 infection and applications to diagnostic and public health strategies",
            doi="10.1371/journal.pbio.3001333",
            description=folded_str("This study followed 68 people (90% male) during the NBA 2019 to 2020 season with frequent RT-qPCR to map SARS-CoV-2 viral RNA trajectories and 46 had acute infections. A single low Ct (<30) strongly indicates acute infection, and a second PCR within 48 hours helps tell whether someone is early (proliferation) or late (clearance) in infection.\n"),
            analytes=dict(
                AN_OPS_SARSCoV2_ct=dict(description=folded_str("SARS-CoV-2 RNA gene copy concentration in combined anterior nares and oropharyngeal swabs. Results were reported in cycle threshold (Ct) values.\n"),
                          specimen=["anterior_nares_swab", "oropharyngeal_swab"],
                          biomarker="SARS-CoV-2",
                          gene_target="N1, N2, RdRp",
                          limit_of_quantification="unknown",
                          limit_of_detection=40,
                          unit="cycle threshold",
                          reference_event="confirmation date"),
                AN_OPS_SARSCoV2_viral=dict(description=folded_str("Estimated SARS-CoV-2 RNA genome copy concentration measured from combined anterior nares and oropharyngeal swab specimens. Viral load was derived by the study authors by transforming the corresponding Ct values using a Bayesian model implemented in Stan, and results are reported in genome copies per milliliter.\n"),
                          specimen=["anterior_nares_swab", "oropharyngeal_swab"],
                          biomarker="SARS-CoV-2",
                          gene_target="N1, N2, RdRp",
                          limit_of_quantification=10 ** 2.65760968215156,
                          limit_of_detection=10 ** 2.65760968215156,  
                          unit="gc/mL",
                          reference_event="confirmation date")          
),
    participants=participant_list,
)



with open("kissler2021viral.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")
    yaml.dump(kissler2021, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close()

```
