# Extracting Data From Vector Graphics

[Tan et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7826084/) examines how SARS-CoV-2-specific immune responses influence disease severity and viral clearance in COVID-19 patients by analyzing virological and immunological changes from symptom onset to recovery or death, with a focus on the roles of T cell responses and antibody production in controlling the infection. Data for the oral swab results were obtained from the combined dataset in the supplementary materials of [Challenger et al. (2022)](https://doi.org/10.1186/s12916-021-02220-0).

First, we import `python` modules needed:
```python
import pandas as pd
import yaml
from shedding_hub import folded_str
```

We extract and clean data from combined datasets:
```python
# Load dataset
df = pd.read_excel("CombinedDataset.xlsx", sheet_name="Viral_Load")
tan2021 = df[df["StudyNum"] == 16].copy()
columns_to_drop = ["Estimated", "SevMax", "Sev1st", "Died", "value", "SevMax3"]

tan2021 = tan2021.drop(columns=columns_to_drop)

# Group by participant and extract measurements
participants = []

for patient_id, group in tan2021.groupby("PatientID"):
    # Check if 'Age' or 'Sex' is not null
    if pd.notnull(group["Age"].iloc[0]) and pd.notnull(group["Sex"].iloc[0]):
        participant = {
            "attributes": {
                "age": int(group["Age"].iloc[0]),
                "sex": "female" if group["Sex"].iloc[0] == "F" else "male",
            },
            "measurements": [],
        }

    for _, row in group.iterrows():
        if row["Ctvalue"] == 38:
            value = "negative"
        else:
            value = row["Ctvalue"]

        measurementN = {
            "analyte": "swab_SARSCoV2_N",
            "time": row["Day"],
            "value": value,
        }
        participant["measurements"].append(measurementN)

    participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.
```python
tan2021Tcell = dict(
    title="Early induction of functional SARS-CoV-2-specific T cells associates with rapid viral clearance and mild disease in COVID-19 patients",
    doi="10.1016/j.celrep.2021.108728",
    description=folded_str(
        "Measured SARS-CoV-2 viral load in the upper respiratory tract, along with SARS-CoV-2-specific antibodies and T cell responses, at multiple time points from acute infection through convalescence or until death.\n"
    ),
    analytes=dict(
        swab_SARSCoV2_N=dict(
            description=folded_str(
                "The swab_SARSCoV2_N analyte refers to the detection of SARS-CoV-2 RNA in patient oral swabs, using RT-PCR cycle threshold (CT) values to calculate the relative quantities of SARS-CoV-2 RNA. 10^(-11) was the detection limit, which was transformed as CT value equals 36.5. Data was obtained from the combined dataset in the supplementary materials of Challenger et al. (2022)\n"
            ),
            limit_of_quantification="unknown",
            limit_of_detection=36.5,
            specimen="throat_swab", #based on the communication with the corresponding author Dr. Antonio Bertoletti at Duke-Nus Medical School (antonio@duke-nus.edu.sg).
            biomarker="SARS-CoV-2",
            unit="cycle threshold",
            reference_event="symptom onset",
        )
    ),
    participants=participants,
)

with open("tan2021early.yaml", "w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(tan2021Tcell, outfile, default_flow_style=False, sort_keys=False)
```
