import pandas as pd
import yaml
from shedding_hub import folded_str

# Load dataset
df = pd.read_excel('CombinedDataset.xlsx', sheet_name='Viral_Load')
Anthony2021 = df[df['StudyNum'] == 16].copy()
columns_to_drop = ['Estimated', 'SevMax', 'Sev1st', 'Died', 'Ctvalue', 'SevMax3']
Anthony2021 = Anthony2021.drop(columns=columns_to_drop)

# Group by participant and extract measurements
participants = []

for patient_id, group in Anthony2021.groupby('PatientID'):
    participant = {
        'attributes': {
            'age': int(group['Age'].iloc[0]),
            'sex': 'female' if group['Sex'].iloc[0] == 'F' else 'male'
        },
        'measurements': []
    }

    for _, row in group.iterrows():
        if row['value'] == 1:
            value = "negative"
        else:
            value = row['value']
        measurementN = {
            'analyte': 'swab_SARSCoV2_N',
            'time': row['Day'],
            'value': value
        }
        participant['measurements'].append(measurementN)

    participants.append(participant)

Anthony2021Tcell = dict(title="Early induction of functional SARS-CoV-2-specific T cells associates with rapid viral clearance and mild disease in COVID-19 patients",
        doi="10.1016/j.celrep.2021.108728",
        description=folded_str('Measured SARS-CoV-2 viral load in the upper respiratory tract, along with SARS-CoV-2-specific antibodies and T cell responses, at multiple time points from acute infection through convalescence or until death.'),

        analytes=dict(swab_SARSCoV2_N=dict(description=folded_str('The swab_SARSCoV2_N metric refers to the detection of SARS-CoV-2 RNA in patient swabs,using RT-PCR cycle counts as a proxy for viral quantity. Data for the nasopharyngeal swab results were obtained from the combined dataset in the supplementary materials of Challenger et al. (2022)'),

                                    limit_of_quantification="unknown",
                                    limit_of_detection="unknown",
                                    specimen="nasopharyngeal_swab",
                                    biomarker="SARS-CoV-2",
                                    gene_target="N, NP, M, ORF7, ORF8, ORF3, S, NSP7, NSP13",
                                    unit="copies/mL",
                                    reference_event="symptom onset",)),
            participants = participants
                                    )

with open("anthony2021early.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(Anthony2021Tcell, outfile, default_flow_style=False, sort_keys=False)
outfile.close()