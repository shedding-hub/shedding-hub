# Extraction for CDC INHERENT Study

[CDC INHERENT](https://github.com/abtassociates/CDC_NHPHRN/blob/main/) , part of the CDC funded Nursing Home Public Health Response Network (NHPHRN), examined SARSCoV2 shedding in nursing home residents and staff. It characterized the viral shedding kinetics including proliferation, peak, and clearance using qRTPCR, antigen testing, genetic sequencing, and culture. The original dataset is published and updated on their [GitHub repository](https://github.com/abtassociates/CDC_NHPHRN/blob/main/).


First, we `import` python modules needed:

```python
import pandas as pd
import yaml
from shedding_hub import folded_str
```


```python
#load the data;
test = pd.read_excel('NHPHRN INHERENT Public Use Files and Data Dictionary.xlsx', sheet_name='Case list test results')

columns_to_keep = [
    'nh', 'outbreak_id', 'study_id', 'type', 'event', 
    'colldate', 'day', 'day_v2', 'antiviral_start', 'antigen','pcr', 'NHConductedPCRTest', 'Ct', 
    'VL', 'logVL', 'VL_status', 'variant', 'lineage'
]
test = test[columns_to_keep]

descriptive = pd.read_excel('NHPHRN INHERENT Public Use Files and Data Dictionary.xlsx', sheet_name='Case descriptive information')

columns_to_keep = [
    'study_id', 'enr_date', 'casedate', 
    'latesdat', 'symptom', 'age', 'gender', 'bio_sex', 
    'race', 'ethnicity', 'paxlovid', 'remdesivir', 
    'molnupiravir', 'monoclonal', 'steroid', 'vaccine'
]

descriptive = descriptive[columns_to_keep]

merged_df = pd.merge(test, descriptive, on='study_id', how='left')
```

Finally, the data is formatted and output as a YAML file.


```python 
participants = []

for patient_id, group in merged_df.groupby('study_id'):
    participant = {'attributes': {}, 'measurements': []}

    if pd.notna(group['age'].iloc[0]):
        participant['attributes']['age'] = int(group['age'].iloc[0])

    vaccine_val = group['vaccine'].iloc[0]
    # A numeric string representing the number of COVID-19 vaccinations received (allowed values: 0-10) or the string "COVID-19 vaccination status unknown" to indicate an unknown status.
    if pd.isna(vaccine_val) or vaccine_val == 'COVID-19 vaccination status unknown':
        participant['attributes']['vaccinated'] = 'unknown'
    else:
        vaccine_numeric = pd.to_numeric(vaccine_val, errors='coerce')
        if pd.notna(vaccine_numeric):
            participant['attributes']['vaccinated'] = vaccine_numeric > 0
            # If a numeric value is provided and it is greater than 0, the participant is considered vaccinated.
        else:
            participant['attributes']['vaccinated'] = 'unknown'

    # Process bio_sex: if missing, set as 'unknown'
    if pd.isna(group['bio_sex'].iloc[0]):
        participant['attributes']['sex'] = 'unknown'
    else:
        sex = ('female' if group['bio_sex'].iloc[0] == 'Female at birth'
               else 'male' if group['bio_sex'].iloc[0] == 'Male at birth'
               else 'unknown')
        participant['attributes']['sex'] = sex

    if pd.isna(group['ethnicity'].iloc[0]):
        participant['attributes']['ethnicity'] = 'unknown'
    else:
        ethnicity = group['ethnicity'].iloc[0]
        ethnicity_val = ('not hispanic' if ethnicity == 'Non-Hispanic'
                         else 'hispanic' if ethnicity == 'Hispanic'
                         else 'unknown')
        participant['attributes']['ethnicity'] = ethnicity_val

    if pd.isna(group['race'].iloc[0]):
        participant['attributes']['race'] = 'unknown'
    else:
        race = group['race'].iloc[0]
        race_val = ('white' if race == 'White'
                    else 'black' if race == 'Black or African American'
                    else 'other' if race in ['Other', 'American Indian or Alaska Native', 'Native Hawaiian or Other Pacific Islander']
                    else 'unknown')
        participant['attributes']['race'] = race_val

    for col in ['lineage', 'variant']:
        unique_vals = group[col].dropna().unique()
        # Some individuals have multiple non nan lineage entries. Below are all:
        # Study ID 131265-32: [nan 'XBB.1.16' 'Unassigned']
        # Study ID 131516-426: [nan 'XBB.1.16' 'Unassigned']
        # Study ID 131516-458: [nan 'XBB.1.16' 'Unassigned']
        # Study ID 131516-512: [nan 'EG.5.1' 'Unassigned']
        # We only used the valid value in those cases, so we filter out "Unassigned"
        filtered_vals = [v for v in unique_vals if v != 'Unassigned']
        if len(filtered_vals) == 0:
            participant['attributes'][col] = 'unknown'
        elif len(filtered_vals) == 1:
            participant['attributes'][col] = filtered_vals[0]
        else:
            error_message = f"Error: Multiple distinct values for {col}: {filtered_vals}"
            print(error_message)
        

    # Process measurements
    for _, row in group.iterrows():
        if row['pcr'] in ['Not tested', 'Not collected'] or pd.isna(row['pcr']):
            continue
        elif row['pcr'] == 'Negative':
            value = "negative"
        elif row['pcr'] == 'Positive':
            # Process VL_status for positive PCR results
            if row['VL_status'] == 'Not eligible for VL testing':
                continue  # Skip this row
            elif row['VL_status'] in ['Undetectable, will not be done', 'VL run, not detected']:
                value = "negative"
            elif row['VL_status'] == 'Less than 10,000 copies/mL':
                value = row["VL"] # values are reported instead of "positive"
            elif row['VL_status'] in ['Within detectable range', 'More than 1 billion copies/mL']:
                value = row["VL"]
            else:
                continue    
        else:
            continue

        measurement = {
            'analyte': 'NP_SARSCoV2_VL',
            'time': row['day'],
            'value': value
        }

        if row['pcr'] in ['Not tested', 'Not collected'] or pd.isna(row['pcr']):
            continue
        elif row['pcr'] == 'Negative':
            value_CT = "negative"
        elif row['pcr'] == 'Positive':
            # Process VL_status for positive PCR results
            if row['VL_status'] == 'Not eligible for VL testing':
                continue  # Skip this row
            elif row['VL_status'] in ['Undetectable, will not be done', 'VL run, not detected']:
                value_CT = "negative"
            elif row['VL_status'] == 'Less than 10,000 copies/mL':
                value_CT = row["Ct"] # values are reported instead of "positive"
            elif row['VL_status'] in ['Within detectable range', 'More than 1 billion copies/mL', "Less than 10,000 copies/mL"]:
                value_CT = row["Ct"]
            else:
                continue    
        else:
            continue

        measurement_CT = {
            'analyte': 'NP_SARSCoV2_CT',
            'time': row['day'],
            'value': value_CT
        }

        participant['measurements'].append(measurement)
        participant['measurements'].append(measurement_CT)

    if participant['measurements']:
        participants.append(participant)
```

```python
CDC = dict(title="Centers for Disease Control and Prevention (CDC) Nursing Home Public Health Response Network (NHPHRN)",
               url="https://github.com/abtassociates/CDC_NHPHRN/blob/main/",
               description=folded_str('The INHERENT study, part of the CDC funded Nursing Home Public Health Response Network (NHPHRN), examined SARSCoV2 shedding in nursing home residents and staff. It characterized the viral shedding kinetics including proliferation, peak, and clearance using qRTPCR, antigen testing, genetic sequencing, and culture. The original dataset is published and updated on their GitHub repository (https://github.com/abtassociates/CDC_NHPHRN/blob/main/).\n'),
               analytes=dict(NP_SARSCoV2_VL=dict(description=folded_str('Viral concentrations were quantified using RTqPCR targeting the N and S genes in nasopharyngeal swab samples with the TaqPath COVID19, Flu A, Flu B Combo Kit. One nasal swab was processed using the TaqPath Combo Kit on the QuantStudio 7 Pro, which provided results for all three analytes. For specimens with a SARSCoV2 cycle threshold (Ct) value below 30, a second analysis was performed using the CDC Influenza SARSCoV2 (Flu SC2) Multiplex Assay on the ABT 7500 Fast Dx. The reference event, "confirmation date", is defined as the date of the initial positive test.\n'),
                                        limit_of_quantification=10000,
                                        limit_of_detection="unknown",
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="S and N", 
                                        unit="gc/mL",
                                        reference_event="confirmation date",),
                            NP_SARSCoV2_CT=dict(description=folded_str('Cycle threshold (Ct) values were measured using RTqPCR targeting the N and S genes in nasopharyngeal swab samples with the TaqPath COVID19, Flu A, Flu B Combo Kit. One nasal swab was processed using the TaqPath Combo Kit on the QuantStudio 7 Pro, which provided results for all three analytes. For specimens with a SARSCoV2 cycle threshold value below 30, a second analysis was performed using the CDC Influenza SARSCoV2 (Flu SC2) Multiplex Assay on the ABT 7500 Fast Dx. The reference event, "confirmation date", is defined as the date of the initial positive test.\n'),
                                        limit_of_quantification="unknown",
                                        limit_of_detection=37,
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="S and N", 
                                        unit="cycle threshold",
                                        reference_event="confirmation date",)),
                participants = participants
                                        )

with open("cdc2024nhphrn.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(CDC, outfile, default_flow_style=False, sort_keys=False)


```