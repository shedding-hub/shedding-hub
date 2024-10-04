# Extraction for Fajnzylber et al. (2020)

[Fajnzylber et al. (2020)](https://doi.org/10.1038/s41467-020-19057-5) quantified SARS-CoV-2 viral load from participants with a diverse range of COVID-19 disease severity, including those requiring hospitalization, outpatients with mild disease, and individuals with resolved infections. Blood was collected from hospitalized participants, non-hospitalized symptomatic individuals seeking care at a respiratory infection clinic, and participants who had recovered from known COVID-19 disease. Nasopharyngeal swabs, oropharyngeal swabs, sputum, and urine were collected from hospitalized participants. Data were obtained from the supplementary materials.

```python
# Functions to add folded blocks and literal blocks;
class folded_str(str): pass
class literal_str(str): pass

def folded_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='>')
def literal_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(folded_str, folded_str_representer)
yaml.add_representer(literal_str, literal_str_representer)

# Load dataset 
fajnzylbersup = pd.read_excel('41467_2020_19057_MOESM3_ESM.xlsx', sheet_name='Inpatient')

# Drop the unnecessary columns from the dataframe
columns_to_drop = ['Use In 1st TP Analysis', 'HTN', 'Active_cancer ', 'CLD', 'BMI']
fajnzylbersup = fajnzylbersup.drop(columns=columns_to_drop)
fajnzylbersup = fajnzylbersup.iloc[:, 0:13]

# Rename viral load columns 
fajnzylbersup.rename(columns={
    'NP_VL': 'Nasopharyngeal_SARSCoV2_N',      
    'OP_PBS_VL': 'Oropharyngeal_PBS_SARSCoV2_N',
    'OP_VTM_VL': 'Oropharyngeal_VTM_SARSCoV2_N',  
    'Sputum_VL': 'Sputum_SARSCoV2_N',            
    'Plasma_VL': 'Plasma_SARSCoV2_N',             
    'Urine_VL': 'Urine_SARSCoV2_N'               
}, inplace=True)

# Convert log10 values to actual values 
columns_to_transform = [
    'Nasopharyngeal_SARSCoV2_N',
    'Oropharyngeal_PBS_SARSCoV2_N',
    'Oropharyngeal_VTM_SARSCoV2_N',
    'Sputum_SARSCoV2_N',
    'Plasma_SARSCoV2_N',
    'Urine_SARSCoV2_N'
]

fajnzylbersup[columns_to_transform] = fajnzylbersup[columns_to_transform].apply(lambda x: 10 ** x)

```

```python
participants = []

# Define mappings for race and ethnicity
def map_race_ethnicity(value):
    race_mapping = {
        'W': 'white',
        'B': 'black',
        'O': 'other',
        'U': 'unknown',
        'H': 'unknown',  # 'H' stands for Hispanic, so race is unknown
        'WH': 'white',
        'BH': 'black'
    }
    ethnicity_mapping = {
        'H': 'hispanic',
        'W': 'not hispanic',
        'B': 'not hispanic',
        'O': 'not hispanic',
        'U': 'unknown',
        'WH': 'hispanic',
        'BH': 'hispanic'

    }
    race = race_mapping.get(value, 'unknown')  # Default to unknown if not in mapping
    ethnicity = ethnicity_mapping.get(value, 'unknown')  # Default to unknown if not in mapping
    return race, ethnicity

# Group by participant and extract measurements
for patient_id, group in fajnzylbersup.groupby('PID'):
    race_value = group['Race_Ethnicity'].iloc[0]
    race, ethnicity = map_race_ethnicity(race_value)
    
    participant = {
        'attributes': {
            'age': int(group['Age'].iloc[0]),
            'sex': 'female' if group['Sex'].iloc[0] == 'F' else 'male',
            'race': race,  # Add race attribute
            'ethnicity': ethnicity  # Add ethnicity attribute
        },
        'measurements': []
    }

    for _, row in group.iterrows():
        for column in group.columns:
            # Skip columns that are not analytes or should not be used
            if column in ['PID', 'Age', 'Sex', 'Race_Ethnicity', 'Sx_Onset_to_SC', 'Diabetes ', 'O2_on_date', 'Pregnancy']:
                continue
            
            value = row[column]
            if pd.notna(value):  # Only process non-NA values
                measurementN = {
                    'analyte': column,  # Use the column name as analyte
                    'time': int(row['Sx_Onset_to_SC']) if pd.notna(row['Sx_Onset_to_SC']) else "unknown",
                    'value': "negative" if value == 10 else value
                }
                participant['measurements'].append(measurementN)
    
    participants.append(participant)

```

The data is formatted and output as a YAML file.

```python
fajnzylber2020sars = dict(title="SARS-CoV-2 viral load is associated with increased disease severity and mortality",
               doi="10.1038/s41467-020-19057-5",
               description=folded_str('The paper quantified SARS-CoV-2 viral load from participants with a diverse range of COVID-19 disease severity, including those requiring hospitalization, outpatients with mild disease, and individuals with resolved infections. Blood was collected from hospitalized participants, non-hospitalized symptomatic individuals seeking care at a respiratory infection clinic, and participants who had recovered from known COVID-19 disease. Nasopharyngeal swabs, oropharyngeal swabs, sputum, and urine were collected from hospitalized participants. Data were obtained from the supplementary materials.\n'),
               analytes=dict(Nasopharyngeal_SARSCoV2_N=dict(description=folded_str("Nasopharyngeal swabs were collected in 3 mL of phosphate-buffered saline. Viral concentrations were quantified using RT-qPCR targeting the N gene in the nasopharyngeal swab samples.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection=40,
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N1", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",),
                             Oropharyngeal_PBS_SARSCoV2_N=dict(description=folded_str("Oropharyngeal swabs were collected in 3 mL of phosphate buffered saline. Viral concentrations were quantified using RT-qPCR targeting the N gene in oropharyngeal swab samples.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection=40,
                                        specimen="oropharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N1", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",),
                             Oropharyngeal_VTM_SARSCoV2_N=dict(description=folded_str("Oropharyngeal swabs were collected in viral transport medium. Viral concentrations were quantified using RT-qPCR targeting the N gene in oropharyngeal swab samples.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection=40,
                                        specimen="oropharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N1", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",),                      
                             Sputum_SARSCoV2_N=dict(description=folded_str("Sputum samples were collected and viral concentrations were quantified using RT-qPCR targeting the N gene.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection=40,
                                        specimen="sputum", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N1", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",),           
                             Plasma_SARSCoV2_N=dict(description=folded_str("Plasma samples were collected and viral concentrations were quantified using RT-qPCR targeting the N gene.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection=40,
                                        specimen="plasma", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N1", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",),
                             Urine_SARSCoV2_N=dict(description=folded_str("Urine samples were collected and viral concentrations were quantified using RT-qPCR targeting the N gene.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection=40,
                                        specimen="urine", 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N1", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",)                      
                                        ),
                participants = participants
                                        )

with open("fajnzylber2020sars.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(fajnzylber2020sars, outfile, default_flow_style=False, sort_keys=False)
outfile.close() 

```