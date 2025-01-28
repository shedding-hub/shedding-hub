# Extraction for Tim et al. (2016)

[Tim et al. (2016)](https://datadryad.org/stash/dataset/doi:10.5061/dryad.1p3kn) measured Influenza A virus within households during December 2008, focusing on viral loads in index cases and their household contacts to identify secondary infections. The study includes oropharynx and nasopharynx samples. The raw data is stored at [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/tim2016individual).

First, we `import` python modules needed:

```python
import os
import yaml
import pandas as pd
import numpy as np
from shedding_hub import folded_str
```

Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/tim2016individual), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
# Read in the CSV file containing data and stored as RawData
RawData = pd.read_csv("data.csv", header=None)
columns_extract = [0, 1, 3, 4, 7, 8, 9] + list(range(15, 28)) + [31]

# Read only the required columns from the dataset
ExtractedData = RawData.iloc[:, columns_extract]

# Rename the columns 
column_names = [
    "HouseholdID", "MemberID", "PCR_confirmed_infection", "DateofSymptomOnset",
    "Age", "Sex", "Vaccination", "Day0", "Day1", "Day2", "Day3", "Day4", 
    "Day5", "Day6", "Day7", "Day8", "Day9", "Day10", "Day11", "Day12", 
    "Influenza_virus_subtype"
]

ExtractedData.columns = column_names
ExtractedData.to_csv("ExtractedData.csv", index=False) 

Tim2016 = pd.read_csv("ExtractedData.csv")
Tim2016 = Tim2016.replace({"Sex": {1: "male", 0: "female"}})
# Replace numerical values in 'Influenza_virus_subtype' with subtype names displayed in README.txt
Tim2016 = Tim2016.replace({"Influenza_virus_subtype": {3: "Pandemic A(H1N1)", 2: "Seasonal A(H3N2)", 1: "Seasonal A(H1N1)", 0: "nsubtyple influenz"}})
# Convert the 'Vaccination' column to boolean values (1 -> True, 0 -> False)
Tim2016["Vaccination"] = Tim2016["Vaccination"].replace({1: True, 0: False}).astype(bool)
#print(Tim2016["Vaccination"].dtype) 
# Add a new column 'Type' to label the sample type for all rows
Tim2016["Type"] = "NPS+OPS"

#Separate index cases and non-index cases into two DataFrames for further processing

#1. Filter for index cases (MemberID == 0 and PCR_confirmed_infection == 1)
index_Tim2016 = Tim2016[(Tim2016['MemberID'] == 0) & (Tim2016['PCR_confirmed_infection'] == 1)] #keep both house id and member id

#2. Filter for non-index cases (MemberID != 0 and PCR_confirmed_infection == 1)
unindex_Tim2016= Tim2016[(Tim2016['MemberID'] != 0) & (Tim2016['PCR_confirmed_infection'] == 1)]

#1. Create index participant list
participant_list_index = []

# Loop through each unique HouseholdID in the filtered dataset
for household_id in pd.unique(index_Tim2016["HouseholdID"]):
    patient_data = index_Tim2016[index_Tim2016["HouseholdID"] == household_id]
    
    # Extract age and sex
    age = int(patient_data['Age'].iloc[0])
    sex = str(patient_data['Sex'].iloc[0])
    vaccination = patient_data['Vaccination'].iloc[0] 
    virus_subtype = str(patient_data['Influenza_virus_subtype'].iloc[0])


    
    # Initialize the measurements list
    measurements = []
    
    # Iterate through each day's value (Day0 to Day12)
    #To adjust the reference date to symptom onset for secondary cases in the households: day - (secondary symptom onset - index symptom onset)
    for day in range(0, 13): #day = day(0-12) - (secondary symptom onset - index symptom onset)
        day_column = f'Day{day}'
        
        for _, row in patient_data.iterrows():
            try:
                value = 10 ** float(row[day_column])
                
                if value == 0.1:
                    value = 'missing'
                    
            except ValueError:
                value = 'missing'
                
            # Remove all the measurements with missing values, samples were not collected on those days
            if row['Type'] == 'NPS+OPS' and value != 'missing':
                measurements.append({
                    "analyte": "NPSOPS",
                    "time": day,
                    "value": value
                })

    # Create the index participant dictionary
    participant_dict_index = {
        "attributes": {
            "householdid": household_id,
            "age": age,
            "sex": sex,
            "vaccination": vaccination,
            "virus_subtype": virus_subtype
        },
        "measurements": measurements
    }
    
    # Append to the participant list
    participant_list_index.append(participant_dict_index)

    # for participant in participant_list_index[:10]:
    #     print(participant)


#2. Create index participant list
participant_list_unindex = []
# Process non-index cases similarly, adjusting days based on symptom onset
for household_id in pd.unique(unindex_Tim2016["HouseholdID"]):  
    patient_data = unindex_Tim2016[unindex_Tim2016["HouseholdID"] == household_id]
    index_case_data = index_Tim2016[index_Tim2016["HouseholdID"] == household_id]

    if index_case_data.empty:
        continue  # Skip households without an index case

    index_symptom_onset = index_case_data['DateofSymptomOnset'].iloc[0]

    for _, participant_row in patient_data.iterrows():
        # Extract participant attributes
        age = int(participant_row['Age'])
        sex = str(participant_row['Sex'])
        vaccinated = str(participant_row['Vaccination'])
        virus_subtype = str(participant_row['Influenza_virus_subtype'])

        # Adjust symptom onset for non-index cases
        participant_symptom_onset = participant_row['DateofSymptomOnset']
        symptom_offset = participant_symptom_onset - index_symptom_onset

        measurements = []

        # Iterate through each day's value (adjusted day)
        for day in range(0, 13):
            adjusted_day = day - symptom_offset  # Adjust day based on onset difference
            day_column = f'Day{day}'

            try:
                value = 10 ** float(participant_row[day_column])
                if value == 0.1:
                    value = 'missing'
            except ValueError:
                value = 'missing'


            if participant_row['Type'] == 'NPS+OPS' and value != 'missing':
                measurements.append({
                    "analyte": "NPSOPS",
                    "time": adjusted_day,
                    "value": value
                })

    
        participant_dict_unindex = {
            "attributes": {
                "householdid": household_id,
                "age": age,
                "sex": sex,
                "vaccination": vaccination,
                "virus_subtype": virus_subtype
            },
            "measurements": measurements
        }

        # Append to the participant list
        participant_list_unindex.append(participant_dict_unindex)

        # for participant in participant_list_unindex[:10]:
        #     print(participant)


# Merge index and non-index participant lists
participant_list = []

# Group participants by householdID using a dictionary
household_participants = {}

# Add index participants to the dictionary
for participant in participant_list_index:
    household_id = participant["attributes"]["householdid"]
    if household_id not in household_participants:
        household_participants[household_id] = []
    household_participants[household_id].append(participant)

# Add unindex participants to the dictionary
for participant in participant_list_unindex:
    household_id = participant["attributes"]["householdid"]  # Ensure consistent key name
    if household_id not in household_participants:
        household_participants[household_id] = []
    household_participants[household_id].append(participant)

# Flatten the dictionary into a single list
for household_id, participants in household_participants.items():
    for participant in participants:
        participant_list.append(participant)

# Remove householdID from attributes in the final list
for participant in participant_list:
    participant["attributes"].pop("householdid", None)  
```

Finally, the data is formatted and output as a YAML file.

```python
tim2016 = dict(
    title= "Individual Correlates of Infectivity of Influenza A Virus Infections in Households",
    doi= "10.1371/journal.pone.0154418",
    description=folded_str('The community-based study reports the transmission of Influenza A virus within households during December 2008, focusing on viral loads in index cases and their household contacts to identify secondary infections. The study includes oropharynx and nasopharynx samples.\n'),   
    analytes=dict(
        NPSOPS=dict(description=folded_str("Influenza A virus RNA gene copy concentration in oropharynx and nasopharynx samples. Paired nasal and throat swabs were pooled in viral transport medium immediately after collection and subsequently processed for quantitative reverse transcription PCR to detect and quantify influenza A virus. The concentration was quantified in gene copies per milliliter. M gene was used as the gene target for the PCR assay, as described in a referenced study (Chan et al., 2008).\n"), 
            specimen=["nasopharyngeal_swab", "oropharyngeal_swab"],
            biomarker="influenza",
            gene_target="M gene", 
            limit_of_quantification="unknown",
            limit_of_detection=900, 
            unit="gc/mL",
            reference_event="symptom onset"
        )
    ),
    participants=participant_list
)

# Helper function to convert numpy types to native Python types
def clean_data(data):
    if isinstance(data, dict):
        return {key: clean_data(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [clean_data(item) for item in data]
    elif isinstance(data, np.bool_):  # Handle numpy booleans
        return bool(data)
    elif isinstance(data, np.generic):  # Check for numpy types
        return data.item()  # Convert numpy to Python native type
    return data

# Clean the `tim2016` data
tim2016_cleaned = clean_data(tim2016)


with open("tim2016individual.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(tim2016_cleaned, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close() 
```    