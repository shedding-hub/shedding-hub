# Extraction for Xing et al. (2020)

[Xing et al. (2020)](https://doi.org/10.1016/j.jmii.2020.03.021) compared the dynamic changes of severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) RNA in respiratory and fecal specimens in children with coronavirus disease 2019 (COVID-19). Only one case with complete demographic information (e.g., age, sex, etc.) and specimen test results was included in the paper. Throat samples were recorded as oropharyngeal samples in the standardized dataset.

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
from shedding_hub import folded_str

```

The raw data was extracted from Figure 4 (see below) in Xing et al. (2020) and entered in the markdown (.md) file.

![image](patient_data.png)

```python
# Define a dictionary containing patient demographic information (Case Number, Sex, Age) from [Table 1] in Xing et al. (2020).
patient_demo_info = {
    'Case1': {'CaseNumber': 1, 'Sex': 'male', 'Age': 1.5}
}

# List case1's Ct_Values for fecal specimen and throat swab along with time (days from admission)
time = [1, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36]
ct_fecal = [12.65, 12.31, 15.81, 13.97, 18, 24.37, 28.29, 34.67, 37.05, 38.76, 40, 40, 40, 40, 40, 40, 40]
ct_throat = [18.17, 23.97, 20.53, 25.91, 28.36, 34.83, 32.17, 38.19, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]

# Initialize an empty list to store participant information
participant_list = []

# Loop through each unique patient in the patient_demo_info dictionary
for case_key, case_info in patient_demo_info.items():
    age = float(case_info['Age'])  # Convert to Python float (since age is 1.5)
    sex = str(case_info['Sex'])  # Convert to Python str

    measurements = []
    # Iterate over each time point
    for index, t in enumerate(time):  # Use enumerate to get both index and value
        # Append oropharyngeal swab data for all time points
            measurements.append({
                "analyte": "oropharyngeal_swab_SARSCoV2",
                "time": t,
                "value": ct_throat[index] if ct_throat[index] < 40 else 'negative'
            })
    for index, t in enumerate(time):
        # Append fecal data starting from the third time point
        if index >= 2 and (index - 2) < len(ct_fecal):  
            measurements.append({
                "analyte": "stool_SARSCoV2",
                "time": t,
                "value": ct_fecal[index - 2] if ct_fecal[index - 2] < 40 else 'negative'
            })

    participant_dict = {
        "attributes": {
            "age": age,
            "sex": sex
        },
        "measurements": measurements
    }
    participant_list.append(participant_dict)  # Append the participant's data to the list

```

Finally, the data is formatted and output as a YAML file.

```python
Xing2020 = dict(title="Prolonged viral shedding in feces of pediatric patients with coronavirus disease 2019",
               doi="10.1016/j.jmii.2020.03.021",
               description=folded_str('This study characterized the dynamic profiles of SARS-CoV-2 shedding in respiratory and fecal specimens in three children with COVID-19. However, complete specimen test value data were available for only one child. A total of 19 oropharyngeal swab specimens and 17 fecal specimens were collected.\n'),
               analytes=dict(stool_SARSCoV2=dict(description=folded_str("Presence of SARS-CoV-2 RNA was detected by RT-PCR in stool samples. The PCR assay simultaneously amplified two target genes of SARS-CoV-2 included open reading frame 1ab (ORF1ab) and nucleocapsid protein (N). A cycle threshold (Ct) value no more than 40 with evident amplification curve was considered as a positive test, and a value of 40 indicated the virus was molecularly undetectable.\n"),
                                                    specimen="stool",
                                                    biomarker="SARS-CoV-2",
                                                    gene_target="ORF1ab, N", 
                                                    limit_of_quantification='unknown', 
                                                    limit_of_detection=40,
                                                    unit="cycle threshold",
                                                    reference_event="hospital admission"), 
                             oropharyngeal_swab_SARSCoV2=dict(description=folded_str("Presence of SARS-CoV-2 RNA was detected by RT-PCR in oropharyngeal_swab samples. The PCR assay simultaneously amplified two target genes of SARS-CoV-2 included open reading frame 1ab (ORF1ab) and nucleocapsid protein (N). A cycle threshold (Ct) value no more than 40 with evident amplification curve was considered as a positive test, and a value of 40 indicated the virus was molecularly undetectable.\n"),
                                              specimen="oropharyngeal_swab",
                                              biomarker="SARS-CoV-2",
                                              gene_target="ORF1ab, N",
                                              limit_of_quantification="unknown",
                                              limit_of_detection=40,
                                              unit="cycle threshold",
                                             reference_event="hospital admission")), 

               participants=participant_list)

with open("xing2020prolonged.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(Xing2020, outfile, default_style=None, default_flow_style=False, sort_keys=False)
```

