# Extraction of Perera et al. (2020) using the data published by Challenger et al. (2022) and Goyal et al. (2020)

[Perera et al. (2020)](https://wwwnc.cdc.gov/eid/article/26/11/20-3219_article) investigated 68 respiratory specimens from 35 coronavirus disease patients in Hong Kong, 32 of whom had mild disease. The specimens comprised 46 combined nasopharyngeal aspirates and throat swabs, 2 nasopharyngeal aspirates alone, 4 combined nasopharyngeal and throat swabs, 3 nasopharyngeal swabs alone, 11 sputum samples, and 2 saliva samples. The patient origin of each sample is unknown. Data on RNA loads were obtained from the ['Culture_probability_data_wild_type.csv'](https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/blob/20210614/data/Culture_probability_data_wild_type.csv) dataset in the GitHub repository of [Jones et al. (2021)](https://www.science.org/doi/10.1126/science.abi5273).


First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
import numpy as np
from shedding_hub import folded_str, literal_str
```

Raw data ([Culture_probability_data_wild_type.csv](https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/blob/20210614/data/Culture_probability_data_wild_type.csv)), which is stored on [Shedding Hub](https://github.com/shedding-hub), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
df = pd.read_csv("Culture_probability_data_wild_type.csv")
ranawaka = df[df["original_study"] == "ranawaka"]
ranawaka['viral_load_extracted'] = 10 ** ranawaka['log10_viral_load_extracted']
```

```python
participants = []

for _, row in ranawaka.iterrows():
    if row['log10_viral_load_extracted'] <= 1:
        val = "negative"
    else:
        val = row['viral_load_extracted']

    measurement = {
        'analyte': 'SARSCoV2_N',
        'time': row['day_from_symptom_onset_rounded'],
        'value': val
    }
    participant = {'measurements': [measurement]}
    participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.

```python
perera2020sars = dict(title="SARS-CoV-2 Virus Culture and Subgenomic RNA for Respiratory Specimens from Patients with Mild Coronavirus Disease",
               doi="10.3201/eid2611.203219",
               description=folded_str("This study investigated 68 respiratory specimens from 35 coronavirus disease patients in Hong Kong, 32 of whom had mild disease. The specimens comprised 46 combined nasopharyngeal aspirates and throat swabs, 2 nasopharyngeal aspirates alone, 4 combined nasopharyngeal and throat swabs, 3 nasopharyngeal swabs alone, 11 sputum samples, and 2 saliva samples. The patient origin of each sample is unknown. Data on RNA loads were obtained from the 'Culture_probability_data_wild_type.csv' dataset in the GitHub repository of Jones et al. (2021).\n"),
               analytes=dict(SARSCoV2_N=dict(description=folded_str("The specimens comprised 46 combined nasopharyngeal aspirates and throat swabs, 2 nasopharyngeal aspirates alone, 4 combined nasopharyngeal and throat swabs, 3 nasopharyngeal swabs alone, 11 sputum samples, and 2 saliva samples. The limit of detection for viral N gene RNA was 10 copies/mL. RNA was extracted using the QIAamp Viral RNA Extraction Kit and then tested by RT quantitative PCR targeting the SARSCoV2 nucleoprotein (N) gene. Serial dilutions of a copy number control plasmid DNA were included in each RT qPCR run to construct a standard curve correlating cycle threshold values with gene copy numbers in the samples.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection= 10,
                                        specimen=["nasopharyngeal_swab", "nasopharyngeal_aspirates", "throat_swabs", "sputum", "saliva"], 
                                        biomarker="SARS-CoV-2", 
                                        gene_target="N", 
                                        unit="gc/mL",
                                        reference_event="symptom onset",)
                                        ),
                participants = participants
                                        )

with open("perera2020sars.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(perera2020sars, outfile, default_flow_style=False, sort_keys=False)

```
