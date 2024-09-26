# Extraction of the raw data for Young et al. (2020) using the data published by Goyal et al. (2020)

[Goyal et al. (2020)](https://www.science.org/doi/10.1126/sciadv.abc7112) developed mathematical models to project multiple therapeutic approaches. The authors used four datasets of SARS-CoV-2 shedding in the absence of effective treatment to develop and validate a mathematical model, which includes a dataset of 11 patients from Singapore ([Young et al. 2020](https://jamanetwork.com/journals/jama/fullarticle/2762688)).

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
import numpy as np

#functions to add folded blocks and literal blocks;
class folded_str(str): pass
class literal_str(str): pass

def folded_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='>')
def literal_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(folded_str, folded_str_representer)
yaml.add_representer(literal_str, literal_str_representer)
```

Raw data ([Viral_Loads](https://github.com/shedding-hub/shedding-hub/blob/main/data/young2020epidemiologic/Viral_Loads.csv)), which is stored on [Shedding Hub](https://github.com/shedding-hub), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
#load the data;
Young2020 = pd.read_csv("Viral_Loads.csv")
#Subset the data for Young et al. 2020;
Young2020 = Young2020[Young2020["cov_study"]==1]
#some data cleaning to match the schema;
Young2020.loc[Young2020["cens"]==1,"VL"]="negative"
```

Finally, the data is formatted and output as a YAML file.

```python
participant_list = [dict(measurements=[dict(analyte="naso_swab_SARSCoV2",
                                             time=int(Young2020.loc[j,"dao"].item()),
                                             value=Young2020.loc[j,"VL"]) for j in Young2020.loc[Young2020["ID"]==i].index]) for i in pd.unique(Young2020["ID"])]

young2020epidemiologic = dict(title="Epidemiologic Features and Clinical Course of Patients Infected With SARS-CoV-2 in Singapore",
                              doi="10.1001/jama.2020.3204",
                              description=folded_str('Clinical, laboratory, and radiologic data were collected, including PCR cycle threshold values from nasopharyngeal swabs and viral shedding in blood, urine, and stool. Clinical course was summarized, including requirement for supplemental oxygen and intensive care and use of empirical treatment with lopinavir-ritonavir. The numbers of positive stool, blood, and urine samples were small. The data was obtained from Goyal et al. 2020 for the nasopharyngeal swab results in 11 patients.\n'),
                              analytes=dict(naso_swab_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration in nasopharyngeal swab samples. The concentration were quantified in genome copies per swab.\n"),
                                                                    specimen="nasopharyngeal_swab",
                                                                    biomarker="SARS-CoV-2",
                                                                    gene_target="N, E, and ORF1lab",
                                                                    limit_of_quantification="unknown",
                                                                    limit_of_detection="unknown",
                                                                    unit="gc/swab",
                                                                    reference_event="confirmation date")),
               participants=participant_list)

with open("young2020epidemiologic.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(young2020epidemiologic, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close()
```
