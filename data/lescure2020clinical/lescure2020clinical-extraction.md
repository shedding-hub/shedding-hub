# Extraction of the raw data for Lescure et al. (2020) using the data published by Goyal et al. (2020)

[Goyal et al. (2020)](https://www.science.org/doi/10.1126/sciadv.abc7112) developed mathematical models to project multiple therapeutic approaches. The authors used four datasets of SARS-CoV-2 shedding in the absence of effective treatment to develop and validate a mathematical model. These data included including a dataset of 4 patients from France ([Lescure et al. 2020](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30200-0/fulltext)).

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

Raw data ([Viral_Loads](https://github.com/shedding-hub/shedding-hub/blob/main/data/lescure2020clinical/Viral_Loads.csv)), which is stored on [Shedding Hub](https://github.com/shedding-hub), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
#load the data;
Lescure2020 = pd.read_csv("Viral_Loads.csv")
#Subset the data for Lescure et al. 2020;
Lescure2020 = Lescure2020[Lescure2020["cov_study"]==4]
#some data cleaning to match the schema;
Lescure2020.loc[Lescure2020["cens"]==1,"VL"]="negative"
```

Finally, the data is formatted and output as a YAML file.

```python
#Patient information from Table 1 and Figure 2 in Lescure et al. (2020)
patients_info = pd.DataFrame(dict(ID=['E1','E3', 'E4', 'E5'],
                                  age=[31, 80, 30, 46],
                                  sex=['male','male','female','female'],
                                  day0=[6,8,2,2]))

participant_list = [dict(attributes=dict(age=patients_info.loc[patients_info["ID"]==i,"age"].item(),
                                         sex=patients_info.loc[patients_info["ID"]==i,"sex"].item(),),
                         measurements=[dict(analyte="naso_swab_SARSCoV2",
                                             time=int(Lescure2020.loc[j,"dao"].item()+patients_info.loc[patients_info["ID"]==i,"day0"].item()),
                                             value=Lescure2020.loc[j,"VL"]) for j in Lescure2020.loc[Lescure2020["ID"]==i].index]) for i in pd.unique(Lescure2020["ID"])]

lescure2020clinical = dict(title="Clinical and virological data of the first cases of COVID-19 in Europe: a case series",
                              doi="10.1016/S1473-3099(20)30200-0",
                              description=folded_str('The authors followed five patients admitted to Bichat-Claude Bernard University Hospital (Paris, France) and Pellegrin University Hospital (Bordeaux, France) and diagnosed with COVID-19 by semi-quantitative RT-PCR on nasopharyngeal swabs. We assessed patterns of clinical disease and viral load from different samples (nasopharyngeal and blood, urine, and stool samples), which were obtained once daily for 3 days from hospital admission, and once every 2 or 3 days until patient discharge. Stool samples only have positive and negative results (currently not included in this data). The data was obtained from Goyal et al. 2020 for the nasopharyngeal swab results in 4 patients.\n'),
                              analytes=dict(naso_swab_SARSCoV2=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration in nasopharyngeal swab samples. The concentration were quantified in genome copies per swab.\n"),
                                                                    specimen="nasopharyngeal_swab",
                                                                    biomarker="SARS-CoV-2",
                                                                    gene_target="RdRp or E",
                                                                    limit_of_quantification=100,
                                                                    limit_of_detection=1,
                                                                    unit="gc/swab",
                                                                    reference_event="symptom onset")),
               participants=participant_list)

with open("lescure2020clinical.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(lescure2020clinical, outfile, default_style=None, default_flow_style=False, sort_keys=False)
outfile.close()
```

```python

```
