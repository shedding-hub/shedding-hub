# Extraction of the raw data for Lescure et al. (2020) using the data published by Goyal et al. (2020)

[Goyal et al. (2020)](https://www.science.org/doi/10.1126/sciadv.abc7112) developed mathematical models to project multiple therapeutic approaches. The authors used four datasets of SARS-CoV-2 shedding in the absence of effective treatment to develop and validate a mathematical model. These data included including a dataset of 4 patients from France ([Lescure et al. 2020](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30200-0/fulltext)).

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd

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
#correct some doa errors; #the variable "dao" is "the day after first positive" used in Goyal et al. (2020).
#The last two records had incorrect days (13 and 15) based on Figure 2 in Lescure et al. (2020). We corrected them with days 14 and 16.
Lescure2020.loc[(Lescure2020["ID"]=="E3") & (Lescure2020["dao"]==13),"dao"]=14
Lescure2020.loc[(Lescure2020["ID"]=="E3") & (Lescure2020["dao"]==15),"dao"]=16
#some data cleaning to match the schema;
Lescure2020.loc[Lescure2020["cens"]==0,"VL"]=10**Lescure2020.loc[Lescure2020["cens"]==0,"VL"]
Lescure2020.loc[Lescure2020["cens"]==1,"VL"]="negative"
#The positive but not quantifiable data are from Figure 2 in the paper;
Lescure2020.loc[(Lescure2020["ID"]=="E3") & (Lescure2020["dao"].isin([14,16])),"VL"]="positive"
Lescure2020.loc[(Lescure2020["ID"]=="E5") & (Lescure2020["dao"]==8),"VL"]="positive"

#The `Viral_Loads.csv` is the data from Goyal et al. (2020), which only included 4 patients from Lescure et al. (2020). 
#Data from Patient 2 were added based on Figure 2 in Lescure et al. (2020).
P2 = {
    "ID": ['E2','E2','E2','E2','E2'],
    "dao": [0,2,5,11,13],
    "VL": ['positive','negative','negative','negative','negative']
}
dat_P2 = pd.DataFrame(P2)

Lescure2020 = pd.concat([Lescure2020,dat_P2])
```

Finally, the data is formatted and output as a YAML file.

```python
#Patient information from Table 1 and Figure 2 in Lescure et al. (2020)
patients_info = pd.DataFrame(dict(ID=['E1','E2','E3', 'E4', 'E5'],
                                  age=[31, 48, 80, 30, 46],
                                  sex=['male','male','male','female','female'],
                                  day0=[6,9,8,2,2])) #"day0" is "the sampling day after onset".

participant_list = [dict(attributes=dict(age=patients_info.loc[patients_info["ID"]==i,"age"].item(),
                                         sex=patients_info.loc[patients_info["ID"]==i,"sex"].item(),),
                         measurements=[dict(analyte="naso_swab_SARSCoV2",
                                            #the variable "dao" is "the day after first positive" used in Goyal et al. (2020) and "day0" is "the sampling day after onset". "dao" + "day0" will be "day after symptom onset".
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
