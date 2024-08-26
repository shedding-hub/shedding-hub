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

#load the data;
Liu2024 = pd.read_csv("Liu2024.csv")

#sort by ID and date;
Liu2024 = Liu2024.sort_values(by=['subject','day_actual'])

#some data cleaning to match the scheme;
Liu2024.loc[Liu2024["Gender"]=="M","Gender"]="male"
Liu2024.loc[Liu2024["Gender"]=="F","Gender"]="female"
Liu2024.loc[Liu2024["Cohort"]=="Breakthrough","Cohort"]=True
Liu2024.loc[Liu2024["Cohort"]=="Unvaccinated","Cohort"]=False
Liu2024.loc[Liu2024["dpcr_result_class_N1"]=="Negative","gc_dryg_N1"]="negative"
Liu2024.loc[Liu2024["dpcr_result_class_PMMoV"]=="Negative","gc_dryg_PMMoV"]="negative"
Liu2024.loc[Liu2024["dpcr_result_class_mtDNA"]=="Negative","gc_dryg_mtDNA"]="negative"

participant_list = [dict(attributes=dict(age=float(Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"Age"]),
                                         sex=Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"Gender"],
                                         vaccinated=Liu2024.loc[Liu2024.loc[Liu2024["subject"]==i].index[0],"Cohort"]),
                         measurements=[dict(analyte="stool_SARSCoV2_N1",
                                             time=int(Liu2024.loc[j,"day_actual"].item()),
                                             value=Liu2024.loc[j,"gc_dryg_N1"]) for j in Liu2024.loc[Liu2024["subject"]==i].index] +
                                       [dict(analyte="stool_PMMoV",
                                             time=int(Liu2024.loc[j,"day_actual"].item()),
                                             value=Liu2024.loc[j,"gc_dryg_PMMoV"]) for j in Liu2024.loc[Liu2024["subject"]==i].index] +
                                       [dict(analyte="stool_mtDNA",
                                             time=int(Liu2024.loc[j,"day_actual"].item()),
                                             value=Liu2024.loc[j,"gc_dryg_mtDNA"]) for j in Liu2024.loc[Liu2024["subject"]==i].index]) for i in pd.unique(Liu2024["subject"])]

liu2024 = dict(title="Longitudinal Fecal Shedding of SARS-CoV-2, Pepper Mild Mottle Virus, and Human Mitochondrial DNA in COVID-19 Patients",
               doi="10.1101/2024.04.22.24305845",
               description=folded_str('The authors measured SARS-CoV-2, pepper mild mottle virus (PMMoV), and human mitochondrial DNA (mtDNA) in longitudinal stool samples collected from 42 COVID-19 patients for up to 42 days after the first sample collection date. Abundances were quantified using Digital PCR assays targeting the N1 gences. The data was contributed by the corresponding author.\n'),
               analytes=dict(stool_SARSCoV2_N1=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration in stool samples. The concentration were quantified in genome copies per dry weight of stool.\n"),
                                                    specimen="stool",
                                                    biomarker="SARS-CoV-2",
                                                    gene_target="N1",
                                                    limit_of_quantification=1000,
                                                    limit_of_detection="unknown",
                                                    unit="gc/dry gram",
                                                    reference_event="confirmation date"),
                             stool_PMMoV=dict(description=folded_str("PMMoV genome copy concentration in stool samples. The concentration were quantified in genome copies per dry weight of stool.\n"),
                                              specimen="stool",
                                              biomarker="PMMoV",
                                              limit_of_quantification="unknown",
                                              limit_of_detection="unknown",
                                              unit="gc/dry gram",
                                              reference_event="confirmation date"),
                             stool_mtDNA=dict(description=folded_str("mtDNA genome copy concentration in stool samples. The concentration were quantified in genome copies per dry weight of stool.\n"),
                                              specimen="stool",
                                              biomarker="mtDNA",
                                              limit_of_quantification="unknown",
                                              limit_of_detection="unknown",
                                              unit="gc/dry gram",
                                              reference_event="confirmation date")),
               participants=participant_list)

with open("Liu2024Longitudinal.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=.schema.yaml\n")
    yaml.dump(liu2024, outfile, default_flow_style=False, sort_keys=False)
outfile.close() 