---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.4
  kernelspec:
    display_name: python_clean
    language: python
    name: python3
---

Article title: Suppression of a SARS-CoV-2 outbreak in the Italian municipality of Vo
This study was conducted in the Italian municipality of Vo. Lockdown was implemented
  after first death of pneumonia was reported. Two surveys and virus tests were conducted
  with the first survey near the start of lockdown and the second one at the end of
  lockdown.

```python
import yaml
import pandas as pd
import numpy as np
import textwrap
import datetime


#functions to add folded blocks and literal blocks;
from shedding_hub import folded_str
```

Data are stored in 2 worksheets in excel file. 'df_demographics' extracts the demographic info of patients; and 'df_pcr' extracts the pcr results of patients.

```python
#Read Excel file into the environment
df_demographics = pd.read_excel('anonymised_data_public_final.xlsx', sheet_name='anonymised_dataset')
df_pcr = pd.read_excel('anonymised_data_public_final.xlsx', sheet_name='RT_PCR_DATA')

```

```python
# Filter the useful columns in df_demographics
df_demographics2 = df_demographics.loc[:,'id':'gender'].join(df_demographics['first_symptoms_date'])

```

```python
#Filter and rename the columns in df_pcr
df_pcr2 = df_pcr.loc[:,[ 'id','RT_PCR_Genome_Equivalents_First_survey_RdRp','RT_PCR_Genome_Equivalents_First_survey_E',
                        'RT_PCR_Genome_Equivalents_Second_survey_RdRp','RT_PCR_Genome_Equivalents_Second_survey_E',
                        'symptomatic_at_first_sampling','symptomatic_at_follow_up',]]
df_pcr2.columns = ['id','RdRp_1','E_1','RdRp_2','E_2','symptom_1','symptom_2']
```

```python
#Join pcr data with the demographic data
df = df_pcr2.merge(df_demographics2,how = 'left',on = 'id')

```

The first survey occurred between 21 and 29 February 2020 and the
second survey occurred on 7 March 2020.

```python
#Recode the gender.
df['gender'] = df['gender'].replace({'M':'male','F':'female'})

```

```python
print(df['id'].duplicated().value_counts()) #No Duplicate entries

#Reshape the dataframe.
df = pd.wide_to_long(df, 
                          stubnames=['RdRp', 'E', 'symptom'], 
                          i=['id','age_group','gender','first_symptoms_date'], 
                          j='Round', 
                          sep='_', 
                          suffix='(1|2)')
df = df.reset_index()
```

Test date were stored as many columns with column names as test date, and the values are test results (Pos, Neg, or NA). We need to extract and reshape them into single field.

```python
#Retrieve the test date and result information.
df_test = df_demographics.loc[:,datetime.datetime(2020, 2, 21, 0, 0):datetime.datetime(2020, 3, 10, 0, 0)].join(df_demographics['id'])
df_test = pd.melt(df_test,id_vars=['id'],var_name = 'test_date',value_name = 'results')
df_test = df_test[~df_test['results'].isnull()]

```

```python
#Filter to id with available PCR results.
df_test = df_test[df_test['id'].isin(df['id'])]
#Convert to date
df_test['test_date'] = pd.to_datetime(df_test['test_date'])
df_test['test_date'].value_counts()
```

Some patients had more than 2 test dates, and we apply some algorithm to assign test dates to pcr values. 
1. The test dates were first grouped into round 1 and 2 by whether it is before 2020/03/06.
2. In each round, the first positive date and first negative date were calculated.
3. For PCR side, if any of the E and RdRp gene shows numbers, then the PCR result is regarded as 'Positive'. In each round, the date will be matched by the PCR result. (For example, if in round 1, the PCR result is 'Positive', the test date will be first positive date in round 1).

```python
#Separate round 1 and round 2 by date before 2020/03/06 or after.
df_test_1= df_test[df_test['test_date'] < datetime.datetime(2020, 3, 6)]
df_test_2 = df_test[df_test['test_date'] >= datetime.datetime(2020, 3, 6)]
```

```python
#Generate the first positive /neg day for both round 1 and round 2   
for type in ['Pos','Neg']:
    temp_m = pd.DataFrame()
    for round in range(1,3):
        if round == 1:
            temp = df_test_1
        if round == 2:
            temp = df_test_2
        temp = temp[temp['results'] == type]
        temp = temp.groupby('id')['test_date'].agg(min).reset_index() #Retrieve the minimum positive day
        temp = temp.rename(columns = {'test_date':'first_'+ type})
        temp['Round'] = round
        temp_m = pd.concat([temp_m,temp],axis = 0,ignore_index=True)
    df = df.merge(temp_m,how = 'left',on=['id','Round'])

```

```python
#Create the PCR results
df['pcr_results'] = 'Neg'
df.loc[~df['RdRp'].isnull() | ~df['E'].isnull(),'pcr_results'] = 'Pos' 
```

```python
#Assign the Test date
df.loc[df['pcr_results'] == 'Pos','test_date'] = df['first_Pos']
df.loc[df['pcr_results'] == 'Neg','test_date'] = df['first_Neg']
 
```

The reference date is determined with both first positive date and first symptom date.
If first symptom date is available, then the reference date is first symptom date.
If not available, we use the first positive date as reference date.
This reference date serve as the start point of the 'time' field of each test.

```python
#Retrieve the first positive date
df_test_pos = df_test[df_test['results'] == 'Pos']
df_test_pos = df_test_pos.groupby('id')['test_date'].agg(min).reset_index()
df_test_pos = df_test_pos.rename(columns = {'test_date' : 'first_pos_date_overall'})
df = df.merge(df_test_pos,on = 'id',how = 'left')
```

```python
#Try to create the Reference date.
df['reference_type'] = 'first_pos'
df.loc[~df['first_symptoms_date'].isnull(),'reference_type'] = 'symptom'

```

```python
#Create the numeric time
df.loc[df['reference_type'] == 'first_pos','time'] = (df['test_date'] - df['first_pos_date_overall']).dt.days
df.loc[df['reference_type'] == 'symptom','time'] = (df['test_date'] - df['first_symptoms_date']).dt.days

```

```python
#Exclude missing time (Which means no test was conducted.).
df = df.loc[~df['time'].isnull(),]
#Change the coding of negative pcr values to meet schema.
for i in ['RdRp','E']:
    df.loc[df[i].isnull(),i] = 'negative'
```

```python
#Enter the yaml writing stage
participant_list = []
for type in ['first_pos','symptom']:
      df_output = df[df['reference_type'] == type]    
      temp = [dict(attributes=dict(age_group=df_output.loc[df_output.loc[df_output["id"]==i].index[0],'age_group'],
                                   sex=df_output.loc[df_output.loc[df_output["id"]==i].index[0],'gender']),
                                   measurements=[dict(analyte='RdRp_'+type,
                                                time=df_output.loc[(df_output['Round'] == j) & (df_output['id'] == i),"time"].item(),
                                                value=df_output.loc[(df_output['Round'] == j) & (df_output['id'] == i),"RdRp"].item()) for j in np.unique(df_output.loc[df_output['id'] == i,'Round'])] +
                                                [dict(analyte="E_"+type,
                                                time=df_output.loc[(df_output['Round'] == j) & (df_output['id'] == i),"time"].item(),
                                                value=df_output.loc[(df_output['Round'] == j) & (df_output['id'] == i),"E"].item()) for j in np.unique(df_output.loc[df_output['id'] == i,'Round'])]) for i in df_output['id']]
      participant_list.extend(temp)


```

```python
lavezzo2020 = dict(title="Suppression of a SARS-CoV-2 outbreak in the Italian municipality of Vo",
               doi="10.1038/s41586-020-2488-1",
               description=folded_str('This study was conducted in the Italian municipality of Vo. Lockdown was implemented after first death of pneumonia was reported. Two surveys and virus tests were conducted with the first survey near the start of lockdown and the second one at the end of lockdown\n'),
               analytes=dict(RdRp_first_pos=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration calculated from evaluation of RdRp gene, and the reference event is first positive day.\n"),
                                                    specimen="oropharyngeal_swab",
                                                    biomarker="SARS-CoV-2",
                                                    gene_target="RdRp",
                                                    limit_of_quantification="unknown",
                                                    limit_of_detection="unknown",
                                                    unit="gc/mL",
                                                    reference_event="confirmation date"),
                             E_first_pos=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration calculated from evaluation of E gene, and the reference event is first positive day.\n"),
                                              specimen="oropharyngeal_swab",
                                              biomarker="SARS-CoV-2",
                                              gene_target="E",
                                              limit_of_quantification="unknown",
                                              limit_of_detection="unknown",
                                              unit="gc/mL",
                                              reference_event="confirmation date"),
                            RdRp_symptom=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration calculated from evaluation of RdRp gene, and the reference event is symptom onset day.\n"),
                                                specimen="oropharyngeal_swab",
                                                biomarker="SARS-CoV-2",
                                                gene_target="RdRp",
                                                limit_of_quantification="unknown",
                                                limit_of_detection="unknown",
                                                unit="gc/mL",
                                                reference_event="symptom onset"),
                             E_symptom=dict(description=folded_str("SARS-CoV-2 RNA genome copy concentration calculated from evaluation of E gene, and the reference event is symptom onset day.\n"),
                                              specimen="oropharyngeal_swab",
                                              biomarker="SARS-CoV-2",
                                              gene_target="E",
                                              limit_of_quantification="unknown",
                                              limit_of_detection="unknown",
                                              unit="gc/mL",
                                              reference_event="symptom onset")),
               participants=participant_list)

```

```python
with open("lavezzo2020suppression.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(lavezzo2020, outfile, default_style=None, default_flow_style=False, sort_keys=False)

```
