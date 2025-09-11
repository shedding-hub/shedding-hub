# Extraction of Jones et al. (2021)

[Jones et al. (2021)](https://www.science.org/doi/10.1126/science.abi5273) analyzed viral load data from 25,381 German cases, including 9519 hospitalized patients, 6110 PAMS cases from walk-in test centers, 1533 B.1.1.7 variant infections, and the viral load time series of 4434 (mainly hospitalized) patients. Viral load results were then combined with estimated cell culture isolation probabilities, producing a clinical proxy estimate of infectiousness. Data were obtained from the ['viral-load-with-negatives.tsv.zip'](https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/blob/20210614/data/viral-load-with-negatives.tsv.zip) dataset in the GitHub repository of [Jones et al. (2021)](https://www.science.org/doi/10.1126/science.abi5273). Samples lacking a known onset date or a detected positive test result were filtered out.

First, we `import` python modules needed:

```python
#import modules;
import yaml
import pandas as pd
from shedding_hub import folded_str, literal_str
```

Raw data ['viral-load-with-negatives.tsv'](https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/blob/20210614/data/viral-load-with-negatives.tsv.zip), which is stored on [Shedding Hub](https://github.com/shedding-hub), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

Each individual has only one positive test at record. 
```python
df = pd.read_csv("/Users/fionaxiao/Desktop/💩/SheddingHub/jones2021estimating/data/viral-load-with-negatives.tsv", sep="\t")
df_filtered = df.groupby('personHash').filter(lambda g: not (g['Ct'] == -1.0).all())
df_filtered['viral_load_extracted'] = 10 ** df_filtered['log10Load']

df_individual_positive = (
    df_filtered.assign(valid=(df_filtered["Ct"] != -1.0))
    .groupby("personHash", as_index=False)
    .agg(valid_Ct_count=("valid", "sum"))
)

df_individual_positive["valid_Ct_count"].unique()
```

We split into two datasets based on if people have onset information. 
df_with_onset → people with valid onset info + reference_date_onset
df_without_onset → people without onset info + reference_date_confirmation

```python
df_filtered['Onset'] = pd.to_datetime(df_filtered['Onset'], format='%Y-%m-%d', errors='coerce')
df_filtered['Date']  = pd.to_datetime(df_filtered['Date'],  format='%Y-%m-%d', errors='coerce')

g = df_filtered.groupby('personHash', group_keys=False)
has_onset = g['Onset'].transform(lambda s: s.notna().any())

# ----------------------
# Dataset 1: With onset
# ----------------------
df_with_onset = df_filtered[has_onset].copy()
df_with_onset['reference_date_onset'] = (df_with_onset['Date'] - df_with_onset['Onset']).dt.days

# -------------------------
# Dataset 2: Without onset
# -------------------------
df_without_onset = df_filtered[~has_onset].copy()

first_positive_date_by_person = (
    df_without_onset.groupby('personHash')
    .apply(lambda sub: sub.loc[sub['Ct'] != -1.0, 'Date'].iloc[0]
           if (sub['Ct'] != -1.0).any() else pd.NaT)
)

df_without_onset['confirmation_date'] = df_without_onset['personHash'].map(first_positive_date_by_person)
df_without_onset['reference_date_confirmation'] = (
    (df_without_onset['Date'] - df_without_onset['confirmation_date']).dt.days
)
```

reference_date is the number of days between the test date and the patient’s reported symptom-onset date. Large positive and negative values mean the swab was taken well before or after they were infected. Those samples are irrelevant, so we drop any record with reference_date < -7 and >30.

```python
example = df_with_onset.loc[df_with_onset['personHash'] == "8bb13ce1d8b40383949ff842074f4040"]
example = example.reset_index(drop=True)
print (example)
```

```python 
# Without onset
df_filtered_without_onset = df_without_onset.loc[
    (df_without_onset['reference_date_confirmation'] >= -7) &
    (df_without_onset['reference_date_confirmation'] <= 30)
].copy()

# With onset
df_filtered_with_onset = df_with_onset.loc[
    (df_with_onset['reference_date_onset'] >= -7) &
    (df_with_onset['reference_date_onset'] <= 30)
].copy()
```

```python

participants = []

for patient_id, group in df_filtered_with_onset.groupby('personHash'):

    sex = ('female' if group['Gender'].iloc[0] == 'F'
                else 'male' if group['Gender'].iloc[0] == 'M'
                else 'unknown')

    participant = {
        'attributes': {
            'age': int(group['Age'].iloc[0]),
            'sex': sex, 
            'lineage': 'B.1.1.7' if group['Gender'].iloc[0] == 1 else 'unknown' 
        },
        'measurements': []
    }

    for _, row in group.iterrows():
        if row['Ct'] == -1:
            value = "negative"
        else:
            value = row['viral_load_extracted']
        measurementN = {
            'analyte': 'swab_SARSCoV2_onset',
            'time': row['reference_date_onset'],
            'value': value
        }
        participant['measurements'].append(measurementN)
    
    participants.append(participant)
    
for patient_id, group in df_filtered_without_onset.groupby('personHash'):

    sex = ('female' if group['Gender'].iloc[0] == 'F'
                else 'male' if group['Gender'].iloc[0] == 'M'
                else 'unknown')

    participant = {
        'attributes': {
            'age': int(group['Age'].iloc[0]),
            'sex': sex, 
            'lineage': 'B.1.1.7' if group['Gender'].iloc[0] == 1 else 'unknown' 
        },
        'measurements': []
    }

    for _, row in group.iterrows():
        if row['Ct'] == -1:
            value = "negative"
        else:
            value = row['viral_load_extracted']
        measurementN = {
            'analyte': 'swab_SARSCoV2_confirmation',
            'time': row['reference_date_confirmation'],
            'value': value
        }
        participant['measurements'].append(measurementN)
    
    participants.append(participant)
```

Finally, the data is formatted and output as a YAML file.

```python
jones2021estimating = dict(title="Estimating infectiousness throughout SARS-CoV-2 infection course",
               doi="10.1126/science.abi5273",
               description=folded_str("This study analyzed viral load data from 25,381 German cases, including 9519 hospitalized patients, 6110 PAMS cases from walk-in test centers, 1533 B.1.1.7 variant infections, and the viral load time series of 4434 (mainly hospitalized) patients. Viral load results were then combined with estimated cell culture isolation probabilities, producing a clinical proxy estimate of infectiousness. Data were obtained from the 'viral-load-with-negatives.tsv' dataset in the GitHub repository of Jones et al. (2021). Samples lacking a known onset date or a detected positive test result were filtered out.\n"),
               analytes=dict(swab_SARSCoV2_onset=dict(description=folded_str("The reference event is symtom onset. Viral load is semiquantitative, estimating RNA copies per entire swab sample, whereas only a fraction of the volume can reach the test tube. The quantification is based on a standard curve and derive a formula in which RT-PCR cycle threshold values are converted to viral loads. Viral load is estimated from the cycle threshold (Ct) value using the empirical formulae 14.159 - (Ct x 0.297) for the Roche Light Cycler 480 system and 15.043 - (Ct x 0.296) for the Roche cobas 6800/8800 systems. An estimated 3 percent of our samples were from the lower respiratory tract. These were not removed from the dataset because of their low frequency and the fact that the first samples for patients are almost universally swab samples. Samples from the lower respiratory tract are generally taken from patients only after intubation, by which point viral loads have typically fallen.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection= "unknown",
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        unit="gc/swab",
                                        reference_event="symptom onset",),
                            swab_SARSCoV2_confirmation = dict(description=folded_str("The reference event is confirmation date. Viral load is semiquantitative, estimating RNA copies per entire swab sample, whereas only a fraction of the volume can reach the test tube. The quantification is based on a standard curve and derive a formula in which RT-PCR cycle threshold values are converted to viral loads. Viral load is estimated from the cycle threshold (Ct) value using the empirical formulae 14.159 - (Ct x 0.297) for the Roche Light Cycler 480 system and 15.043 - (Ct x 0.296) for the Roche cobas 6800/8800 systems. An estimated 3 percent of our samples were from the lower respiratory tract. These were not removed from the dataset because of their low frequency and the fact that the first samples for patients are almost universally swab samples. Samples from the lower respiratory tract are generally taken from patients only after intubation, by which point viral loads have typically fallen.\n"),
                                        limit_of_quantification="unknown",
                                        limit_of_detection= "unknown",
                                        specimen="nasopharyngeal_swab", 
                                        biomarker="SARS-CoV-2", 
                                        unit="gc/swab",
                                        reference_event="confirmation date",)
                                        ),
                participants = participants
                                        )

with open("jones2021estimating.yaml","w") as outfile:
    outfile.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(jones2021estimating, outfile, default_flow_style=False, sort_keys=False)

```
