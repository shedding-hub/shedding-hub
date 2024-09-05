# Extraction for Natarajan et al. (2022)

[Natarajan et al. (2022)](https://doi.org/10.1016/j.medj.2022.04.001) quantified fecal shedding of different SARS-CoV-2 gene targets in feces for 113 individuals, and their data are published across several files in a [GitHub repository](https://github.com/alex-dahlen/lambda_fecal_shedding/). For reproducibility, we consider commit `7affa71`. They obtained diverse data using different sample collection kits, different gene targets, and different quantification methods.

```python
import pandas as pd
import yaml
import shedding_hub as sh
```

```python
# We load baseline data including demographics, the sample index mapping samples to
# participants.
base_url = "https://github.com/alex-dahlen/lambda_fecal_shedding/raw/7affa71/"
baseline = pd.read_csv(
    base_url + "Source%20Data/Source%20data_modified.xlsx%20-%20Baseline%20data.csv",
    index_col="Participant ID",
)
index = pd.read_csv(
    base_url + "Source%20Data/Source%20data_modified.xlsx%20-%20Index%20of%20clinical"
    "%20stool%20samples.csv",
    skiprows=3,
    parse_dates=["Date of sample collection", "Date of subject enrollment"],
    dayfirst=False,
)

# The raw data is grouped by assay type: sgRNA (subgenomic RNA by RT-qPCR), gRNA
# (genomic RNA by RT-qPCR), ddPCR (genomic RNA).
suffix_by_assay = {
    "gRNA-RT-qPCR": "Raw%20RT_qPCR%20gRNA%20data",
    "gRNA-ddPCR": "Raw%20ddPCR%20gRNA%20data",
    "sgRNA-RT-qPCR": "Raw%20RT-qPCR%20sgRNA%20data",
}

# Combine all the different datasets into one large raw dataset.
parts = []
for key, value in suffix_by_assay.items():
    raw = pd.read_csv(
        f"{base_url}Source%20Data/Source%20data_modified.xlsx%20-%20{value}.csv",
        skiprows=3,
    )
    # Subgenomic RNA was only targeting the N1 gene as described in the section "RT-qPCR
    # quantification of RNA".
    raw["analyte"] = (raw["Target gene"] if "Target gene" in raw else "N1") + "-" + key
    parts.append(raw)
raw = pd.concat(parts)
```

```python
# Merge the data and ensure we didn't gain or lose any samples. Then update the analyte
# definition.
merged = pd.merge(index, raw, on="RNA sample ID", how="outer")
assert merged.shape[0] == raw.shape[0]
merged["analyte"] = merged.analyte + "-" + merged["Stool preservative"]

participants = []
for subject_id, subset in merged.groupby("Subject study ID"):
    subject = baseline.loc[subject_id]

    participants.append(
        {
            "attributes": {
                "age": int(subject.Age),
                # This coding is based on manually comparing the {0, 1} coding of `Sex`
                # here with the `male` column in the `demographics.xlsx` spreadsheet.
                "sex": "female" if subject.Sex else "male",
                "study_arm": (
                    "test" if subject["Arm of study"] == "Lambda" else "control"
                ),
            },
            "measurements": [
                {
                    "analyte": row.analyte,
                    "value": (
                        "negative"
                        if (value := row["Viral RNA concentration (copies/μL)"]) == 0
                        else 1e3 * value
                    ),
                    "time": (
                        "unknown"
                        if pd.isnull(
                            days := (
                                row["Date of sample collection"]
                                - row["Date of subject enrollment"]
                            ).days
                        )
                        else days
                    ),
                    "sample_id": row["RNA sample ID"],
                }
                for _, row in subset.iterrows()
            ],
        }
    )

```

```python
description = sh.normalize_str(
    """
    Fecal samples were collected from 113 participants "in a randomized controlled study
    of Peg-interferon lambda-1a versus a placebo control for the treatment of mild to
    moderate COVID-19." A total of 7,563 measurements were taken for 679 samples using
    seven distinct assays; assays were typically run in duplicates.

    Samples were
    initially collected using OMNIGene GUT collection tubes (OG) but subsequently
    replaced with Zymo DNA/RNA shield fecal collection tubes (ZY) because of better
    performance. There are consequently 14 analytes following the naming convention
    `{target gene E, RdRP, N1, or N2}-{genomic RNA (gRNA) or subgenomic RNA
    (sgRNA)}-{quantification method RT-qPCR or ddPCR}-{collection tube OG or ZY}`.

    The reference event for temporal offsets in days is the day of enrollment.
    """
)
data = {
    "title": "Gastrointestinal symptoms and fecal shedding of SARS-CoV-2 RNA suggest "
    "prolonged gastrointestinal infection",
    "doi": "10.1016/j.medj.2022.04.001",
    "description": description,
    "analytes": {
        "E-gRNA-RT-qPCR-OG": {
            "description": "Concentration of genomic RNA of the E gene quantified "
            "using RT-qPCR for a sample collected using the OMNIGene GUT collection "
            "tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "E-gRNA-RT-qPCR-ZY": {
            "description": "Concentration of genomic RNA of the E gene quantified "
            "using RT-qPCR for a sample collected using the Zymo DNA/RNA shield fecal "
            "collection tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "RdRP-gRNA-RT-qPCR-OG": {
            "description": "Concentration of genomic RNA of the RdRp gene quantified "
            "using RT-qPCR for a sample collected using the OMNIGene GUT collection "
            "tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "RdRP-gRNA-RT-qPCR-ZY": {
            "description": "Concentration of genomic RNA of the RdRp gene quantified "
            "using RT-qPCR for a sample collected using the Zymo DNA/RNA shield fecal "
            "collection tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "N1-gRNA-RT-qPCR-OG": {
            "description": "Concentration of genomic RNA of the N1 gene quantified "
            "using RT-qPCR for a sample collected using the OMNIGene GUT collection "
            "tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "N1-gRNA-RT-qPCR-ZY": {
            "description": "Concentration of genomic RNA of the N1 gene quantified "
            "using RT-qPCR for a sample collected using the Zymo DNA/RNA shield fecal "
            "collection tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "N2-gRNA-RT-qPCR-OG": {
            "description": "Concentration of genomic RNA of the N2 gene quantified "
            "using RT-qPCR for a sample collected using the OMNIGene GUT collection "
            "tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "N2-gRNA-RT-qPCR-ZY": {
            "description": "Concentration of genomic RNA of the N2 gene quantified "
            "using RT-qPCR for a sample collected using the Zymo DNA/RNA shield fecal "
            "collection tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "N1-gRNA-ddPCR-OG": {
            "description": "Concentration of genomic RNA of the N1 gene quantified "
            "using ddPCR for a sample collected using the OMNIGene GUT collection "
            "tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "N1-gRNA-ddPCR-ZY": {
            "description": "Concentration of genomic RNA of the N1 gene quantified "
            "using ddPCR for a sample collected using the Zymo DNA/RNA shield fecal "
            "collection tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "E-gRNA-ddPCR-OG": {
            "description": "Concentration of genomic RNA of the E gene quantified "
            "using ddPCR for a sample collected using the OMNIGene GUT collection "
            "tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "E-gRNA-ddPCR-ZY": {
            "description": "Concentration of genomic RNA of the E gene quantified "
            "using ddPCR for a sample collected using the Zymo DNA/RNA shield fecal "
            "collection tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "N1-sgRNA-RT-qPCR-OG": {
            "description": "Concentration of sub-genomic RNA of the N1 gene quantified "
            "using RT-qPCR for a sample collected using the OMNIGene GUT collection "
            "tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
        "N1-sgRNA-RT-qPCR-ZY": {
            "description": "Concentration of sub-genomic RNA of the N1 gene quantified "
            "using RT-qPCR for a sample collected using the Zymo DNA/RNA shield fecal "
            "collection tube.",
            "specimen": "stool",
            "limit_of_detection": 1000,
            "limit_of_quantification": "unknown",
        },
    },
    "participants": participants,
}

with open("../data/natarajan2022.yaml", "w") as fp:
    fp.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(data, fp, sort_keys=False)
```
