# Extraction for Han et al. (2020)

[Han et al. (2020)](https://doi.org/10.1093/cid/ciaa447) reports SARS-CoV-2 viral loads in different specimen types for a neonate and her mother diagnosed on 2020-03-20. The study includes nasopharyngeal, oropharyngeal, stool, plasma, saliva, and urine samples. The viral load in the respiratory specimens gradually decreased with time and was undetectable after 17 days from the onset of symptoms.

First, we `import` python modules needed:

```python
import yaml
import pandas as pd
import numpy as np
```
Raw data, which is stored on [Shedding Hub](https://github.com/shedding-hub/shedding-hub/tree/main/data/han2020sequential), will be loaded and cleaned to match the most updated [schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml).

```python
limit_of_detection = 5.7e3
configs = [
    {
        "filename": "./mother_project.json",
        "attributes": {
            "age": "unknown",
            "sex": "female",
        },
    },
    {
        "filename": "./neonate_project.json",
        "attributes": {
            "age": 0.074, #0.074 represent 27 days as participant is a 27-day-old neonate
            "sex": "female",
        },
    },
]
analytes = {
    "nasopharynx_E": {
        "description": "SARS-CoV-2 RNA gene copy concentration in nasopharynx samples.",
        "specimen": "nasopharyngeal_swab",
        "biomarker": "SARS-CoV-2",
        "gene_target": "E",
        "limit_of_quantification": "unknown",
        "limit_of_detection": limit_of_detection,
        "unit": "gc/mL",
        "reference_event": "symptom onset"
    },
    "oropharynx_E": {
        "description": "SARS-CoV-2 RNA gene copy concentration in oropharynx samples.",
        "specimen": "oropharyngeal_swab",
        "biomarker": "SARS-CoV-2",
        "gene_target": "E",
        "limit_of_quantification": "unknown",
        "limit_of_detection": limit_of_detection,
        "unit": "gc/mL",
        "reference_event": "symptom onset"
    },
    "NPSOPS_E": {
        "description": "SARS-CoV-2 RNA gene copy concentration in oropharynx and "
        "nasopharynx samples.",
        "specimen": ["nasopharyngeal_swab", "oropharyngeal_swab"],
        "biomarker": "SARS-CoV-2",
        "gene_target": "E",
        "limit_of_quantification": "unknown",
        "limit_of_detection": limit_of_detection,
        "unit": "gc/mL",
        "reference_event": "symptom onset"
    },
    "plasma_E": {
        "description": "SARS-CoV-2 RNA gene copy concentration in plasma samples.",
        "specimen": "plasma",
        "biomarker": "SARS-CoV-2",
        "gene_target": "E",
        "limit_of_quantification": "unknown",
        "limit_of_detection": limit_of_detection,
        "unit": "gc/mL",
        "reference_event": "symptom onset"
    },
    "saliva_E": {
        "description": "SARS-CoV-2 RNA gene copy concentration in saliva samples.",
        "specimen": "saliva",
        "biomarker": "SARS-CoV-2",
        "gene_target": "E",
        "limit_of_quantification": "unknown",
        "limit_of_detection": limit_of_detection,
        "unit": "gc/mL",
        "reference_event": "symptom onset"
    },
    "urine_E": {
        "description": "SARS-CoV-2 RNA gene copy concentration in urine samples.",
        "specimen": "urine",
        "biomarker": "SARS-CoV-2",
        "gene_target": "E",
        "limit_of_quantification": "unknown",
        "limit_of_detection": limit_of_detection,
        "unit": "gc/mL",
        "reference_event": "symptom onset"
    },
    "sputum_E": {
        "description": "SARS-CoV-2 RNA gene copy concentration in sputum samples.",
        "specimen": "sputum",
        "biomarker": "SARS-CoV-2",
        "gene_target": "E",
        "limit_of_quantification": "unknown",
        "limit_of_detection": limit_of_detection,
        "unit": "gc/mL",
        "reference_event": "symptom onset"
    },
    "stool_E": {
        "description": "From manuscript: \"Viral RNA was detected using the PowerChek "
        "2019-nCoV real-time polymerase chain reaction kit (Kogene Biotech, Seoul, "
        "Korea) for amplification of the E gene and the RNA-dependent RNA polymerase "
        "region of the ORF1b gene, and quantified with a standard curve that was "
        "constructed using in vitro transcribed RNA provided from the European Virus "
        "Archive.\" However, data reported in Figure 1 explicitly refer to the E gene "
        "target.",
        "specimen": "stool",
        "biomarker": "SARS-CoV-2",
        "gene_target": "E",
        "limit_of_quantification": "unknown",
        "limit_of_detection": limit_of_detection,
        "unit": "gc/mL",
        "reference_event": "symptom onset"
    },
}
```

```python
data = {
    "title": "Sequential Analysis of Viral Load in a Neonate and Her Mother Infected "
    "With Severe Acute Respiratory Syndrome Coronavirus 2",
    "doi": "10.1093/cid/ciaa447",
    "description": "The study reports SARS-CoV-2 viral loads in different specimen "
    "types for a neonate and her mother diagnosed on 2020-03-20. The study includes "
    "nasopharyngeal, oropharyngeal, stool, plasma, saliva, and urine samples. Viral "
    "loads were extracted manually from Figure 1 using "
    "[WebPlotDigitizer](https://automeris.io).",
    "analytes": analytes,
    
}
for config in configs:
    with open(config["filename"]) as fp:
        part = yaml.safe_load(fp)

    measurements = []
    for subset in part["datasetColl"]:
        analyte = subset["name"]
        assert analyte in analytes, (analyte, config["filename"])
        for measurement in subset["data"]:
            time, log10 = measurement["value"]
            value = 10 ** log10
            measurements.append({
                "analyte": analyte,
                "value": value if value >= limit_of_detection else "negative",
                "time": round(time)
            })
    
    data.setdefault("participants", []).append({
        "attributes": config["attributes"],
        "measurements": measurements,
    })

with open("han2020sequential.yaml", "w") as fp:
    fp.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump(data, fp, default_style=None, default_flow_style=False, sort_keys=False)
```
