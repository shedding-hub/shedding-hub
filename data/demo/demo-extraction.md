```python
!pip install openpyxl
```

```python
import pandas as pd
import yaml
```

```python
raw = pd.read_excel("/Users/till/Downloads/CombinedDataset.xlsx")
raw = raw[raw.StudyNum == 17]
raw
```

```python
for patient, subset in raw.groupby("PatientID"):
    pass
```

```python
subset
```

```python
with open("demo.yaml", "w") as fp:
    yaml.dump({
        "doi": "askdsk",
        "title": "akaka",
        "participants": [],
    }, fp)
```

```python
yaml.dump({
        "doi": "askdsk",
        "title": "akaka",
        "participants": [],
    })
```

```python

```
