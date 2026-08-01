# Getting started

```bash
pip install shedding-hub
```

## Load a dataset

Every dataset lives in the repository as validated YAML and loads by its
`[author][year][word]` identifier.

```python
import shedding_hub as sh

data = sh.load_dataset("woelfel2020virological")
data["dataset_id"]
```

## Summarise it

`calc_shedding_summary` turns the raw measurements into one row per
participant — first/last positive time, duration, peak value, clearance.

```python
summary = sh.calc_shedding_summary(data, specimen="sputum")
summary.head()
```

## Plot it

```python
fig = sh.plot_time_course(data, specimen="sputum")
fig.savefig("time_course.png")
```

That is the whole loop: load, summarise, plot. The [reference](reference/datasets.md)
carries every argument these three take, and the [tutorial](tutorial.md) goes
on to simulate shedding for synthetic individuals rather than only describing
the studies already collected.
