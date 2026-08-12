---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.7
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- --8<-- [start:body] -->
# Simulating shedding for synthetic individuals

A short tour of the workflow an agent-based model needs: **browse the
fitted estimates → pick one → simulate individuals**.

Nothing here refits anything. The catalog ships precomputed, so every
cell below runs in seconds.

Run top to bottom with *Run All*.

```python
import warnings

import numpy as np
import shedding_hub as sh

catalog = sh.load_shedding_catalog()
print(f"{len(catalog.table)} fits over "
      f"{catalog.table.dataset_id.nunique()} datasets, "
      f"models {sorted(catalog.table.model.unique())}")
```

## 1. What is in the catalog

One row per fit, summarised by its median individual. `peak_day` and
`half_life_days` are the interpretable summaries; the full population
mean and covariance live on the fit objects.

```python
catalog.table[
    ["dataset_id", "biomarker", "specimen", "reference_event",
     "unit", "model", "n_subjects", "peak_day", "half_life_days"]
].head(8)
```

Analytes the fitter **refused** are recorded too, with a reason — so a
study missing from the table is a decision, not a silent gap.

```python
catalog.skipped["reason"].value_counts()
```

## 2. See the choice

Picking a fit by hand means agreeing on five keys — biomarker,
specimen, reference event, unit and model — and those keys cut the
catalog into 82 groups, 71 of which hold a single study. So the hard
part is not *combining* estimates; it is choosing among ones that are
not comparable.

`shedding_options` lists every group that can actually be built, best
first. `rank` 1 is what `shedding_for` will return.

```python
options = sh.shedding_options("SARS-CoV-2", "stool", catalog=catalog)
options
```

The ordering encodes a judgement, in strict precedence:

1. **a reference event that can be placed on an infection timeline** —
   `symptom onset` beats `enrollment`, because the date someone joined
   a study has no fixed relation to when they were infected;
2. **the unit most studies report** — `gc/mL` and `gc/dry gram` measure
   different things, so this is settled on weight of evidence rather
   than as a side effect of which unit happened to support a richer
   model;
3. **a model that resolves the rise** — `gamma_shifted`, then `gamma`,
   then `exponential`;
4. **weight of evidence** — studies, then subjects, then measurements.

Rule 4 coming last is deliberate and is the one real cost: within a
settled clock and unit, a 2-study `gamma` is preferred to a 3-study
`exponential`, because the exponential pins peak shedding to the
reference event by construction — a structural misstatement rather
than a matter of precision.


## 3. Make the choice

`shedding_for` takes rank 1 and builds it. It always returns a
`SheddingEnsemble` — single-component when only one study matched — so
you get one type and one code path however many studies backed it.

```python
source = sh.shedding_for("SARS-CoV-2", "stool", catalog=catalog)

print(source.selection)
print("\ncomponents:", len(source.fits))
print("study -> analyte:", source.selection.analytes)
print("weights:", source.weights.round(3).tolist())
```

`selection` is data, not just a printed line, so a model run can record
*why* it used what it used. `passed_over` is the rest of the ranked
table.

```python
print("picked:", source.selection.picked)
print("\nreason:", source.selection.reason)
source.selection.passed_over.head(3)
```

## 4. Simulate a cohort

Each individual gets its own parameters drawn from the fitted
population, so the spread across agents is between-person variation,
not noise around one average curve.

`method="mixture"` (the default) draws a study per agent, which keeps
genuine disagreement between studies visible instead of averaging it
away — note `source_dataset_id` in the output.

```python
traj = sh.simulate_shedding(
    source,
    n_individuals=500,
    times=np.arange(1, 31),
    seed=42,
)
traj.head()
```

```python
print("rows:", len(traj))
print("attrs:", {k: traj.attrs[k] for k in
                 ("time_origin", "reference_event_class", "unit", "model")})
print("\nagents per source study:")
print(traj.groupby("source_dataset_id").individual_id.nunique())
```

```python
fig = sh.plot_simulated_shedding(traj, source=source)
fig  # ending the cell on the figure renders it
```

The band is the central 90% of the simulated cohort and the line its
median. For a cohort-level quantity — total load into a sewer, say —
sum the individual trajectories rather than scaling the median, which
is not the mean trajectory.

```python
# What an agent-based model actually consumes: one row per agent per day.
wide = traj.pivot(index="individual_id", columns="time", values="value")
print("shape (agents x days):", wide.shape)

daily_total = wide.sum(axis=0)
print("\ncohort total load, gc/mL summed over 500 agents:")
print(daily_total.head(8).map(lambda v: f"{v:,.0f}").to_string())
```

## 5. Time origin — the part to get right

An agent is infected, not "symptom onset-ed". `incubation_period`
shifts the clock back from the study's reference event to infection.

```python
infection_time = sh.simulate_shedding(
    source,
    n_individuals=500,
    times=np.arange(0, 31),
    incubation_period=5.0,   # days from infection to symptom onset
    seed=42,
)
print("time_origin:", infection_time.attrs["time_origin"])
```

### The pre-symptomatic window is not free

The `gamma` model is undefined before its own reference event, so once
the clock starts at infection, every day before symptom onset comes
back `NaN`. That is the model declining to extrapolate, and it is the
epidemiologically interesting window.

```python
defined = (
    infection_time.groupby("time").log10_value
    .apply(lambda s: float(np.isfinite(s).mean()))
)
print("fraction of agents with a defined concentration, by day since infection:")
print(defined.head(9).to_string())
```

`gamma_shifted` exists to fit an onset rather than assume one — but it
does not conjure data that was never collected. For SARS-CoV-2 in
stool the fitted onset is **after** symptom onset, because stool
sampling in these studies began late:

```python
shifted = sh.shedding_for(
    "SARS-CoV-2", "stool", catalog=catalog, model="gamma_shifted"
)
print(shifted.selection)
print("\nfitted onset t0, in days from symptom onset:",
      round(float(shifted.fits[0].median_params[3]), 2))
```

So for this biomarker and specimen there is **no fit that supports a
pre-symptomatic stool signal**, and the honest options are to simulate
on the study's own clock (leave `incubation_period` unset), or to
start `times` at the reference event. Fabricating the early window
would be inventing data the studies never collected.


### Reference events that cannot be anchored

Only `symptom onset` earns a genuine `"infection"` origin. Shifting a
fit anchored to an administrative date — when someone enrolled, or was
confirmed — asserts a relationship that does not exist, so it warns and
records what it actually did.

```python
administrative = sh.shedding_for(
    "SARS-CoV-2", "stool", catalog=catalog, reference_event="enrollment"
)

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    shifted_admin = sh.simulate_shedding(
        administrative, n_individuals=20, times=np.arange(1, 10),
        incubation_period=5.0, seed=1,
    )

print("time_origin:", shifted_admin.attrs["time_origin"])
print("class      :", shifted_admin.attrs["reference_event_class"])
print("\nwarning:", caught[-1].message)
```

## 6. Not taking the default

The ranking is a documented opinion, not a verdict. There is a ladder
of control, from nudging the pick to bypassing it entirely.

### 6.1 Pin a key

Pass the key you want fixed and the ranking applies only to what is
left. **Any fit attribute works** — not just the three that are ranked.

```python
choices = [
    dict(),                                        # the default
    dict(model="exponential"),
    dict(unit="gc/dry gram"),
    dict(reference_event="enrollment"),
    dict(unit="gc/dry gram", model="gamma"),       # pin two at once
    dict(dataset_id="woelfel2020virological"),     # a single named study
    dict(gene_target="N1"),
]
for choice in choices:
    picked = sh.shedding_for("SARS-CoV-2", "stool", catalog=catalog, **choice)
    label = str(choice) if choice else "(default)"
    print(f"{label:<48} -> {picked.selection.picked['unit']:<12} "
          f"{picked.model:<14} {len(picked.fits)} study(ies)")
```

Check `shedding_options` first to see what is available — every row it
lists can actually be built, so anything pinned from that table will
work. A combination that does not exist raises rather than quietly
substituting something near it:

```python
try:
    sh.shedding_for(
        "SARS-CoV-2", "stool", catalog=catalog,
        unit="gc/mL", model="gamma_shifted", reference_event="symptom onset",
    )
except ValueError as error:
    print("ValueError:", error)
```

### 6.2 Change how the studies are combined

Independent of *which* group you get. `weights` defaults to
`"n_subjects"`, so larger studies count for more; `method` defaults to
`"mixture"`, which draws a study per agent and keeps genuine
disagreement between studies visible. `"moment"` collapses the
components into a single Gaussian instead.

```python
for weights in ["n_subjects", "equal"]:
    src = sh.shedding_for(
        "SARS-CoV-2", "stool", catalog=catalog,
        model="exponential", weights=weights,
    )
    print(f"weights={weights:<12} -> {src.weights.round(3).tolist()}")

src = sh.shedding_for("SARS-CoV-2", "stool", catalog=catalog, model="exponential")
print("\ncomponent studies :", [f.dataset_id for f in src.fits])
print("component subjects:", [f.n_subjects for f in src.fits])
```

An explicit array is accepted too — but components are ordered by
`dataset_id`, so read `ensemble.components` before relying on the
positions.


### 6.3 Restrict which studies contribute

**`dataset_ids` (plural) is not a `shedding_for` key.** It would be
read as an attribute filter, match nothing, and raise. To hand-pick a
subset of studies, go through the catalog instead:

```python
chosen = catalog.ensemble(
    biomarker="SARS-CoV-2", specimen="stool",
    reference_event="symptom onset", unit="gc/mL", model="exponential",
    dataset_ids=["woelfel2020virological", "zuo2020alterations"],
)
print("components:", [f.dataset_id for f in chosen.fits])
print("weights   :", chosen.weights.round(3).tolist())
```

That is deliberately strict as well: naming a study that is not in the
group raises, rather than silently shrinking the ensemble to whatever
happened to match.

```python
try:
    catalog.ensemble(
        biomarker="SARS-CoV-2", specimen="stool",
        reference_event="symptom onset", unit="gc/mL", model="exponential",
        dataset_ids=["not_a_study"],
    )
except ValueError as error:
    print("ValueError:", str(error)[:160], "...")
```

### 6.4 Bypass the picker entirely

`catalog.select` returns one fit, which `simulate_shedding` accepts
directly. `make_ensemble` builds one from any fits at all — including
ones you fitted yourself on data that is not in this repository.

```python
one = catalog.select(
    dataset_id="woelfel2020virological", analyte="stool", model="gamma"
)
print(f"{one.dataset_id} {one.model}: peak day {one.peak_day:.2f}, "
      f"half-life {one.half_life_days:.2f}, n={one.n_subjects}")

by_hand = sh.make_ensemble(
    [one, catalog.select(
        dataset_id="lui2020viral", analyte="2019-nCoV_N1", model="gamma"
    )],
    weights="equal",
)
print("hand-built components:", [f.dataset_id for f in by_hand.fits])
```

### Whatever you choose, it stays auditable

`selection.reason` names the rule that decided and `passed_over` is
everything it beat, so an override is a recorded choice rather than
just a different answer. Note that a hand-built ensemble has no
`selection` — nothing chose it but you.

```python
override = sh.shedding_for("SARS-CoV-2", "stool", catalog=catalog,
                           model="exponential")
print("reason      :", override.selection.reason)
print("it beat     :", len(override.selection.passed_over), "other group(s)")
print("hand-built  :", by_hand.selection)
```

## 7. Using an estimate without this package

`docs/shedding_parameters.json` carries, per fit, the population mean
and covariance, the measurement-error SD and the censoring limit —
everything `simulate_shedding` needs. So an estimate can be reused from
any language, without installing this package and without refitting.

```python
import json, pathlib

path = pathlib.Path(sh.__file__).parent.parent / "docs" / "shedding_parameters.json"
records = json.loads(path.read_text())
print(records["n_fits"], "records over", records["n_datasets"], "datasets")

example = next(
    r for r in records["fits"]
    if r["dataset_id"] == "woelfel2020virological" and r["model"] == "gamma"
)
print(json.dumps(
    {k: example[k] for k in
     ("dataset_id", "model", "parameters", "summary",
      "measurement_error_sd", "censoring_limit_log10")},
    indent=2,
))
```

## What an estimate does not support

Worth reading before using any of this in anger — the full account is
in `docs/modeling-methods.md`:

- **Point estimates only.** There is no posterior uncertainty on the
  parameters; the spread you see across agents is between-person
  variation, not uncertainty in the fit.
- **Two-stage fitting overstates between-subject variance**, so cohorts
  are somewhat wider than the studies were.
- **`exponential` puts the peak at the reference event by
  construction.** Read `peak_log10` beside `median_first_observed_day`:
  where sampling began late, the peak is a backward extrapolation to a
  time nobody was observed at.
- **Do not compare models by AIC.** They are fitted to different
  observation sets — `gamma` discards every reading at or before the
  reference event, `gamma_shifted` keeps the detected ones — so the
  numbers are not commensurable. Compare `n_measurements` first.
<!-- --8<-- [end:body] -->
