# Getting started

```bash
pip install shedding-hub
```

Python 3.10 or newer is required.

This page walks the whole arc: load a curated study, describe it, fit a
shedding curve to it, and simulate synthetic individuals from fitted
parameters. The [tutorial](tutorial.md) goes deeper on the last step.

## Load a dataset

Every dataset lives in the repository as validated YAML and loads by its
`[author][year][word]` identifier.

```python
import shedding_hub as sh

data = sh.load_dataset("woelfel2020virological")
data["dataset_id"]
```

Pin it to a commit or a pull request when an analysis has to stay
reproducible — `sh.load_dataset("woelfel2020virological", ref="259ca0d")`.

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

## Fit a shedding curve

Descriptive summaries stop at what was observed. To get a *curve* — a peak
height, a time to peak, a decay rate — fit one.

`fit_shedding_model` fits by censored maximum likelihood, so non-detects
contribute what they actually say ("below this limit") rather than being
dropped or replaced by a substituted value. Three models are available:
`exponential` for post-peak decay only, `gamma` for a rise and a decay, and
`gamma_shifted` when shedding starts some time after the reference event.

```python
fit = sh.fit_shedding_model(data, analyte="stool", model="gamma")

fit.converged      # did the optimiser actually finish
fit.peak_day       # days from the reference event to peak
fit.peak_log10     # peak height, log10 of the analyte's unit
fit.half_life_days # decay half-life after the peak
fit.n_subjects     # participants the fit is built on
```

The per-participant estimates are on `fit.subject_params`, one row each, with
a `degenerate` flag marking the ones whose parameters collapsed onto a bound
and are therefore excluded from the population summary:

```python
fit.subject_params.head()
```

The fitter refuses work it cannot do rather than returning a confident
number. Ask for `gamma` on a series where sampling began after the peak and it
raises `SheddingDataError` naming the problem — that study's `sputum` series
has only 44% of subjects showing a rise, against a required 50% — and tells
you to use `exponential` instead.

Check the fit against the data it came from. It takes the dataset too, so the
observations and the curve land on the same page:

```python
fig = sh.plot_fit_diagnostic(fit, data)
```

![plot_fit_diagnostic](images/plot_fit_diagnostic.png)

### Seeing what the fit implies, not just what it fitted

The median individual says nothing about spread, which is most of what
separates a usable fit from an unusable one. Widen the band to the full range
of a simulated cohort and the page stops describing the curve and starts
describing the population it would generate — which is what a simulation drawn
from this fit will actually produce:

```python
fig = sh.plot_fit_diagnostic(
    fit,
    data,
    band_quantiles=(0.0, 1.0),        # the whole simulated range, not its middle 90%
    band_inner_quantiles=(0.025, 0.975),  # with the central 95% dashed inside it
    band_sets_ylim=True,              # let that range widen the axis
    x_from_fitted=True,               # and keep discarded readings from stretching time
)
```

![plot_fit_diagnostic_range](images/plot_fit_diagnostic_range.png)

A range is not a fixed property of the population: it grows with
`n_simulated`, because drawing more individuals reaches further into the tails.
Read it as what the fit considers *possible*, and the dashed lines as where the
mass actually sits. This is the view the
[dataset pages](https://shedding-hub.github.io/) carry, fitted under a stricter
extrapolation gate — pass `max_peak_above_observed=2` to `fit_shedding_model`
to reproduce it.

## Fit cycle-threshold data

**78 of the repository's 216 analytes report cycle thresholds** rather than
concentrations — 31 datasets and 12,698 measurements. They are fitted by the
same function, with no conversion and no assumption about PCR efficiency:

```python
data = sh.load_dataset("wang2020fecal")
fit = sh.fit_shedding_model(data, analyte="stool_SARSCoV2_N", model="gamma")

fit.value_type     # 'ct'
fit.peak_day       # days to peak — comparable with a concentration fit
fit.ct_cutoff      # the Ct at or above which a reading counts as a non-detect
```

![plot_fit_diagnostic_ct](images/plot_fit_diagnostic_ct.png)

Because `Ct = α − β·log10 C` is affine, a gamma curve in concentration is a
gamma curve in Ct, and the analyte is fitted on `CT_REFERENCE − Ct`. The
unknown slope `β` multiplies both `a0` and `b0`, so it **cancels in `b0/a0`**:

```python
fit.comparable_with(concentration_fit)
# ('peak_day', 't0', 'rise_days')
```

Peak time, onset and rise duration therefore compare across the two scales.
Decay rate, half-life and peak height do not — they carry `β`, so a Ct fit's
`half_life_days` is *not* a viral half-life. Heights are reported as Ct rather
than log10, and the axis carries real Ct numbers.

Two deliberate limits: the shipped catalog is concentration-only, so Ct fits
need `fit_shedding_models(..., value_types=("ct",))` to build; and
`simulate_shedding` **refuses a Ct source**, because turning cycles back into
concentrations needs a standard curve the studies do not report. See
[modeling methods](modeling-methods.md) for how far the cross-scale comparison
actually gets you — the algebra is exact, the estimates less so.

## Use fits that are already computed

You do not have to fit anything. The package ships a catalog of **159
converged fits over 53 studies**, each one a censored ML fit with its sample
size, censoring rate and diagnostics attached.

```python
catalog = sh.load_shedding_catalog()
catalog.table.head()
```

`catalog.table` is a DataFrame with one row per fit — `dataset_id`,
`biomarker`, `specimen`, `unit`, `reference_event`, `model`, `n_subjects`,
`pct_censored`, `peak_day`, `peak_log10`, `half_life_days`, `aic` — so you can
filter it however you like:

```python
t = catalog.table
t[(t.biomarker == "SARS-CoV-2") & (t.specimen == "stool")]
```

## Choose a source

That filter returns several rows, because several studies measured the same
biomarker in the same specimen — in different units, against different
reference events, in different populations. Picking between them is a real
decision, so the package makes it explicit rather than silently choosing.

`shedding_options` ranks the candidates and shows its work:

```python
opts = sh.shedding_options("SARS-CoV-2", "stool")
opts[["reference_event", "event_class", "unit", "model", "n_studies", "n_subjects", "rank"]].head(4)
```

```text
reference_event event_class        unit         model  n_studies  n_subjects  rank
  symptom onset    landmark       gc/mL         gamma          2          16     1
  symptom onset    landmark       gc/mL   exponential          3          26     2
  symptom onset    landmark gc/dry gram gamma_shifted          1          29     3
  symptom onset    landmark gc/dry gram         gamma          1          30     4
```

Ranking prefers a reference event you can actually anchor a simulation to
(`symptom onset` over `enrollment`), then the unit most studies agree on, then
the model that resolves the most structure, then sample size.

`shedding_for` takes the top row and builds something you can simulate from:

```python
source = sh.shedding_for("SARS-CoV-2", "stool")

source.unit             # 'gc/mL'
source.reference_event  # 'symptom onset'
source.selection.reason # 'its model resolves the rise'
source.selection.passed_over
```

It defaults to a `mixture` ensemble — each simulated individual draws its
parameters from one contributing study, preserving between-study spread. Pass
`method="moment"` for a single pooled Gaussian instead. Section 6 of the
[tutorial](tutorial.md) covers overriding the choice entirely.

## Simulate a cohort

```python
sim = sh.simulate_shedding(source, n_individuals=500, times=range(0, 31))
sim.head()
```

One row per individual per time point — exactly what an agent-based or
wastewater model consumes:

| column | meaning |
| --- | --- |
| `individual_id` | synthetic individual, `0` to `n_individuals - 1` |
| `time` | days from the time origin |
| `log10_value` | simulated concentration, log10 of `source.unit` |
| `value` | the same on the natural scale |
| `detected` | whether it clears the assay's detection limit |
| `source_dataset_id` | which study supplied that individual's parameters |

Inter-individual variability is drawn from the fitted population covariance,
so the cohort spreads the way the studies did rather than tracing one average
curve. `sim.attrs` records the provenance — `time_origin`, `model`, `unit`,
`biomarker`, `specimen` — so a saved simulation stays interpretable.

```python
fig = sh.plot_simulated_shedding(sim)
```

!!! warning "Day 0 means the reference event, not infection"
    Times are relative to the source's `reference_event`. For
    `symptom onset`, day 0 is symptom onset — a simulated individual is
    already shedding before it. Pass `incubation_period=` to shift onto an
    infection time origin. The
    [tutorial](tutorial.md) works through why this matters and when it cannot
    be done at all.

## Where to go next

- **[Tutorial](tutorial.md)** — the simulation path end to end, including overriding the source selection.
- **[Modeling methods](modeling-methods.md)** — what the fitted estimates mean, and what they do not support. Read this before publishing anything built on them.
- **[Reference](reference/datasets.md)** — every argument of every public function, with a worked example.
