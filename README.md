# Shedding Hub [![Shedding Hub](https://github.com/shedding-hub/shedding-hub/actions/workflows/build.yaml/badge.svg)](https://github.com/shedding-hub/shedding-hub/actions/workflows/build.yaml) [![DOI](https://zenodo.org/badge/836912278.svg)](https://doi.org/10.5281/zenodo.15052772)

The Shedding Hub collates data and statistical models for biomarker shedding (such as viral RNA or drug metabolites) in different human specimen (such as stool or sputum samples). Developing wastewater-based epidemiology into a quantitative, reliable epidemiological monitoring tool motivates the project.

Datasets are extracted from appendices, figures, and supplementary materials of peer-reviewed studies. Each dataset is stored as a [`.yaml`](https://en.wikipedia.org/wiki/YAML) file and validated against our [data schema](https://github.com/shedding-hub/shedding-hub/blob/main/data/.schema.yaml) to verify its integrity.

## 📊 Getting the Data

You can obtain the data by [downloading it from GitHub](https://github.com/shedding-hub/shedding-hub/tree/main/data). We also provide a [convenient Python package](http://pypi.org/project/shedding-hub/) so you can download the most recent data directly in your code or obtain a specific version of the data for reproducible analysis. Install the package by running `pip install shedding-hub` from the command line; it requires Python 3.10 or newer. The example below downloads the [data from Wölfel et al. (2020)](https://shedding-hub.github.io/datasets/woelfel2020virological.html) as of the commit [`259ca0d`](https://github.com/shedding-hub/shedding-hub/commit/259ca0d).

```python
>>> import shedding_hub as sh

>>> sh.load_dataset('woelfel2020virological', ref='259ca0d')
{'title': 'Virological assessment of hospitalized patients with COVID-2019',
 'doi': '10.1038/s41586-020-2196-x',
 ...}

```

You can also check whether a paper is already in the dataset collection using the `check_dataset` function.

```python
>>> sh.check_dataset(doi='10.1038/s41586-020-2196-x')
True

>>> sh.check_dataset(title='Virological assessment of hospitalized patients with COVID-2019')
True

```

## 📈 Analyzing the Data

The package provides statistical summaries and visualization tools to analyze shedding patterns across studies.

### Statistical Summaries

Calculate per-participant shedding statistics including duration, peak values, and clearance status.

```python
>>> data = sh.load_dataset('woelfel2020virological', ref='259ca0d')
>>> summary = sh.calc_shedding_summary(data, specimen='sputum')
>>> list(summary.columns)  # doctest: +NORMALIZE_WHITESPACE
['participant_id', 'biomarker', 'specimen', 'value_type', 'reference_event',
 'first_positive_time', 'last_positive_time', 'shedding_duration', 'peak_value',
 'peak_time', 'n_positive', 'n_negative', 'n_total', 'clearance_status', 'clearance_time']

```

Analyze detection rates over time with confidence intervals.

```python
>>> detection = sh.calc_detection_summary(data, specimen='sputum', time_bin_size=7)
>>> list(detection.columns)
['time', 'n_tested', 'n_positive', 'n_negative', 'proportion', 'ci_lower', 'ci_upper']

```

Compare shedding patterns across multiple datasets.

```python
>>> data1 = sh.load_dataset('woelfel2020virological', local='./data')
>>> data2 = sh.load_dataset('kimse2020viral', local='./data')
>>> comparison = sh.compare_datasets([data1, data2], specimen='sputum', value='concentration')
>>> list(comparison.columns)  # doctest: +NORMALIZE_WHITESPACE
['dataset_id', 'n_participants', 'n_measurements', 'pct_positive',
 'median_shedding_duration', 'iqr_shedding_duration', 'median_peak_value',
 'iqr_peak_value', 'median_peak_time', 'pct_cleared', 'median_clearance_time']

```

### Simulating Shedding

Simulate shedding trajectories for synthetic infected individuals — intended for
agent-based models of wastewater surveillance. Browse the catalog of fitted
estimates, pick one study or an ensemble across studies, then simulate.

Full documentation — a tutorial, the modeling methods, and a generated API
reference covering every public name with a worked example — is at
**[shedding-hub.readthedocs.io](https://shedding-hub.readthedocs.io/)**.

[`examples/simulating-shedding.md`](https://github.com/shedding-hub/shedding-hub/blob/main/examples/simulating-shedding.md) walks through the whole workflow, including how
to override the default choice and what the estimates do *not* support. It refits
nothing and runs in seconds. It is a [jupytext](https://jupytext.readthedocs.io)
notebook, like the extraction scripts, so open it directly in Jupyter or convert
it first:

```bash
jupytext --to ipynb examples/simulating-shedding.md
jupyter lab examples/simulating-shedding.ipynb
```

*The examples in this section run in CI as doctests against the shipped catalog
(`shedding_hub/data/shedding_catalog.yaml`), deliberately not skipped, so that
documentation drift fails loudly: if a future catalog rebuild gates out the
`woelfel2020virological` stool gamma fit used below, update the
`dataset_id`/`analyte`/`model`/ensemble filters to a fit that still exists
rather than re-adding `+SKIP`.*

```python
>>> import numpy as np
>>> import shedding_hub as sh
>>> catalog = sh.load_shedding_catalog()
>>> catalog.table[['dataset_id', 'specimen', 'model', 'peak_day']].head()  # doctest: +SKIP
>>> fit = catalog.select(
...     dataset_id='woelfel2020virological', analyte='stool', model='gamma'
... )
>>> traj = sh.simulate_shedding(
...     fit, n_individuals=100, times=np.arange(0, 30), seed=42
... )
>>> list(traj.columns)
['individual_id', 'time', 'log10_value', 'value', 'detected', 'source_dataset_id']

```

Picking a fit by hand means naming six keys — biomarker, specimen, reference
event, unit, value type and model — and those keys cut the catalog into 82
groups, 71 of which hold a single study. To see the choice, and to have it
made for you:

```python
>>> import shedding_hub as sh
>>> options = sh.shedding_options(biomarker='SARS-CoV-2', specimen='stool')
>>> list(options.columns)
['biomarker', 'specimen', 'reference_event', 'event_class', 'unit', 'value_type', 'n_unit_studies', 'model', 'n_studies', 'n_subjects', 'n_measurements', 'rank']
>>> source = sh.shedding_for('SARS-CoV-2', 'stool')
>>> source.selection.picked['event_class']
'landmark'

```

`shedding_for` takes rank 1 from `shedding_options`, preferring a reference event
that can be placed on an infection timeline, then the unit most studies report
for that biomarker and specimen, then a model that resolves the rise, then the
weight of evidence. Pass `model=`, `unit=` or `reference_event=` to pin any of
them, and read `source.selection` for what was chosen and what it beat.

`source.selection` is a `Selection`: `picked` (the winning group's keys and
counts), `passed_over` (the rest of the ranked table), `reason` (the rule that
decided it) and `analytes` (which analyte was taken from each study offering
more than one). `str(selection)` summarises all of it in one line:

```python
>>> print(sh.shedding_for('SARS-CoV-2', 'stool').selection)
picked symptom onset / gc/mL / gamma (2 study/studies, 16 subjects); ...

```

Three models are available. `exponential` is a pure decay from the reference
event. `gamma` rises and falls after it. `gamma_shifted` is the same rise and
fall with a fitted onset `t0`, so its support starts when shedding started
rather than at the reference event:

```
c(t) = c0 * (t - t0)**b0 * exp(-a0 * (t - t0)),   t > t0
```

Both rise-and-fall models are only fitted where a rise was actually observed —
at least half of a study's subjects must have their peak reading later than
their first sample, since otherwise the rise parameter is unidentifiable.
`catalog.skipped` records why anything is missing.

`gamma_shifted` exists because `gamma` is undefined at `t <= 0` and therefore
discards every reading there, including **26,023 detected measurements at
exactly the reference event**. It is offered only for analytes that have a
detected reading at or before their reference event; without one, `t0` has
nothing to locate and merely absorbs curve shape.

**Do not compare `gamma` and `gamma_shifted` by AIC.** Where `gamma_shifted` is
admitted it is fitted to *more observations* than `gamma` — 2072 against 1679
for `kissler2021viral` — and AIC is only comparable across models fitted to the
same data. The choice between them is made by data availability, which is what
the gate above encodes, not by fit statistic. Comparing `exponential` against
either is likewise only meaningful when `n_measurements` matches.

Pass `incubation_period` to express times as days since infection rather than
days since the study's reference event:

```python
>>> traj = sh.simulate_shedding(
...     fit, n_individuals=100, times=np.arange(0, 30),
...     incubation_period=5.0, seed=42
... )
>>> traj.attrs['time_origin']
'infection'

```

The gamma curve is undefined at or before the reference event, so any row
whose time falls there comes back as `NaN` with `detected=False` — including,
when `incubation_period` shifts the timeline, every early `times` entry that
still falls within the incubation window. `pandas` skips `NaN` when summing,
so aggregating simulated load across a cohort is safe by default, but account
for it if you do your own arithmetic on `log10_value` or `value`.

To pool evidence across studies, build an ensemble. Each simulated individual is
drawn from one contributing study, so between-study variation is preserved:

```python
>>> ensemble = catalog.ensemble(
...     biomarker='SARS-CoV-2', specimen='stool',
...     reference_event='symptom onset', unit='gc/mL', model='gamma',
... )
>>> traj = sh.simulate_shedding(
...     ensemble, n_individuals=1000, times=np.arange(0, 30), seed=42
... )

```

Estimates come from a censored maximum-likelihood fit, so `negative`
measurements inform the fit rather than being discarded. Because the two-stage
fit does not shrink individual estimates toward the population mean, simulated
cohorts are somewhat more dispersed than reality.

That over-dispersion matters when a few agents can dominate a total. Pass
`dispersion` below 1 to scale the between-subject covariance by `dispersion ** 2`,
narrowing the cohort's spread while leaving its centre and correlation structure
alone:

```python
>>> traj = sh.simulate_shedding(
...     fit, n_individuals=1000, times=np.arange(0, 30),
...     dispersion=0.7, seed=42,
... )

```

The default of `1.0` simulates the fitted population exactly as estimated. There
is no automatic way to choose a lower value — it is a judgement about how much of
the fitted spread is real estimation noise rather than genuine heterogeneity. What
makes shrinkage the only direction offered is that the two-stage bias runs one
way: the fitted spread is too wide, never too narrow.

For how the estimates are produced, what they mean, and where they should not
be trusted, see [docs/modeling-methods.md](https://github.com/shedding-hub/shedding-hub/blob/main/docs/modeling-methods.md). `make
parameters` writes the fitted parameters for every dataset to
`docs/shedding_parameters.json` (reusable: each record carries the population
mean and covariance, so it can be simulated from without refitting) and
`docs/shedding_parameters.csv` (flat, one row per fit).

### Visualization

Plot individual shedding trajectories over time.

```python
>>> fig = sh.plot_time_course(data, specimen='sputum')

```

Visualize aggregate shedding patterns with mean or median trajectories and confidence bands.

```python
>>> fig = sh.plot_mean_trajectory(data, specimen='sputum', central_tendency='median')

```

Generate heatmaps to compare shedding across participants.

```python
>>> fig = sh.plot_shedding_heatmap(data, specimen='sputum')

```

Create Kaplan-Meier clearance curves for survival analysis.

```python
>>> from shedding_hub.viz import plot_clearance_curve
>>> fig = plot_clearance_curve(data, specimen='sputum')

```

Compare the fitted curves themselves across studies. Each panel holds one
`(biomarker, specimen, unit, reference_event)` group, since curves disagreeing on
either axis cannot honestly be overlaid; colour identifies the study and
linestyle the model. The stretch of each curve before a study's median first
observation is faded, marking where the curve is functional form rather than
measurement — for half the catalog's gamma fits, that is the entire rise phase.

```python
>>> fig = sh.plot_catalog_fits(
...     catalog, biomarker='SARS-CoV-2', unit='gc/mL',
...     reference_event='symptom onset',
... )

```

To judge a single fit rather than compare several, plot it against the data
behind it. The points come from the fitter's own view of the dataset, so the page
shows exactly what the fit saw — the same censoring limit, the same excluded
subjects — and censored readings are drawn on the limit rather than dropped. The
estimated parameters and the fit's context sit in the legend.

```python
>>> fig = sh.plot_fit_diagnostic(fit, data1)

```

Behind the observations it shades the central 5–95% of a simulated cohort drawn
from the fitted population, because the median individual alone says nothing about
whether the *spread* is right — which is most of what separates a usable fit from
an unusable one. Pass `dispersion` to see the narrowed cohort, or
`show_band=False` to omit it.

Measurements the fitter discarded are drawn too, marked as excluded rather than as
data the curve should explain. The gamma model is undefined at `t <= 0`, so
readings there are dropped — 391 of them for `kissler2021viral` — and a page that
drew nothing would imply the study never sampled before its reference event.

`make review` renders every fit in the catalog this way into a single
`shedding_catalog_review.pdf`, one page each. `make review_range` writes a second
PDF shading the full range of the simulated cohort rather than its central 90%,
with the 95% interval drawn inside it as two dashed lines and the y axis widened
to fit — what each fit considers *possible* rather than typical. A range is also a
property of how many individuals were drawn, which is why each page names its
draw count.

To judge how much the over-extrapolation gate matters, rebuild the catalog at a
different threshold and render it:

```bash
python scripts/build_shedding_catalog.py --max-peak-above-observed 2 \
    --output shedding_catalog_gate2.yaml
python scripts/build_catalog_review.py --catalog shedding_catalog_gate2.yaml \
    --output shedding_catalog_review_gate2.pdf
```

## 🤝 Contributing

Thank you for contributing your data to the Shedding Hub and supporting wastewater-based epidemiology! If you hit a bump along the road, [create a new issue](https://github.com/shedding-hub/shedding-hub/issues/new) and we'll sort it out together.

We use [pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) to add and update data, allowing for review and quality assurance. Learn more about the general workflow [here](https://docs.github.com/en/get-started/using-github/github-flow). To contribute your data, follow these easy steps (if you're already familiar with pull requests, steps 2 and 3 are for you):

1. Create a [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) of the Shedding Hub repository by clicking [here](https://github.com/shedding-hub/shedding-hub/fork) and [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the fork to your computer. You only have to do this once.
2. Create a new `my_cool_study/my_cool_study.yaml` file in the [`data`](https://github.com/shedding-hub/shedding-hub/tree/main/data) directory and populate it with your data. See [here](https://github.com/shedding-hub/shedding-hub/blob/main/data/woelfel2020virological/woelfel2020virological.yaml) for a comprehensive example from [Wölfel et al. (2020)](https://www.nature.com/articles/s41586-020-2196-x). A minimal example for studies with a single analyte (e.g., SARS-CoV-2 RNA concentration in stool samples) is available [here](https://github.com/shedding-hub/shedding-hub/blob/main/tests/examples/valid_single_analyte.yaml), and a minimal example for studies with multiple analytes (e.g., crAssphage RNA concentration in stool samples and caffeine metabolites in urine) is available [here](https://github.com/shedding-hub/shedding-hub/blob/main/tests/examples/valid_multiple_analytes.yaml).
3. Optionally, if you have a recent version of [Python](https://www.python.org) installed, you can validate your data to ensure it has the right structure before contributing it to the Shedding Hub.
    - Run `pip install -r requirements.txt` from the command line to install all the Python packages you need.
    - Run `pytest` from the command line to validate all datasets, including the one you just created.
4. Create a new [branch](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches) by running `git checkout -b my_cool_study`. Branches let you isolate changes you are making to the data, e.g., if you're simultaneously working on adding multiple studies–much appreciated! You should create a new branch from the `main` branch for each dataset you contribute; see [here](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) for more information.
5. Add your changes by running `git add data/my_cool_study/my_cool_study.yaml` and commit them by running `git commit -m "Add data from Someone et al. (20xx)."`. Feel free to pick another commit message if you prefer.
6. Push the dataset to your fork by running `git push origin my_cool_study`. This will send the data to GitHub, and the output of the command will include a line `Create a pull reuqest for 'my_cool_study' on GitHub by visiting: https://github.com/[your-username]/shedding-hub/pull/new/my_cool_study`. Click on the link and follow the next steps to create a new pull request.

Congratulations, you've just created your first pull request to contribute a new dataset! We'll now [review the changes](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/about-pull-request-reviews) you've made to make sure everything looks good. Once any questions have been resolved, we'll [merge your changes](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/merging-a-pull-request) into the repository. You've just contributed your first dataset to help make wastewater-based epidemiology a more quantitative public health monitoring tool–thank you!

## 🔄 What happens after a dataset is merged

Each dataset gets a page on [the website](https://shedding-hub.github.io/), and each of its analytes gets a figure there: the fitted shedding curve where a model could be fitted, and the measurements alone where none could. Those figures live in [`figures/`](https://github.com/shedding-hub/shedding-hub/tree/main/figures) in this repository, because the website has no Python to generate them with — it copies them out of the archive it already downloads for the dataset YAML.

Merging a dataset therefore sets off two things at once, and **one step in the middle is manual**:

```
merge a dataset PR into main
        │
        ├── github-pages.yaml  (any push to main)
        │      └── tells the website to rebuild, immediately
        │          → the new dataset's page goes live, with no figure yet
        │
        └── refresh-figures.yaml  (pushes touching data/**)
               └── ~45 min: rebuilds both gate-2 catalogs, re-renders every
                   figure, and opens a "figures: refresh the dataset pages" PR
                      │
                      └── a maintainer merges it        ← the manual step
                             │
                             └── the website rebuilds again, now with figures
```

A page with no figure is expected in the meantime and is not a fault: the layout omits the figure when the dataset has no entry yet, so the page is complete apart from the illustration.

**Reviewing the figures pull request.** Nearly every figure will show as changed, but not all of those changes are cosmetic, and the difference between the two matters.

Rendering itself is deterministic: measured across a full refresh, all 88 observations-only figures came back **pixel-for-pixel identical**, differing only in PNG encoding bytes. Every visible difference comes from *refitting*. Of the 194 fitted figures, 192 differed, by a median of 7% of their pixels.

Most of those are slight, but the tail is not. The largest in that refresh, `mijatovicrustempasic2017shedding / rotavirus vaccine_stool / exponential`, moved its half-life from **2.28 days to 1.97 days** — 14% — while its AIC changed by 0.1. Both fits explain the data equally well; the likelihood is simply flat along that direction, so the optimizer can settle anywhere along it. That is a property of the study's data rather than of the machine, and it is a caution worth carrying into any estimate drawn from a weakly-identified analyte.

So read the diff at two levels:

- **Added or removed figures** are the clearest signal: an analyte gained or lost a converged fit, changing what the catalog supports.
- **Among changed figures, read the legend numbers rather than the picture.** Parameters that shift by a few percent mark analytes whose fits are weakly identified. Note that `DRIFT_RTOL` in `tests/test_parameter_export.py` does *not* cover this — it is the noise floor for re-exporting one catalog, and its own comment says a genuine refit moves values by percent rather than by `1e-4`.

**Regenerating by hand.** `refresh-figures.yaml` only watches `data/**`, so a change to the fitting or plotting code does not raise a pull request on its own. Trigger one deliberately with `gh workflow run "Refresh Dataset Figures"`, or rebuild locally:

```bash
make catalog_gate2      # concentration fits, 2 log10 extrapolation gate
make catalog_ct_gate2   # cycle-threshold fits, 2 cycles — a stricter gate, see docs/modeling-methods.md
make figures            # 282 figures plus figures/index.json
```
