# Catalog fit diagnostics — design

**Date:** 2026-07-29
**Status:** implemented 2026-07-29, as designed. Extended the same day with a
simulated-cohort band — see
[2026-07-29-simulation-dispersion-design.md](2026-07-29-simulation-dispersion-design.md),
which added `show_band`, `dispersion`, `band_quantiles` and `n_simulated` to
`plot_fit_diagnostic`. The band answers the question the median individual cannot:
whether the *spread* of the fitted population covers the data.

## Problem

`plot_catalog_fits` draws every fitted curve in the catalog, but it draws them
alone. A curve with no data behind it cannot be judged: it shows what the model
says, never whether the model agrees with what the study measured.

Reviewing a fit means seeing three things at once — the observations, the curve,
and the parameter estimates — for each of the 83 fits individually. Nothing in
the package does that. `plot_time_course` draws a study's raw measurements but
knows nothing about fits; `plot_simulated_shedding` can overlay observations, but
on a simulated cohort rather than on the fit's own median individual, and only
for one fit at a time by hand.

The goal is one page per fit, all 83 in a single document to page through.

## What the points must be

Points come from `prepare_observations(dataset, analyte, model)` — the fitter's
own data preparation — not from the dataset's raw measurements.

This is the decision the whole design rests on. `prepare_observations` resolves
the censoring limit the fit actually used, drops subjects with too few usable
readings, and drops subjects with no positive reading at all. Plotting raw
measurements instead would put points on the page that the curve was never asked
to explain, and a reviewer would spend their attention on a disagreement that
does not exist.

It has a consequence worth stating: `min_observations` defaults to the number of
per-subject parameters, 2 for exponential and 3 for gamma, so a subject with
exactly two readings is retained on the exponential page and absent from the
gamma page for the same analyte. The two pages genuinely fit different data, and
seeing that is the point.

## Rendering

One fit per page, titled with `dataset_id / analyte / model`.

**Observed positives** are scattered at their measured times and values. Each
subject's points are joined by a faint line, but only when the fit retains at
most 40 subjects: the largest study has 455, where the lines would be a hairball
that hides the very scatter they are meant to organize.

**Censored observations are drawn, as open downward triangles on the censoring
limit.** They cannot be omitted. Across the catalog a median of 40% of
measurements are censored, 25 of 83 fits are above 50%, and one reaches 91%;
a page showing only the positives would make a mostly-undetected analyte look
like a clean decay through a handful of points. They sit *on* the limit and point
down because that is all the assay established — a value somewhere below it.

**The censoring limit** is a dotted horizontal line, labelled with its value.

**The fitted curve** is the median individual, `log10_concentration(model,
fit.median_params, t)` — the same object `plot_catalog_fits` draws. The stretch
before `median_first_observed_day` is faded to `alpha=0.35`, reusing that
function's convention so the two views agree about what is extrapolation.

**Both axes come from the data.** Unlike `plot_catalog_fits`, which has no
observations and must derive a horizon from the decay rate, here the observations
exist and define the range: x from 0 to just past the last observed day, y from
just below the censoring limit to just above the highest observation.

## The legend carries the estimates

Two blocks, because the raw parameters and their interpretable transforms answer
different questions and the repository has an established position that the raw
ones are not self-explanatory:

- **Estimated parameters** — `a0`, `b0`, `c0`, whichever the model has.
- **Derived summaries** — peak day, peak log10, half-life, sigma.
- **Fit context** — subjects, measurements, % censored, AIC, and whether the
  optimizer converged.

Context sits beside the estimates because a parameter value cannot be judged
without it: 6.8 log10 at peak means one thing from 52 subjects and another from
3, and `converged=False` disqualifies the rest of the block.

## API

Added to `shedding_hub/viz.py`, exported from `shedding_hub`:

```python
def plot_fit_diagnostic(
    fit: SheddingFit,
    dataset: dict,
    *,
    figsize: tuple[float, float] = (9, 6),
    max_subject_lines: int = 40,
) -> Figure
```

`dataset` is a loaded dataset dictionary, from `load_dataset`. The fit already
knows its `analyte` and `model`, so nothing else need be passed.

Passing a dataset that does not contain the fit's analyte raises `ValueError`
naming both, rather than drawing a curve over an empty panel — the same
reasoning as `plot_catalog_fits` refusing an empty filter match.

## The document

`scripts/build_catalog_review.py` loads the shipped catalog, loads each of the 29
datasets once, and renders all 83 pages into `shedding_catalog_review.pdf` with
`matplotlib.backends.backend_pdf.PdfPages`, ordered by dataset, analyte, then
model, so a study's two models face each other. A `review` target is added to the
Makefile beside the existing `catalog` target.

The PDF is not committed. It is a regenerable binary that would be rewritten on
every catalog rebuild, and `make review` reproduces it; it is added to
`.gitignore` instead.

## Testing

Built with TDD, in `tests/test_viz.py`, against a real fit of the synthetic
dataset the `make_synthetic_dataset` fixture builds — a genuine fit/data pair
rather than a stub, since the point of this plot is the relationship between the
two.

- Observed positives appear at their measured coordinates.
- Censored observations are drawn, on the censoring limit, with a distinct
  marker from the positives.
- A subject excluded by the fitter contributes no points, which is the
  observations-come-from-the-fitter rule.
- Subject-joining lines appear below the threshold and are suppressed above it.
- The legend text carries each estimated parameter name and the derived
  summaries.
- The faded segment appears when `median_first_observed_day` falls inside the
  observed range.
- A dataset lacking the fit's analyte raises `ValueError`.
- The figure has one axis and is closed in the pyplot state, matching the
  convention in `shedding_peak.py` and `shedding_simulate.py`.

## Out of scope

**No per-subject fitted curves.** The shipped catalog does not serialize
`subject_params` — a loaded fit carries `None` — so drawing each subject's own
fitted curve would mean refitting all 83 analytes at review time. The population
median curve against the observed scatter is the review this delivers.

**No goodness-of-fit statistics.** Residual plots, QQ plots and posterior
predictive checks are a larger piece of work with their own design questions.
This page is for looking, and `aic` is already on it.

**No changes to the fitter.** This design consumes `prepare_observations` and
`SheddingFit` exactly as they are.
