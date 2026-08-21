# Taming simulated shedding outliers — design

**Date:** 2026-07-29
**Status:** implemented 2026-07-29

## Problem

Simulating a cohort from a fitted population produced individuals no biology
supports. From `woelfel2020virological` stool exponential — 8 subjects whose own
day-0 levels span 4.0 to 9.0 log10 gc/mL — a 10,000-agent cohort reached
**10^19 gc/mL**, and the top 0.1% of agents carried essentially *all* of the
cohort's total shed load. For an agent-based model that sums load across a
population, a handful of impossible agents decide the answer.

Two independent causes, which compound.

## Cause 1: the exponential model's height coordinate (the larger one)

`POPULATION_COORDS["exponential"]` was `("log_a0", "log_c0")`. Since
`log10 c(0) = c0 / ln(10)`, modelling `log c0` as normal makes the *log10*
concentration lognormal, and therefore the concentration a **double exponential**
of the draw, `e^(e^θ)`. A 3.5-sigma draw moved the day-0 level from 6.1 to 18.8
log10.

The gamma model never had this problem: `to_population_coords` already gives it
`peak_log10` as a coordinate, so its log10 height is normal and its concentration
merely lognormal. The fix is to give the exponential model the same treatment —
`("log_a0", "peak_log10")`, the peak being at `t = 0`.

The original docstring justified leaving the exponential model in `log(params)`
on the grounds that its subjects form a compact cloud, which is true: the
*curved-ridge* argument that motivated the gamma reparameterization genuinely
does not apply to it. But that argument is about the median individual, and the
tail is a separate question that was never asked.

Measured effect on woelfel stool, day-0 concentration across 50,000 agents:

| | median | 99th | 99.9th | max |
|---|---|---|---|---|
| before | 6.08 | 12.84 | 16.55 | 27.46 |
| after | 6.36 | 11.13 | 12.68 | 15.43 |

**The trade-off is real and accepted.** The level summary moves from the
geometric mean of `c0` to the arithmetic mean of the subjects' log10 heights.
Against the subjects' own pointwise median curve over a 10-fit sample, the new
coordinate is closer on only 3 of 10 (median RMS gap 1.10 vs 1.16 log10; mean
2.36 vs 1.85). It is accepted because the costs are not the same size: the median
individual moves by a fraction of a log10, the simulated maximum by twelve orders
of magnitude, and simulation is what the catalog is for.

This changes the meaning of every exponential fit's `population_mean`, so it
requires `make catalog`. The coordinate names are recorded per fit and checked on
load, so a stale catalog fails loudly rather than being misread.

## Cause 2: two-stage over-dispersion

`Σ_empirical ≈ Σ_between + E[within-subject estimation error]`. Two-stage
estimation does not shrink individual estimates toward the population mean, so
the fitted covariance absorbs each subject's estimation error and every simulated
cohort is too dispersed — the more so the fewer observations per subject. Even
after the coordinate fix, woelfel stool's 99.9th percentile is 12.68 against its
subjects' maximum of 9.00.

The principled correction is variance deconvolution (Global Two-Stage): subtract
the mean within-subject estimation covariance, then project back to PSD. That
needs per-subject uncertainty the fitter does not currently estimate, and it
cannot be applied at simulation time because the catalog deliberately does not
serialize `subject_params`. **Deferred**, and recorded here as the principled fix
rather than left implicit.

What ships instead is an explicit dial.

### `dispersion`

```python
simulate_shedding(fit, n_individuals=1000, times=t, dispersion=0.7)
```

Scales the covariance by `dispersion ** 2`, so the cohort's spread scales by
`dispersion` while its centre and correlation structure are untouched. `1.0`
(default) simulates the fitted population exactly as estimated, so no existing
call changes behaviour. `0.0` gives every agent the median individual.

It is honest about being a judgement, not an estimate: there is no automatic way
to choose it. What makes it defensible rather than a fudge is that the two-stage
bias only ever runs one way — the fitted spread is too wide, never too narrow —
so shrinkage is the only direction worth offering.

Measured on woelfel stool, day-0 concentration:

| dispersion | 99th | max | top 1% share of load |
|---|---|---|---|
| 1.0 | 11.13 | 15.43 | 98% |
| 0.7 | 9.70 | 12.71 | 82% |
| 0.5 | 8.75 | 10.90 | 50% |

At 0.5 the 99th percentile (8.75) sits just under the subjects' observed maximum
of 9.00.

**Truncation was considered and rejected** as the primary dial. Measured before
the coordinate fix, cutting draws at 3 Mahalanobis SD left the top-1% load share
at 99%; only a 1.5 SD cut brought the tail into range, and that discards a third
of the cohort and flattens the genuine kinetic heterogeneity the user wants to
keep. Truncation fights the symptom; with the transform fixed, a mild scale
factor does the work.

`dispersion` is threaded through `SheddingFit.sample_params` and
`SheddingEnsemble.sample_params`. Under `method="mixture"` it shrinks each
component around its own mean, leaving the between-*study* spread of those means
intact — which is the right behaviour, since between-study heterogeneity is not
an estimation artefact.

## The diagnostic band

`plot_fit_diagnostic` gains `show_band=True`, shading the central 5–95% of a
simulated cohort behind the observations, plus `dispersion`, `band_quantiles` and
`n_simulated`. The median individual alone says nothing about whether the
*spread* is right, which is most of what separates a usable fit from an unusable
one — and over-dispersion of this size is invisible until it is drawn.

The band deliberately does **not** expand the y axis. It can reach far past any
observation, and letting it set the range would squash the data into a sliver; it
is clipped instead, so an over-dispersed fit reads as a band filling the panel
rather than as a flattened plot. Its seed is fixed, so a page redrawn from the
same fit is identical.

## Testing

TDD throughout.

- The exponential height coordinate is linear in log10 concentration: a
  coordinate `k` maps to a curve peaking at exactly `k` log10. Under the old
  coordinates a height of 6.0 produced 175 log10, which is the bug in one
  assertion.
- Exponential coordinates round-trip exactly.
- `dispersion=1.0` is byte-identical to omitting it; below 1 narrows the cohort;
  the spread scales by `dispersion`; `0.0` yields the median individual exactly;
  the median is preserved; negatives raise.
- The band is present by default, absent with `show_band=False`, and narrows with
  `dispersion`.

Four existing tests encoded `median_params == exp(population_mean)` for the
exponential model and were rewritten against `from_population_coords`, that
identity now holding only for `a0`.

## Amendment: over-extrapolated subjects (same day)

Reported immediately after the coordinate change shipped: on
`fajnzylber2020sars` most observations sat *below* the median-individual curve.
They did — 87% of them on the nasopharyngeal analyte, with a −2.40 log10 median
residual.

**The coordinate change did not create this, but it exposed it.** The exponential
model's peak is at `t = 0` by definition, and these studies begin sampling days to
weeks later, so a subject with a steep fitted decay is extrapolated backwards
through many half-lives. One nasopharyngeal subject with a 0.30-day half-life,
first sampled on **day 30**, implies `10**33` gc/mL — a hundred half-lives of
extrapolation — against its own highest reading of `10**2.7` and an analyte
maximum of `10**5.5`. While the model was summarized in `log c0`, the logarithm
compressed such subjects enough to hide them; summarizing in `peak_log10` makes
the coordinate linear in log10, so one of them dominates the mean outright.

They were never harmless: under the old coordinates Plasma was already 80% below
its curve with a −2.52 residual. Across the repository, 66 subjects in 34 of 92
fits imply peaks more than 3 log10 above anything their analyte recorded, reaching
`10**47`.

`_over_extrapolated_subjects` therefore excludes a subject from the population
summary when its implied peak exceeds the analyte's highest observed
concentration by more than `_MAX_PEAK_ABOVE_OBSERVED = 3.0` log10 — a
thousandfold. It composes with the existing degenerate flag: the subject stays in
`subject_params`, flagged and warned about, and if too few survive
`require_estimable_population` refuses the fit.

The threshold is referenced to the data rather than fixed, because what counts as
absurd depends on what the assay can see. Results:

| fit | before | after | median residual |
|---|---|---|---|
| Nasopharyngeal | 87% below | **48%** | −2.40 → **+0.08** |
| Oropharyngeal_PBS | 84% below | **47%** | −0.81 → **+0.02** |
| Sputum | 53% below | 53% | −0.35 → −0.17 |
| tsang2016 NPSOPS (440 subj) | 55% | 55% | −0.27 → −0.27 |

The catalog goes from 83 fits to 81. `Plasma_SARSCoV2_N` and `Urine_SARSCoV2_N`
are now refused as `degenerate_fit`, which is the honest outcome: 3 of plasma's 4
subjects were extrapolation artifacts and its median individual claimed
`10**13.8` gc/mL in plasma.

**This fixed centring, not over-dispersion.** `woelfel2020virological` stool has
no gated subjects and its simulated tail is unchanged (99.9th percentile 12.68,
top 1% of agents still carrying 98% of load). The `dispersion` dial remains the
tool for that, and variance deconvolution the principled fix.

## Out of scope

**Variance deconvolution** — the principled fix for cause 2, described above.

**Re-examining the gamma coordinates.** They already avoid this failure mode.
The bimodal-peak problem on `woelfel2020virological` stool gamma is a different
pathology (a unimodal Gaussian summarizing a bimodal population) that no
reparameterization addresses.
