# Shedding-curve modelling: methods

How the fitted estimates in `shedding_hub/data/shedding_catalog.yaml` are
produced, what they mean, and where they should not be trusted.

The catalog currently holds **144 fits over 48 studies and 94 analytes**: 93
exponential, 32 gamma, 19 gamma_shifted, across 10 biomarkers.

## 1. Models

All three are expressed on the log10 scale, which is the scale the likelihood is
evaluated on. `c0` is a natural-log-scale intercept, so a log10 concentration is
`c0 / ln(10)`.

| model | curve | parameters | peak |
|---|---|---|---|
| `exponential` | `c(t) = c0 · e^(−a0·t)` | `a0`, `c0` | `t = 0` |
| `gamma` | `c(t) = c0 · t^b0 · e^(−a0·t)` | `a0`, `b0`, `c0` | `t = b0/a0` |
| `gamma_shifted` | `c(t) = c0 · (t−t0)^b0 · e^(−a0·(t−t0))` | `a0`, `b0`, `c0`, `t0` | `t = t0 + b0/a0` |

`a0` is the decay rate, giving a half-life of `ln(2)/a0` — exact for the
exponential model, asymptotic for the other two once `t` is well past the peak.
`b0` governs the rise. `t0` is the onset of shedding.

**Choosing between them is not an AIC comparison.** The models are fitted to
different observation sets, and AIC only compares models fitted to the same
data. `gamma` discards every reading at `t ≤ 0`, where its curve is undefined;
`gamma_shifted` keeps the detected ones. On `kissler2021viral` that is 2072
observations against 1679, and the resulting AICs are not commensurable. The
selection rule is data availability, encoded in the gates below, not fit
statistic. Compare `n_measurements` before comparing `aic`.

`gamma_shifted` exists because the reference event is not the same event across
studies — the catalog spans symptom onset, enrollment, confirmation date,
vaccination and hospital admission — and because `gamma` was discarding **26,023
detected measurements at exactly `t = 0`**. Its `t0` is the quantity that makes
those five reference events commensurable.

## 2. Estimation

**Censored maximum likelihood.** Roughly 37% of measurements in the repository
are reported `negative`. They enter the likelihood as left-censored
observations — contributing the information that the concentration was below the
limit — rather than being dropped, which biases decay rates slow and inflates
simulated late-phase shedding.

The censoring limit is the analyte's declared limit of quantification, then its
limit of detection, used as declared. Reported positives *below* that limit are
kept as observed data: a number the assay reported is still a measurement, and
the limit describes only the value a `negative` is known to lie under.

**Two-stage, fitted jointly.** One optimisation per analyte over every subject's
parameters plus a single shared measurement-error SD `sigma`, by L-BFGS-B on
`theta`. For the positive-scale parameters `theta = log(parameter)`, which makes
positivity automatic; `gamma_shifted`'s `t0` is a time on the whole real line and
is carried untransformed, with its constraint imposed as a bound (§4).

**Population summary.** Subjects are summarised as a multivariate normal, but
not in their log-parameters — see §3. The summary uses only non-degenerate
subjects (§4).

## 3. Population coordinates

A population summary averages subjects and treats the result as a Gaussian. That
is only defensible in coordinates where the subjects form a compact,
roughly-elliptical cloud, and the natural log-parameters are not those
coordinates for any of the three models.

| model | coordinates |
|---|---|
| `exponential` | `log_a0`, `peak_log10` |
| `gamma` | `log_a0`, `log_peak_day`, `peak_log10` |
| `gamma_shifted` | `log_a0`, `log_rise_days`, `peak_log10`, `t0` |

A cycle-threshold fit names its height `peak_cycles` instead — it is cycles
below `CT_REFERENCE`, not a log10 concentration. The temporal coordinates are
unchanged, because they mean the same thing on either scale.

**For the rise-and-fall models, `c0` is not separately meaningful.** It is the
concentration at `t = 1`, so a plausible value depends entirely on `b0`. On
`woelfel2020virological` stool, `b0` spans 0.024 to 6.80 while `c0` counter-varies
from 19.6 down to 2.5, `corr(log b0, log c0) = −0.63`. The subjects lie on a
curved ridge, and a coordinate-wise average lands off it: the naive summary
peaked at 3.33 log10, *below every one of the six subjects it summarised*.

**For the exponential model the reason is the tail, not the ridge.** Since
`log10 c(0) = c0/ln(10)`, modelling `log c0` as normal makes the log10
concentration lognormal and the concentration a *double* exponential of the
draw. Simulated day-0 concentrations reached a 99.9th percentile of 10^16.6 and a
worst draw of 10^25, against eight real subjects topping out at 10^9, with the top
0.1% of a cohort carrying essentially all of its shed load. Taking `peak_log10`
as a coordinate makes a draw `k` units above the mean land exactly `k` log10
above it, and the worst draw falls to 10^15.4.

This costs a little median-individual accuracy — measured over a 10-fit sample
the new coordinate was closer on only 3 of 10 — and is accepted because the two
costs differ by orders of magnitude and simulation is what the catalog is for.

`median_params` is the population mean mapped back to parameters. Because each
coordinate is normal, that is the *median* individual, not a compromise. Note it
is not the population's mean trajectory: to aggregate load across a cohort,
simulate rather than scaling it up.

## 4. What is excluded, and why

### Readings

| rule | applies to | reason |
|---|---|---|
| `t < −5` days | all models | Every measurement in the repository earlier than about day −3 is censored, and across all 71 datasets the earliest detected reading anywhere is day −5, so no measured value is ever discarded by this rule. Under a decay-only model those censored points are near-impossible and distort the fit: `tsang2016individual` NPSOPS had its pre-event readings a median 4.01 log10 below its own curve. |
| `t ≤ 0` | `gamma` | `ln(t)` is undefined. Not a judgement — there is nothing to evaluate. |
| `t ≤ 0` and censored | `gamma_shifted` | The curve dives toward −∞ near `t0`, so "below the limit" there is explained for free and `t0` becomes a support parameter pulled onto its own bound. A *detected* reading at the same time is kept, and repels `t0` instead. |

### Subjects

A subject is excluded from the population summary — but kept in
`subject_params`, flagged — when its fit is an artifact rather than an estimate:

- **collapsed**: any parameter at or below 0.01
- **pinned**: any parameter at the optimizer's bound of `e^25`
- **runaway decay**: implied half-life below 0.1 days
- **over-extrapolated**: implied peak more than **3 log10** above the highest
  concentration the analyte ever recorded

The last two thresholds are calibrated in log10 units and applied to whatever
scale the fit is on, so on a cycle-threshold fit they are stricter by roughly
the standard-curve slope. See *Cycle-threshold fits* below.

The last of these is judged against the data rather than a fixed bound, because
what counts as absurd depends on what the assay can see. It matters: when the gate was
introduced, 66 subjects across 34 fits implied peaks beyond that ceiling,
reaching 10^47. One
`fajnzylber2020sars` nasopharyngeal subject with a 0.30-day half-life, first
sampled on **day 30**, implied 10^33 gc/mL — a hundred half-lives of backward
extrapolation — while its own highest reading was 10^2.7. Averaging those in left
87% of that analyte's observations *below* its own median-individual curve;
excluding them brought it to 48%.

For `gamma_shifted`, `t0` is additionally bounded half a day below each subject's
earliest reading, so a curve is never undefined at its own observations.

### Analytes

`catalog.skipped` records a reason for every refusal. Current counts:

| reason | n | meaning |
|---|---|---|
| `ct_units` | 207 | Cycle thresholds. Fittable directly with `fit_shedding_model`, but excluded from the catalog: their heights are cycles below `CT_REFERENCE`, not log10 concentrations, so an ensemble cannot average them together |
| `no_rise_observed` | 63 | Fewer than 50% of subjects peaked later than their first reading, leaving `b0` unidentifiable |
| `no_pre_event_readings` | 61 | `gamma_shifted` only: no detected reading at or before day 0, so `t0` has nothing to locate |
| `too_few_subjects` | 50 | No subject cleared the per-subject minimum |
| `degenerate_fit` | 31 | Too few subjects survived the checks above |
| `too_few_subjects_for_population` | 24 | Too few subjects to estimate a covariance |
| `non_pathogen_biomarker` | 12 | crAssphage, PMMoV, mtDNA — indicators, not shed pathogens |
| `no_positive_measurements` | 4 | Nothing ever detected |
| `no_data_after_reference_event` | 1 | Sampling stopped at day 0, so no post-event trajectory exists |

### Cycle-threshold fits

`Ct = α − β·log10 C` is affine, so a gamma curve in concentration is a gamma
curve in Ct. Cycle-threshold analytes are fitted on `CT_REFERENCE − Ct` —
cycles below a fixed reference of 40, which rises with viral load exactly as a
log10 concentration does.

Because `a0` and `b0` both come back multiplied by the unknown standard-curve
slope, that slope cancels in `b0/a0`. **Peak time, onset and rise duration are
therefore invariant as algebra, with no assumption about PCR efficiency** —
where the two fits converge to corresponding solutions they agree to about
0.01 days. Decay rate, half-life and peak height are not invariant:
`SheddingFit.comparable_with` says which is which for any pair of fits.

That invariance is a statement about the parameters, not a promise about the
estimates. Measured on the studies that report both scales, the two fitted
population peak times are **not interchangeable**; see *How far invariance
actually gets you* below before comparing one against the other.

Onset is invariant for a different reason from the other two, and it is worth
being exact about which. Rise duration is `b0/a0` and peak time is that ratio,
so for them the slope cancels. `t0` is a *time*: the transform rescales and
offsets the response axis and leaves the time axis untouched, so an onset in
days is the same number on either scale — there is no slope in it to cancel.

A Ct fit's height coordinate is named `peak_cycles`, not `peak_log10`, and
`catalog.table` carries a `value_type` column, so a height in cycles cannot be
read as a log10 concentration by accident. `simulate_shedding` **refuses a Ct
source** for the same reason: turning cycles into a concentration needs that
assay's standard curve, which the studies do not report.

### What is *not* invariant across scales: which subjects are kept

Two of the subject-level gates in §4 are calibrated in log10 units and applied
unchanged to the Ct response, which is in cycles:

| gate | calibrated as | on the Ct scale |
|---|---|---|
| over-extrapolation margin | 3 log10 above the highest observation | 3 *cycles*, ≈ 0.9 log10 at a slope of 3.3 |
| minimum half-life | 0.1 days | ≈ 0.1 days on a decay rate carrying the slope, so ≈ 0.33 days in concentration terms |

Both are therefore roughly **β times stricter** on the Ct scale, β being the
standard-curve slope (about 3.3 for a fully efficient assay). They are
deliberately not rescaled: dividing by a nominal slope would put an assumed PCR
efficiency back into the fitting path, and this package fits Ct directly
precisely so that no efficiency is ever assumed.

The consequence is that **subject retention is not scale-invariant even though
the parameter mathematics is.** It is not a marginal effect: on
`cdc2024nhphrn`, 40% of the subjects retained on the Ct scale sit within twice
the minimum half-life, against 2.8% on the concentration scale, and 20 subjects
are excluded as degenerate on Ct alone against 4 on concentration alone.

### How far invariance actually gets you

Three studies report both a Ct and a concentration analyte on the same specimen
and subjects. `teunis2015shedding` supports only the exponential model, whose
`peak_day` is pinned to 0 by construction, so it tests nothing; the evidence is
`kissler2021viral` and `cdc2024nhphrn`. Comparing subjects retained as
non-degenerate on *both* scales, so that retention differences are excluded:

| study | model | n | median \|Δ peak day\| | IQR | max |
|---|---|---|---|---|---|
| `kissler2021viral` | gamma | 33 | 0.46 d | 0.07–1.66 | 105.9 |
| `kissler2021viral` | gamma_shifted | 35 | 0.57 d | 0.24–0.86 | 4.9 |
| `cdc2024nhphrn` | gamma | 51 | 0.19 d | 0.05–0.60 | 4.0 |
| `cdc2024nhphrn` | gamma_shifted | 48 | 1.03 d | 0.45–1.50 | 2.7 |

The algebra is not what fails. `a0_ct/a0_conc` should equal the standard-curve
slope, and `kissler2021viral` reports one of −3.61: among the nine subjects
whose ratio lands within 2% of it, the median peak-time difference is
**0.013 days**. Between 2% and 10% it is 0.38 days, and beyond 10% it is
1.97 days — a rank correlation of 0.83 between ratio error and peak
disagreement. What diverges is the *optimiser*, on subjects too weakly
identified to pin down: `kissler2021viral` is 91% censored.

Two cautions on top of that. Under `gamma` the disagreement is directional —
Ct peaks earlier than concentration in both studies — which convergence noise
alone does not explain, and the over-extrapolation gate above is a plausible
but unconfirmed cause. And in both studies the concentration column is derived
from Ct through a standard curve rather than measured independently, so these
numbers test the transform's algebra, not two independent assays.

**Practical guidance.** Treat a Ct fit's peak time as comparable with a
concentration fit's in kind, not in value. It is sound for ranking and for
order-of-magnitude reasoning, and for well-identified subjects it is accurate;
it is not a drop-in substitute at the study level, and a difference of half a
day between scales is ordinary rather than a sign of error.

## 5. Simulation

`simulate_shedding` draws each individual's coordinates from
`MVN(population_mean, population_cov)`, maps them back to parameters, and
evaluates that individual's own curve. `incubation_period` shifts the clock from
the reference event to infection; `include_measurement_error` adds `N(0, sigma)`
assay noise, off by default because an agent-based model wants the concentration
the host is shedding, not what an assay would report.

**Simulated cohorts are over-dispersed, always in the same direction.** Two-stage
estimation does not shrink individual estimates toward the population mean, so
`population_cov` carries within-subject estimation error on top of true
between-subject variance. The bias only runs one way, which is what makes the
`dispersion` argument a defensible correction rather than a fudge: it scales the
covariance by `dispersion²`, preserving the centre and the correlation structure.
There is no automatic way to choose it — it is a judgement about how much of the
fitted spread is real.

On `woelfel2020virological` stool, day-0 concentration:

| dispersion | 99th pct | max | top 1% share of load |
|---|---|---|---|
| 1.0 | 11.13 | 15.43 | 98% |
| 0.7 | 9.70 | 12.71 | 82% |
| 0.5 | 8.75 | 10.90 | 50% |

The principled fix is variance deconvolution (Global Two-Stage): subtract the
mean within-subject estimation covariance. It is **not implemented** — it needs
per-subject uncertainty the fitter does not estimate and the catalog does not
serialize.

### Reference events are not interchangeable

The catalog spans seven reference events, in three classes. `inoculation` and
`vaccination` *are* the exposure, so nothing separates them from time zero.
`symptom onset` is a natural-history landmark, a defined and documented offset
from infection. `enrollment`, `confirmation date`, `hospital admission` and
`treatment` are administrative: they record when a subject entered a study, was
tested, or was admitted, which depends on testing behaviour and health-system
access rather than on their infection.

Only the landmark class earns an infection time origin. `simulate_shedding`
warns when an incubation period is applied to either of the others, and records
`time_origin` as `"<event>_shifted"` rather than `"infection"`. `shedding_for`
prefers the classes that can be anchored, which is why it will pass over a
better-supported fit measured from an administrative date.

These are two different ideas and the ranking uses the weaker one: exposure and
landmark events can both be *placed* on an infection timeline, which is what
`shedding_for` prefers, but only a landmark *earns*
`attrs["time_origin"] == "infection"` — an exposure event already is the
exposure, so there is no incubation period to shift back through.

## 6. Reading an estimate honestly

- **`peak_log10` is evaluated at the peak, which for the exponential model is
  `t = 0` by definition.** Read it beside `median_first_observed_day`: when that
  is well above zero, the peak is a backward extrapolation to a time most
  subjects were never observed at, not a measured concentration.
- **A high `pct_censored` is not a defect.** It means the analyte is rarely
  detected, and a low `peak_log10` there is the fitter reporting honestly.
- **`converged=False` disqualifies the rest of the row.** No fit in the shipped
  catalog currently carries it. The five that used to — all
  `natarajan2022gastrointestinal` — were not failures of the model but of the
  evaluation budget: L-BFGS-B stopped for want of allowance while still
  descending. The optimizer now continues a round that ended only for that
  reason, and all five converge. One had been badly short rather than
  marginally so, gaining 76 log-likelihood units when allowed to finish.
- **`n_subjects` counts subjects fitted, not subjects summarised.**
  `n_degenerate_subjects` is the difference.

## 7. Reference tables

`make parameters` writes both, from the shipped catalog, fitting nothing:

- `docs/shedding_parameters.json` — one record per fit, each carrying the
  population mean and covariance, `sigma` and the censoring limit. Enough to
  simulate from without this package and without refitting.
- `docs/shedding_parameters.csv` — the flat browsing view, one row per fit.
  Drops the covariance, which is a k×k matrix per fit and does not belong in a
  spreadsheet column.

In Python, `load_shedding_catalog().table` is the same flat view as a DataFrame,
and `catalog.select(...)` returns a fit ready to simulate.

## 8. Design records

Each decision above was recorded when it was made:

- `docs/superpowers/specs/2026-07-18-shedding-simulation-design.md` — the fitting
  and simulation pipeline
- `docs/superpowers/specs/2026-07-27-catalog-fit-plots-design.md` — `plot_catalog_fits`
- `docs/superpowers/specs/2026-07-29-catalog-fit-diagnostics-design.md` — `plot_fit_diagnostic`
- `docs/superpowers/specs/2026-07-29-simulation-dispersion-design.md` — the
  exponential coordinate change, `dispersion`, and the over-extrapolation gate
