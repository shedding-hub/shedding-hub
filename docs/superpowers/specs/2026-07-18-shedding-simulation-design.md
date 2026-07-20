# Shedding Simulation from Fitted Models — Design

**Date:** 2026-07-18
**Status:** Draft (awaiting review)

## Problem

A user working on wastewater surveillance simulation asked for a way to simulate
biomarker shedding as a function of time since infection, driven by real Shedding
Hub data. The concrete need: an agent-based model that infects an agent and then
needs that agent's shedding trajectory — sampled from a realistic distribution of
individual trajectories, not a single average curve.

Nothing in the package supports this today. The existing modules
(`shedding_duration`, `shedding_peak`, `stats`, `viz`) all *describe* observed
data; none produce new synthetic individuals.

The statistical groundwork already exists but lives in R. The site repo
(`shedding-hub.github.io/tutorials`) has a hierarchical Bayesian Rstan workflow
fitting an exponential-decay and a gamma model with a censored likelihood, and a
JAGS implementation of the Teunis within-host model. This work ports the first
two models to Python, fits them across the repository, and adds the simulation
step.

The intended user journey is: **browse a table of available fitted estimates →
pick one study, or an ensemble across studies → simulate individuals.**

## Scope

**In:** exponential and gamma models; censored maximum-likelihood fitting per
analyte; a browsable catalog of estimates shipped with the package; simulation
from a single study, from a user-chosen subset of studies, or from every matching
study via a cross-study ensemble; one plotting helper.

**Out (deliberately):** the Teunis two-compartment model with an estimated
shedding-onset offset; MCMC/posterior inference; Ct-scale modelling. Each is a
plausible follow-up and the interfaces are designed not to block them.

## Models

Both are taken directly from the Rstan tutorial and expressed on the log10 scale,
which is the scale the likelihood is evaluated on.

**Exponential decay** — pure decay, appropriate when sampling starts at or after
peak shedding (the common case in shedding studies):

```
c(t) = c0 * exp(-a0 * t)
log10 c(t) = (c0 - a0 * t) / ln(10)
```

**Gamma** — rise then fall, peaking at `t = b0 / a0`:

```
c(t) = c0 * t^b0 * exp(-a0 * t)
log10 c(t) = (c0 + b0 * ln(t) - a0 * t) / ln(10)
```

Following the tutorial's parameterization, `c0` is on the natural-log scale (so
log10 concentration at `t = 0` is `c0 / ln(10)`), and all of `a0`, `b0`, `c0` are
strictly positive. Per-subject parameters are therefore modelled on the log
scale:

```
theta_i = (log a_i, log c_i)            exponential
theta_i = (log a_i, log b_i, log c_i)   gamma
theta_i ~ MVN(mu, Sigma)
```

`theta_i ~ MVN(mu, Sigma)` is the population model. It is what makes simulation
possible: drawing a new `theta` from it yields a new plausible individual.

`t = 0` is the dataset's own reference event (usually symptom onset), matching the
tutorial's stated assumption that shedding begins at symptom onset.

The exponential model is fitted for every analyte. **The gamma model is fitted
only where a rise is actually observed** — see the rise gate below — so the
catalog carries both models where both are identifiable, and users compare them
on AIC rather than being forced into one.

## The fitting unit is the analyte

A fit is estimated per **(dataset, analyte, model)**. Biomarker, specimen, and
reference event are *selection keys* carried alongside, not the grouping level.

This matters because 30 strata in the repository contain more than one analyte
for the same (biomarker, specimen, reference event), for three different reasons:

- **Ct paired with concentration** (`kim2020viral`, `gutierrez2021nosocomial`,
  `cdc2024nhphrn`) — resolved automatically by rejecting Ct analytes.
- **Different gene targets on the same subjects** (`arts2023longitudinal` measures
  both N and ORF1a in stool, same unit) — pooling these would enter each subject
  twice with correlated repeat measures.
- **Genuinely different exposures** (`cowley2017rotavirus` has three vaccine
  doses; `jacobsen2022differentiation` has RV1 and RV5) — these should never be
  merged; they are distinguished by the schema's `dose` and `vaccine_type`
  fields.

Fitting per analyte sidesteps all three. Each fit records `analyte`,
`dataset_id`, `biomarker`, `specimen`, `reference_event`, `unit`, and, where
present, `gene_target`, `dose`, and `vaccine_type`.

`unit` is treated as a hard constraint everywhere: estimates in `gc/mL`,
`gc/dry gram`, and `pfu/mL` are not comparable and are never combined.

## Fitting

A single joint maximum-likelihood fit over all per-subject parameters and one
shared measurement-error standard deviation `sigma`, using the same censored
likelihood as the Stan code:

```
maximize  sum_observed  log N(y_ik | m_i(t_ik), sigma)
        + sum_censored  log Phi((L - m_i(t_ik)) / sigma)
```

where `m_i` is the model's log10 prediction for subject `i`, `y_ik` are observed
log10 concentrations, and `L` is the log10 censoring limit. The censored term is
the direct analogue of Stan's `target += normal_lcdf(censlim | gene_cen, sig_obs)`.

Handling censoring this way is the single most important correctness property
here: **37% of all measurements in the repository are `negative`** (23,426 of
62,487). Dropping them — the naive alternative — biases estimated decay rates
toward slower decay and inflates simulated late-phase shedding, which is exactly
the quantity a wastewater model is most sensitive to.

Optimization is `scipy.optimize.minimize` with `L-BFGS-B` over the concatenated
parameter vector, with `sigma` optimized as `log sigma` for positivity. Subject
parameters are initialized from a per-subject ordinary least squares fit to the
uncensored points.

The population distribution is then estimated from the fitted per-subject
parameters:

```
mu    = mean(theta_i)
Sigma = cov(theta_i, ddof=1)
```

### Why two-stage, and what it costs

The Stan model estimates individual and population parameters jointly, so
individual estimates shrink toward the population mean. This two-stage MLE does
not shrink. Consequently `Sigma` absorbs within-subject estimation error and
**overestimates true between-subject variance**, the more so when subjects have
few observations. Simulated cohorts will be somewhat more dispersed than reality.

This is an accepted, documented trade-off, not an oversight: it buys a fast,
deterministic, dependency-light fit. It is recorded in the docstring and in the
Limitations section below, and it is the main reason to add a Bayesian backend
later. `SheddingFit` carries an explicit `method` field so a future
`method="bayes"` can populate the same structure from posterior draws.

## Data preparation and edge cases

Measurements are extracted per participant for the analyte being fitted, then:

- **Ct-unit analytes are rejected.** 52 analytes in the repository use
  `cycle threshold`, which is inversely related to concentration and already on a
  log scale; neither model applies. In catalog building they are skipped with a
  recorded reason; in a direct `fit_shedding_model` call they raise `ValueError`.
- **`negative` values become censored observations** at the analyte's declared
  `limit_of_quantification` (falling back to `limit_of_detection`). The declared
  limit is used **as-is**, and every reported positive is kept as observed data —
  **including positives below the limit**. A value the assay reported below its
  limit of quantification is still a measurement (detected, if less precisely),
  so it is used rather than discarded or re-coded as censored; the limit
  describes only the value a `negative` is known to lie below. The
  observed-value likelihood term never references the limit, so an observed
  positive sitting below it is well defined. 22 analytes carry at least one such
  sub-limit positive (e.g. woelfel stool has two below its LOQ of 100); all are
  kept.
- **Censoring-limit fallback** applies only when **neither** limit is declared
  (both are the literal string `unknown`): then it falls back to just below the
  smallest observed positive and warns, so any `negative` still sits below the
  resolved limit. A declared limit is otherwise always honored, even if it sits
  above some or all observed positives.

  The schema also allows a **per-measurement** `limit_of_quantification`. It is
  deliberately **not** implemented: every one of the 271 measurements declaring
  one belongs to `arts2023longitudinal`'s `stool_crAssphage` analyte, which is
  now excluded as a non-pathogen indicator (below). No fitted analyte declares a
  per-measurement limit, so a single scalar limit per fit is correct and the
  likelihood stays simpler. Revisit if a future dataset declares them on a
  pathogen analyte.
- **The gamma rise gate.** The gamma model's `b0` controls the rise, so it is
  unidentifiable when a study starts sampling after peak shedding — which is the
  common case in shedding studies. Fitting it anyway produces a meaningless
  `peak_day` and `peak_log10`. Confirmed by profile likelihood during
  implementation: where sampling is entirely post-peak the likelihood is
  monotone in `b0` until another parameter binds, so no amount of better
  initialization helps.

  A subject "observes a rise" if it has at least 3 usable positive observations
  at `t > 0` and its maximum observed value occurs later than its first
  observation. If fewer than **50%** of an analyte's subjects observe a rise, the
  gamma fit is refused with reason `no_rise_observed`. The exponential model is
  unaffected — post-peak sampling is precisely where it applies.

  The fraction is published as a `pct_subjects_with_rise` column so the judgment
  is auditable rather than buried in a threshold. In practice the gate refuses
  38 analyte/model combinations. After the later degeneracy and population-size
  gates the shipped catalog retains 23 gamma fits.
- **Non-pathogen indicator biomarkers are rejected.** `crAssphage`, `PMMoV`, and
  `mtDNA` are fecal-strength and normalization markers, not pathogens shed by
  infected people — they have no time-since-infection trajectory, so a shedding
  curve fitted to one is meaningless. This excludes exactly 4 analytes across 2
  studies (`arts2023longitudinal` stool_PMMoV and stool_crAssphage,
  `liu2024longitudinal` stool_PMMoV and stool_mtDNA); neither study loses its
  SARS-CoV-2 analytes. Vaccine-strain biomarkers such as `rotavirus vaccine`
  stay in scope: live-attenuated shedding after vaccination is a real trajectory
  with `vaccination` as its reference event.
- **Qualitative positives are dropped with a warning.** `positive`,
  `weak positive`, `strong positive` and `inconclusive` carry no numeric value and
  cannot enter a normal likelihood. This is 172 measurements repository-wide
  (~0.3%), so the loss is negligible and silence would be worse.
- **`time: unknown` rows are dropped.**
- **Non-positive times under the gamma model are dropped with a warning.**
  `ln(t)` is undefined at `t <= 0`. Under the model's own assumption — shedding
  starts at the reference event — these observations predate shedding onset and
  carry no information. The exponential model is defined for all `t` and keeps
  them.
- **Subjects with too few usable observations are excluded with a warning**
  (`min_observations`, defaulting to the number of per-subject parameters: 3 for
  gamma, 2 for exponential). Because `sigma` is shared across subjects rather
  than estimated per subject, a subject does not need residual degrees of freedom
  of its own. If no subject survives, raise `ValueError`.
  This also cleanly excludes cross-sectional data: `jones2021estimating` holds
  25,378 participants averaging ~1.2 samples each and cannot support per-subject
  trajectory fitting at all.

Every exclusion is counted and surfaced on the fit object
(`n_excluded_subjects`, `n_dropped_measurements`) rather than being silently
applied.

## The catalog

`SheddingCatalog` is the browse-and-select surface.

**`catalog.table`** — a `DataFrame`, one row per (dataset, analyte, model):

| column | meaning |
| --- | --- |
| `dataset_id`, `analyte` | which fit |
| `biomarker`, `specimen`, `reference_event`, `unit` | selection keys |
| `gene_target`, `dose`, `vaccine_type` | disambiguators, where present |
| `model` | `exponential` or `gamma` |
| `n_subjects`, `n_measurements`, `pct_censored` | how much data backs it |
| `a_median`, `b_median`, `c_median` | model parameters of the median individual |
| `sigma` | measurement error SD (log10) |
| `peak_day`, `peak_log10`, `half_life_days` | interpretable summaries |
| `aic`, `converged` | fit quality |

**Everything in the table describes the median individual, and says so.** Because
`theta = log(a, b, c)` is normal, the parameters themselves are lognormal, so
`exp(mu)` is exactly their median — not their mean. Reporting these as medians is
therefore the accurate label rather than a compromise, and it sidesteps the
Jensen's-inequality trap of implying a mean. The one thing that remains true and
is documented: the median individual's *trajectory* is not the population's mean
trajectory, so a user summing simulated load across a cohort should simulate,
not scale up this row.

The interpretable columns are what make the table selectable — a modeler picking a
row wants "peaks day 4.2 at 6.8 log10 gc/mL, declines with a 1.5 day half-life",
not raw log-parameters. For the gamma model `peak_day = b_median / a_median`; the
exponential model is monotone so `peak_day` is `0` and `peak_log10` is the value
at the reference event. `half_life_days = ln(2) / a_median` describes the
late-phase decline for both.

The table is the browsing surface and carries these human-readable medians; the
full `mu` and `Sigma` needed to actually simulate live on the `SheddingFit`
objects (and in the shipped YAML), so no precision is lost to the summary.

**`catalog.skipped`** — a `DataFrame` of analytes that could not be fitted, with a
reason. Without this a user cannot tell whether a study is missing because it is
unsuitable or because of a bug. The full set, with counts in the shipped catalog:

| reason | count | meaning |
| --- | --- | --- |
| `ct_units` | 104 | cycle-threshold analyte; neither model applies |
| `no_rise_observed` | 38 | gamma refused — the rise is unidentifiable (see gate above) |
| `too_few_subjects` | 24 | no subject had enough observations |
| `degenerate_fit` | 16 | too few subjects survived collapse/runaway detection |
| `too_few_subjects_for_population` | 16 | fewer subjects than parameters, so `Sigma` is not estimable |
| `non_pathogen_biomarker` | 8 | fecal-strength or normalization marker |
| `unknown_analyte` | 0 | the requested analyte is not in the dataset |
| `unexpected_error` | 0 | anything else, with the original message preserved |

The shipped catalog holds **84 fits** (61 exponential, 23 gamma) across 55
datasets, with 206 analyte/model combinations refused.

**`catalog.select(**keys)`** — returns exactly one `SheddingFit`. If the keys match
zero rows, or more than one, raise `ValueError` listing the candidates and the
columns that would disambiguate them. Never silently pick.

**`catalog.ensemble(**keys, dataset_ids=None, weights=..., method=...)`** —
returns a `SheddingEnsemble` over the matching fits, optionally restricted to a
named subset of studies (see below).

### Precomputed and shipped

The catalog is built offline and shipped as package data at
`shedding_hub/data/shedding_catalog.yaml`, loaded by `load_shedding_catalog()`.
YAML keeps it human-diffable and reuses the `pyaml` dependency; at 84 fits of
small float vectors the file is about 166 KB, which is negligible in a wheel.

- `scripts/build_shedding_catalog.py` regenerates it, wired to a `make catalog`
  target alongside the existing `extraction` targets.
- A CI check asserts every dataset in `data/` appears in **either** `table` or
  `skipped`, so adding a dataset without regenerating fails loudly. Datasets that
  are entirely unfittable (Ct-only, or cross-sectional like
  `jones2021estimating`) legitimately have no `table` rows, so checking `table`
  alone would wrongly fail; they must still be accounted for in `skipped`. This is
  a cheap coverage check, not a full refit — refitting everything on every CI run
  would be needlessly slow, and coverage catches the drift that actually happens.
- `pyproject.toml` gains a `package-data` entry so the YAML ships in the wheel.

Users can always build their own catalog from datasets we do not host via
`fit_shedding_models(...)`, including private data.

## The ensemble

### Choosing what goes into it

The user controls which studies contribute, at three levels of specificity:

```python
# 1. one study only
fit = cat.select(dataset_id="woelfel2020virological", analyte="stool",
                 model="gamma")

# 2. a chosen subset of studies
ens = cat.ensemble(
    biomarker="SARS-CoV-2", specimen="stool", reference_event="symptom onset",
    unit="gc/mL", model="gamma",
    dataset_ids=["woelfel2020virological", "wang2020fecal"],
)

# 3. every matching study
ens = cat.ensemble(
    biomarker="SARS-CoV-2", specimen="stool", reference_event="symptom onset",
    unit="gc/mL", model="gamma",
)

# or assemble explicitly, mixing catalog fits with your own fresh ones
ens = sh.make_ensemble([fit_a, fit_b, my_own_fit], weights="equal")
```

`dataset_ids` narrows an otherwise-matching filter; `make_ensemble` takes fits
directly and is the escape hatch for combining catalog estimates with fits from
private data. A **single-component ensemble is legal** and behaves identically to
the underlying fit, so a user can write one code path and vary only how many
studies feed it.

`make_ensemble` still enforces the compatibility rules below — matching unit,
reference event, biomarker, specimen, and model — because those are correctness
constraints, not conveniences of the catalog.

### Combining them

Two methods, sharing the same component fits:

**`method="mixture"` (default)** — each simulated individual first draws a study,
then draws `theta` from that study's MVN:

```
s     ~ Categorical(w_1..w_S)
theta ~ MVN(mu_s, Sigma_s)
```

This preserves between-study heterogeneity, stays multimodal when studies
genuinely disagree, assumes nothing beyond the per-study fits, and lets every
simulated agent be traced to a source study (reported as a `source_dataset_id`
column in the simulation output). `weights` defaults to `n_subjects`, with
`"equal"` available to stop one large study dominating.

**`method="moment"`** — collapses to a single Gaussian by moment matching:

```
mu_ens    = sum_s w_s * mu_s
Sigma_ens = sum_s w_s * Sigma_s + cov_s(mu_s)     # within + between
```

Unimodal by construction, but it still accounts for between-study variance and
gives one tidy MVN to report or hand to another tool.

**One fit per study.** After filtering, if a single study contributes more than
one analyte (`arts2023longitudinal` would contribute both N and ORF1a), raise
`ValueError` listing the candidates and suggesting a narrowing key such as
`gene_target`. Silently averaging gene targets, or silently picking one, would put
an arbitrary scientific choice inside the package.

**Units must match.** A filter spanning mixed units raises, naming them.

**Seeing what went in.** `ens.components` is a DataFrame with one row per
contributing fit, using the same columns as the catalog table, so the user can
inspect exactly which studies and medians the ensemble rests on. A single
median-individual row is well defined only under `method="moment"` (where it is
`exp(mu_ens)`) and is exposed there as `ens.median_params`. For a mixture there is
no closed-form median of the mixed distribution, so rather than report a
misleading one, the docs point users to simulate and take empirical quantiles of
the result — which is the operation they actually want anyway.

## Simulation

`simulate_shedding` accepts either a `SheddingFit` or a `SheddingEnsemble` — both
expose the same `sample_params(rng, n)` interface — draws `n_individuals`
parameter vectors, evaluates the model at the requested times, and returns a tidy
DataFrame:

```
individual_id  time  log10_value  value  detected  [source_dataset_id]
```

- **Measurement error is off by default.** An agent-based model wants the true
  shed concentration; assay noise is a property of sampling, not of the host.
  `include_measurement_error=True` adds `N(0, sigma)` on the log10 scale for users
  who want to emulate observed data.
- **`detected`** flags whether the value is at or above the censoring limit, so
  users can reproduce detection-rate style outputs.
- **Values below the limit are reported as-is** (not clipped) with
  `detected=False`, so downstream mass-balance calculations stay correct. Under
  the gamma model `t <= 0` yields `-inf`; those rows are returned as `NaN` with
  `detected=False`.
- **`seed`** drives a `numpy.random.Generator` so runs are reproducible.

### Time since infection

Fitting happens in the dataset's native reference-event time. `simulate_shedding`
accepts `incubation_period`, which shifts the origin:

- `None` (default) — `time` is days since the fit's reference event.
- a scalar — every individual shares one incubation period.
- an array of length `n_individuals`, or a callable `(rng, n) -> array` — per
  individual, which adds realistic timing variability across the cohort.

When supplied, `time` in the output means **days since infection**, and the curve
is evaluated at `time - incubation_i`. The offset is never assumed: with no
`incubation_period` the output is explicitly in reference-event time, and
`reference_event` is carried on the fit and echoed in the result's `attrs` so it
cannot be misread.

Note that reference events differ in how far they sit from infection —
`symptom onset` (90 analytes), `confirmation date` (24), `enrollment` (18),
`vaccination` (9), `hospital admission` (4) — which is exactly why
`reference_event` is a selection key and an ensemble is only formed within one.

## API

```python
import numpy as np
import shedding_hub as sh

cat = sh.load_shedding_catalog()
cat.table.query("biomarker == 'SARS-CoV-2' and specimen == 'stool'")

# one study ...
source = cat.select(
    dataset_id="woelfel2020virological", analyte="stool", model="gamma",
)

# ... a chosen subset of studies ...
source = cat.ensemble(
    biomarker="SARS-CoV-2", specimen="stool",
    reference_event="symptom onset", unit="gc/mL", model="gamma",
    dataset_ids=["woelfel2020virological", "wang2020fecal"],
)

# ... or every matching study
source = cat.ensemble(
    biomarker="SARS-CoV-2", specimen="stool",
    reference_event="symptom onset", unit="gc/mL", model="gamma",
)
source.components      # which studies contribute, and their medians

# identical from here regardless of which was chosen
traj = sh.simulate_shedding(
    source, n_individuals=1000, times=np.arange(0, 30),
    incubation_period=5.0, seed=42,
)

fig = sh.plot_simulated_shedding(traj, source=source)
```

Fitting fresh, including on private data:

```python
data = sh.load_dataset("woelfel2020virological", local="./data")
fit = sh.fit_shedding_model(data, analyte="stool", model="gamma")
cat = sh.fit_shedding_models([data], models=("exponential", "gamma"))
```

`SheddingFit` is a dataclass carrying `model`, `method`, `population_mean`,
`population_cov`, `sigma`, `subject_params`, `censoring_limit`, the stratum keys,
`converged`, `aic`, and the exclusion counts. `SheddingFit`, `SheddingEnsemble`,
and `SheddingCatalog` all have `to_dict()`/`from_dict()` so a fit can be saved to
YAML and reloaded — the property that lets an ABM fit once and simulate across
many runs without refitting.

## Errors and validation

`ValueError` for: missing dataset keys; Ct units in a direct fit call; unknown
`model`; no subject meeting `min_observations`; `select()` matching zero or many
rows; `ensemble()` or `make_ensemble()` spanning mixed units, reference events,
biomarkers, specimens, or models, or drawing two analytes from one study;
`dataset_ids` naming a study with no matching fit (rather than silently dropping
it, which would quietly shrink the ensemble); `n_individuals < 1`; an
`incubation_period` array whose length differs from `n_individuals`; a
non-positive-definite `Sigma` (possible with very few surviving subjects — the
message says so and suggests an ensemble instead).

`UserWarning` for: dropped qualitative/unknown-time/non-positive-time
measurements; excluded subjects; censoring-limit fallback; optimizer
non-convergence (fit still returned, `converged=False`).

## Testing

New file `tests/test_simulate.py`, following the existing test modules.

- **Parameter recovery.** Simulate from known `mu`/`Sigma`/`sigma`, fit, assert
  recovery within tolerance. Both models. This is the core correctness test.
- **Censoring correctness.** On synthetic data with a known decay rate and an
  artificially high limit, assert the censored fit recovers the truth
  substantially better than a fit that drops censored points. This pins the
  property the whole design rests on.
- **Tutorial agreement.** Fit the exponential model to woelfel subject 3 stool
  (14 positives, 3 negatives at t=20,22,23, LOQ=100 → log10 = 2 — the tutorial's
  exact setup) and assert the estimates land near the tutorial's reported
  posterior (`a0 ≈ 0.74`, `c0 ≈ 20.37`, `sigma ≈ 0.92`). Priors there are flat, so
  the MLE should be close to the posterior mean. This validates the Python port
  against the published R workflow.
- **Ensemble behaviour.** Mixture draws respect weights and populate
  `source_dataset_id`; moment-matched covariance equals within-plus-between on a
  hand-computable two-study example; mixed units raise; two analytes from one
  study raise.
- **Ensemble membership.** `dataset_ids` restricts components to exactly the
  named studies; an unmatched name raises; `make_ensemble` accepts fits not from
  the catalog; a single-component ensemble produces the same distribution as the
  underlying fit (same seed, same draws), which is what lets users keep one code
  path across all three selection levels; `components` lists the contributing
  fits.
- **Median reporting.** Table `a_median`/`b_median`/`c_median` equal `exp(mu)`,
  and `peak_day` for a gamma row equals `b_median / a_median`.
- **Catalog.** `select()` raises on ambiguous and on empty matches, listing
  candidates; the shipped catalog loads, round-trips through
  `to_dict`/`from_dict`, and covers every dataset in `data/` (the CI staleness
  check); `skipped` records a reason for each unfitted analyte.
- **Reproducibility.** Same seed produces identical output; different seeds do
  not.
- **Incubation shift.** Scalar, per-individual array, and callable forms each
  shift the curve as expected; `None` leaves output in reference-event time.
- **Edge cases.** Ct analyte rejected; subject with all-censored measurements;
  subjects below `min_observations`; gamma with `t <= 0`; qualitative values.
- **Plot** returns a `Figure` and closes it, matching the convention in
  `shedding_peak.py`.

## Build order

The pieces layer cleanly, and each stage is independently testable:

1. **Models and censored fitting** for a single analyte — `fit_shedding_model`,
   plus the recovery, censoring, and tutorial-agreement tests. This is the
   scientific core; nothing else is worth building if it does not recover known
   parameters.
2. **Simulation** from a single fit, including the incubation shift.
3. **Catalog** — `fit_shedding_models`, the table with derived columns,
   `select()`, `skipped`, serialization, the build script and `make catalog`.
4. **Ensemble** — mixture and moment methods on top of catalog fits, plus the
   three selection levels (one study, a chosen subset, everything matching) and
   `make_ensemble` for explicitly assembled fits.
5. **Plotting, README, and the CI coverage check.**

## Limitations

Recorded in module docstrings so users encounter them:

1. Two-stage fitting does not shrink individual estimates, so between-subject
   variance is overestimated (see Fitting above).
2. **The gamma model's `b0` is mildly downward-biased at realistic sampling
   densities** — about **0.15 log units** at ~14 observations per subject
   (measured over six seeds, range 0.02–0.27). This is ordinary finite-sample
   maximum-likelihood bias: the estimator is consistent and the bias vanishes as
   sampling density rises. Since `peak_day = b0 / a0`, the catalog's `peak_day`
   runs slightly early for sparsely-sampled studies. A test pins the direction
   and mechanism so it cannot change silently.

   *This figure was originally measured at 0.55 log units. Most of that was not
   bias at all but the parameter-collapse artifact described below, which was
   inflating the apparent effect nearly fourfold. The residual 0.15 is the real
   finite-sample bias.*
3. **Optimizing on the log scale makes zero an absorbing state, and the fitter
   guards against it rather than being immune.** Because parameters are
   optimized as `theta = log(param)`, the chain rule gives
   `dL/d(theta) = param * dL/d(param)`, so the gradient vanishes as a parameter
   approaches zero: a parameter that drifts small can never climb back.
   Discovered when the first repository-wide build produced 37 fits with
   half-lives over a year — some at exactly `ln(2)/1e-6` — because the
   initialization clipped non-positive least-squares seeds straight onto the
   floor. Initialization now never starts at the floor, and subjects whose
   parameters still collapse are detected by magnitude, excluded from `mu`/`Sigma`
   (but kept in `subject_params` with a `degenerate` flag), and counted in
   `n_degenerate_subjects`. A fit with fewer than two surviving subjects is
   refused with reason `degenerate_fit`. **Every synthetic test passed while this
   was shipping 278-day half-lives**, so real-data regression tests now guard it.
4. **`peak_log10` is a population median over draws, not the value at
   `exp(mu)`.** Coordinate-wise averaging of log-parameters across a correlated
   ridge yields a parameter vector no real subject has. `peak_day` and
   `half_life_days` are unaffected — each is a ratio or transform of a single
   lognormal, so the value at `exp(mu)` genuinely is the population median — but
   `peak_log10` is a nonlinear function of all three parameters and was landing
   below almost every subject's own peak. It is therefore computed as the median
   over 10,000 draws from `MVN(mu, Sigma)` with a fixed seed.
5. Point estimates only — no parameter uncertainty propagates into simulations.
   Cohort spread reflects between-individual variation, not estimation
   uncertainty.
6. Both models assume shedding begins at the reference event. Datasets with
   substantial pre-onset sampling (`kissler2021densely`, `kissler2021viral`) are
   poorly served; the Teunis onset-offset model is the future answer.
7. The mixture ensemble represents between-study heterogeneity but does not
   *explain* it — differences in assay, matrix, and population are conflated.
8. The exponential model cannot represent a rise and will mis-fit datasets
   sampled from before peak. Compare AIC against the gamma fit before choosing.

## Files

- `shedding_hub/simulate.py` — models, fitting, catalog, ensemble, simulation,
  plotting. If this grows past a comfortable size, split fitting from simulation;
  the catalog/ensemble types stay with the fitting side.
- `shedding_hub/__init__.py` — export `fit_shedding_model`,
  `fit_shedding_models`, `load_shedding_catalog`, `make_ensemble`,
  `simulate_shedding`, `plot_simulated_shedding`, `SheddingFit`,
  `SheddingEnsemble`, `SheddingCatalog`.
- `shedding_hub/data/shedding_catalog.yaml` — shipped precomputed estimates.
- `scripts/build_shedding_catalog.py` — regenerates the catalog.
- `Makefile` — `catalog` target.
- `pyproject.toml` — add `scipy` and `numpy` to `dependencies` (`numpy` is
  currently only transitive via pandas/matplotlib but is used directly), plus a
  `package-data` entry for the catalog.
- `tests/test_simulate.py` — new tests.
- `README.md` — a short section under "Analyzing the Data" with a doctest-safe
  example (README doctests run in CI).
