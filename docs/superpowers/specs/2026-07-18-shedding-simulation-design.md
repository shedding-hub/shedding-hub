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
two models to Python and adds the simulation step.

## Scope

**In:** exponential and gamma models; censored maximum-likelihood fitting;
population-level parameter distribution; simulation of synthetic individuals;
one plotting helper.

**Out (deliberately):** the Teunis two-compartment model with an estimated
shedding-onset offset; MCMC/posterior inference; cross-study pooling with
study-level random effects; Ct-scale modelling. Each is a plausible follow-up and
the parameter interface is designed not to block them.

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

Measurements are extracted per participant per analyte, after filtering by
`biomarker` and `specimen`, then:

- **Ct-unit analytes are rejected.** 52 analytes in the repository use
  `cycle threshold`, which is inversely related to concentration and already on a
  log scale; neither model applies. Raises `ValueError` naming the offending
  analyte and telling the caller to select a concentration analyte.
- **Mixed units are rejected** for the same reason, with the units listed.
- **`negative` values become censored observations** at the analyte's
  `limit_of_quantification`, falling back to `limit_of_detection`.
- **Censoring-limit fallback.** Either limit may be the literal string
  `unknown`, and a declared limit can exceed the smallest observed positive value
  (the tutorial hit this and hand-set `censlim = 1.96` against a declared limit of
  2). When the resolved limit is unknown or not strictly below the smallest
  observed positive, fall back to just below that smallest positive value and
  warn. A per-measurement `limit_of_quantification` overrides the analyte-level
  value when present.
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

Every exclusion is counted and surfaced on the fit object (`n_excluded_subjects`,
`n_dropped_measurements`) rather than being silently applied.

## Simulation

`simulate_shedding` draws `n_individuals` parameter vectors from
`MVN(mu, Sigma)`, exponentiates to the natural scale, evaluates the model at the
requested times, and returns a tidy DataFrame.

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

- `None` (default) — `time` is days since the dataset's reference event.
- a scalar — every individual shares one incubation period.
- an array of length `n_individuals`, or a callable `(rng, n) -> array` — per
  individual, which adds realistic timing variability across the cohort.

When supplied, `time` in the output means **days since infection**, and the curve
is evaluated at `time - incubation_i`. The offset is never assumed: with no
`incubation_period` the output is explicitly in reference-event time, and
`reference_event` is carried on the fit object and echoed in the result
attributes so it cannot be misread.

## API

Flat functions returning plain objects, matching the package's existing
`calc_*`/`plot_*` style.

```python
import numpy as np
import shedding_hub as sh

data = sh.load_dataset("woelfel2020virological", local="./data")

fit = sh.fit_shedding_model(data, model="gamma", specimen="stool")
fit.population_mean      # mu over log-params
fit.population_cov       # Sigma
fit.sigma                # measurement error SD (log10)
fit.subject_params       # DataFrame, one row per subject

traj = sh.simulate_shedding(
    fit, n_individuals=1000, times=np.arange(0, 30),
    incubation_period=5.0, seed=42,
)
#    individual_id  time  log10_value  value  detected

fig = sh.plot_simulated_shedding(traj, fit=fit, observed=data)
```

`fit_shedding_model` also accepts a list of datasets, pooling their subjects into
one population (no study-level random effect — a documented simplification).

`SheddingFit` is a dataclass carrying `model`, `method`, `population_mean`,
`population_cov`, `sigma`, `subject_params`, `censoring_limit`, `reference_event`,
`unit`, `specimen`, `biomarker`, `dataset_ids`, `converged`, and the exclusion
counts (`n_excluded_subjects`, `n_dropped_measurements`). It has
`to_dict()` / `from_dict()` so a fit can be saved to YAML and reloaded — the
property that lets an ABM fit once and simulate across many runs without
refitting.

## Errors and validation

`ValueError` for: missing dataset keys; Ct or mixed units; unknown `model`; no
subject meeting `min_observations`; `n_individuals < 1`; an `incubation_period`
array whose length differs from `n_individuals`; a non-positive-definite `Sigma`
(possible when very few subjects survive — the message says so and suggests
pooling datasets).

`UserWarning` for: dropped qualitative/unknown-time/non-positive-time
measurements; excluded subjects; censoring-limit fallback; optimizer
non-convergence (fit still returned, `converged=False` on the object).

## Testing

New file `tests/test_simulate.py`, following the existing test modules.

- **Parameter recovery.** Simulate from known `mu`/`Sigma`/`sigma`, fit, assert
  recovery within tolerance. Run for both models. This is the core correctness
  test.
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
- **Reproducibility.** Same seed produces identical output; different seeds do
  not.
- **Incubation shift.** Scalar, per-individual array, and callable forms each
  shift the curve as expected; `None` leaves output in reference-event time.
- **Edge cases.** Ct dataset rejected; subject with all-censored measurements;
  subjects below `min_observations`; gamma with `t <= 0`; qualitative values;
  `to_dict`/`from_dict` round-trip.
- **Plot** returns a `Figure` and closes it, matching the convention in
  `shedding_peak.py`.

## Limitations

Recorded in module docstrings so users encounter them:

1. Two-stage fitting does not shrink individual estimates, so between-subject
   variance is overestimated (see Fitting above).
2. Point estimates only — no parameter uncertainty propagates into simulations.
   Cohort spread reflects between-individual variation, not estimation
   uncertainty.
3. Both models assume shedding begins at the reference event. Datasets with
   substantial pre-onset sampling (`kissler2021densely`, `kissler2021viral`) are
   poorly served; the Teunis onset-offset model is the future answer.
4. Pooled fits ignore between-study heterogeneity.
5. The exponential model cannot represent a rise and will mis-fit datasets
   sampled from before peak.

## Files

- `shedding_hub/simulate.py` — new module: models, fitting, simulation, plotting.
- `shedding_hub/__init__.py` — export `fit_shedding_model`, `simulate_shedding`,
  `plot_simulated_shedding`, `SheddingFit`.
- `pyproject.toml` — add `scipy` and `numpy` to `dependencies`. `numpy` is
  currently only transitive via pandas/matplotlib but is used directly.
- `tests/test_simulate.py` — new tests.
- `README.md` — a short section under "Analyzing the Data" with a doctest-safe
  example (README doctests run in CI).
