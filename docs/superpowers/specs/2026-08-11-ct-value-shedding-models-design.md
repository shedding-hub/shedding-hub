# Modelling cycle-threshold data — design

**Date:** 2026-08-11
**Status:** approved, not yet implemented

## Problem

`prepare_observations` refuses every cycle-threshold analyte outright, raising
`SheddingDataError` with `reason: "ct_units"`:

> Cycle-threshold values are inversely related to concentration and already on a
> log scale, so neither shedding model applies.

Both clauses of that sentence are true. The conclusion drawn from them is not,
and it costs the repository a third of its evidence:

| | studies | measurements | subjects |
|---|---|---|---|
| concentration | 51 | 53,426 | 27,295 |
| cycle threshold | 27 | 12,537 | 1,005 |

Sixty-nine analytes across 27 studies — SARS-CoV-2 (48), rotavirus vaccine (9),
mpox (5), norovirus (3), influenza (3), rotavirus (1) — are unfittable today.
By study count that is over half of what the concentration side offers.

## Why the existing models already apply

Cycle threshold is an *affine* function of log10 concentration:

```
Ct = α − β · log10(C)
```

where β is the standard-curve slope (≈3.32 cycles per 10-fold at 100% PCR
efficiency) and α is an assay-specific intercept. The shedding models are
written on the log10 scale already, so a gamma curve in concentration *is* a
gamma curve in Ct — reflected and rescaled. What is missing is a transform of
the response, not a new model family.

## Decision: model Ct directly, anchored at a fixed reference

The response variable is

```
CT_REFERENCE = 40.0                   # fixed for every Ct analyte
depth(t) = CT_REFERENCE − Ct(t)       "cycles below the reference"
```

fitted with the existing `gamma`, `gamma_shifted` and `exponential` forms,
unchanged.

Negating Ct is what makes the response increase with viral load, so the curve is
a peak rather than a trough and no sign flip enters the model itself. Adding a
constant offset keeps fitted levels positive, which matters because
`theta_to_params` enforces `c0 > 0` through `exp(θ)`; 40 sits above the observed
Ct median of 31.0 and above 95% of all readings, so fitted peak heights — which
occur at *low* Ct — are comfortably positive.

The reference is a **single constant across all analytes**, not each study's own
detection cutoff. Recorded cutoffs range from 37 to 41, so anchoring per study
would make two studies measuring identical samples report peak heights differing
by up to 4 cycles purely from the anchoring convention. A fixed reference puts
every Ct fit's height on one scale across all 27 studies.

Each study's own cutoff still does the job it should: it sets that analyte's
censoring limit (below). Nothing is lost by not anchoring to it.

### Alternatives rejected

**Convert to pseudo-log10 via an assumed slope**, `(Ct_ref − Ct)/β₀` with
β₀ = 3.32. Recovers `a0` on the true log10-per-day scale, but imports an
efficiency assumption into every fit that almost no study reports and that would
have to be defended in every downstream paper. Since peak time is invariant to
the choice anyway (see below), the assumption buys only the decay rate, and buys
it on credit.

**Fit raw `−Ct`.** Mathematically identical to the chosen approach up to the
offset, but the fitted level would have to be around −31, and the optimizer
parameterizes `c0` as `exp(θ) > 0`. It would require weakening a positivity
constraint that is otherwise doing useful work.

**A native Ct model carrying α and β as parameters.** Not identifiable. Without
a standard curve, a single study's Ct time series cannot separate the assay
intercept from the subject's shedding level.

**Anchoring at each study's own detection cutoff.** Puts non-detects at exactly
zero, which is tidy, but buys that tidiness with the cross-study height
incomparability described above. The censoring machinery resolves a limit per
analyte already, so a varying limit costs nothing while a varying anchor costs
comparability.

## What is comparable across value types

Let the true standard-curve slope be β. Fitting `depth` returns
`a0* = β·a0` and `b0* = β·b0`, so:

| quantity | meaning for a Ct fit | comparable to concentration fits |
|---|---|---|
| `peak_day` = `b0/a0` | days to peak shedding | **yes — exactly** |
| `t0` (gamma_shifted) | onset of shedding | **yes** |
| rise duration `b0/a0` | onset-to-peak interval | **yes** |
| peak height | cycles below reference (= 40 − min Ct) | no — α is assay-specific, but comparable across all Ct studies |
| `a0` | decay in Ct units per day | only after dividing by an assumed β |
| half-life `ln2/a0` | — | no, for the same reason |

The first three rows are ratios, so β cancels identically. This is the whole
scientific case for the feature: **peak time transfers between value types with
no assumption about PCR efficiency at all.** The magnitude parameters do not,
and the API must say so rather than leaving it to a docstring.

A corollary constrains which models are worth running. The exponential model
peaks at `t = 0` by construction, so its peak time carries no information, and
its only other parameter is `a0`, which lands in the assumption-dependent tier.
An exponential fit to Ct data yields nothing assumption-free. The two gamma
variants are where the value is.

## Censoring falls out unchanged

A non-detect means `Ct ≥ cutoff_analyte`, hence

```
depth ≤ CT_REFERENCE − cutoff_analyte
```

which is left-censored at an analyte-specific limit — exactly the shape of the
left-censored normal likelihood already in `shedding_fit.py`, and exactly the
kind of per-analyte value `_resolve_censoring_limit` already resolves. No new
likelihood, and the ~37% of repository measurements recorded as `negative` keep
contributing an inequality rather than being discarded.

The limit sits near zero: exactly zero where the cutoff is 40, negative where an
assay runs to 41, positive where it stops at 37. Nothing in the likelihood cares
which side of zero it falls on.

## Resolving the censoring limit

Every one of the 69 Ct analytes records a detection limit, and most record it
numerically:

- 48 analytes: numeric `limit_of_detection`
- 6 further analytes: non-numeric LOD but numeric `limit_of_quantification`
- 15 analytes: neither is numeric

For the last group, fall back to the maximum observed Ct and warn, following the
pattern `_resolve_censoring_limit` already uses. Observed Ct across the
repository runs 9.5 to 44.3 (median 31.0, 95th percentile 38.9, 99th 39.9),
consistent with the recorded cutoffs of 37–41, so the fallback is sound where it
is needed.

## API changes

`SheddingFit` gains `value_type: Literal["concentration", "ct"]`, the
`ct_reference` the height is measured against, the analyte's resolved
`ct_cutoff`, and per-parameter comparability tags, so the tiering above is
machine-readable. Recording the reference rather than assuming it means a fit
serialized under one convention cannot be silently misread under another, and it
lets a height be turned back into a minimum Ct as `ct_reference − height`.

`POPULATION_COORDS` needs a value-type-aware height name — `peak_cycles` rather
than `peak_log10` for Ct fits. That dict exists precisely so a catalog written
under an older convention fails loudly instead of being silently misread, so
this extends an intent already present rather than adding a new one.

`prepare_observations` stops raising `ct_units` and applies the transform
instead. `fit_shedding_model`'s mathematics is untouched.

## Visualisation

`plot_fit_diagnostic` on a Ct analyte draws the Ct axis **inverted**, so low Ct
sits high and the fit reads as a shedding peak like the concentration plots,
with the cutoff drawn as a reference line and non-detects placed on it.

## Verification

The repository can validate the central claim on real data rather than only in
simulation. Seven studies report both a Ct and a concentration analyte for the
same biomarker and specimen, and they pair at the subject and timepoint level:

| study | subjects with both | matched same-timepoint pairs |
|---|---|---|
| `kissler2021viral` | 68 | 225 |
| `teunis2015shedding` | 70 | 161 |
| `gutierrez2021nosocomial` | 4 | 4 |

`kissler2021viral` records 2,406 timepoints carrying both analytes, but only
**225** of them have both analytes numerically detected. The rest are
qualitative `negative` strings on one or both, which cannot be regressed. The
figure that matters for verification is 225.

plus `cdc2024nhphrn`, `kim2020viral`, `rouphael2025effective` and
`shetty2024influenza`. The verification plan is therefore:

1. **Peak-time agreement.** Fit both analytes of a paired study and check that
   the peak times agree within uncertainty. This tests the invariance claim
   directly, and it is the claim the feature rests on.
2. **Empirical β.** Regress Ct on log10 concentration over the 225 numerically
   detected `kissler2021viral` pairs to estimate the study's effective slope,
   then check that the fitted `a0*/a0` ratio matches it. This tests the scaling
   claim quantitatively.
3. **Sign-flip guard.** A dedicated test that a Ct fit recovers a peak, not a
   trough. Omitting the flip produces an inverted curve that converges happily
   and is entirely wrong, which is the one failure mode here that would not
   announce itself.
4. **Parameter recovery** on simulated Ct data with known parameters.

## Out of scope

Ensembles and the catalog must **refuse to mix value types** until someone
decides what averaging heights across them would mean; today it would not be
meaningful. This design adds the refusal as a guard and leaves the question
open.

Also out of scope: converting Ct fits to absolute concentration (needs a
standard curve the studies do not report), and adding standard-curve fields to
the data schema.
