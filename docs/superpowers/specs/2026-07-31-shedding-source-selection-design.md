# Choosing what to simulate from — design

**Date:** 2026-07-31
**Status:** Draft (awaiting review)

## Problem

The shipped catalog holds 126 fits over 40 datasets. The intended user journey is
*browse a table → pick one study or an ensemble → simulate*, and the browsing
step is where it stalls. An agent-based modeller who wants SARS-CoV-2 shedding in
stool has to name five keys — biomarker, specimen, reference event, unit and
model — before `build_ensemble` will return anything, and naming any one of them
wrongly yields `No fits match`. Nothing tells them what the admissible
combinations are.

The obvious framing of this — "there are many datasets per biomarker, how do I
combine them?" — turns out to be the wrong one. `build_ensemble` already
combines fits, by mixture or by moment matching, and it pools only within a group
agreeing on all five keys. Measured over the shipped catalog, those keys cut 126
fits into **82 groups, 71 of which are a single study**. Only 11 groups hold more
than one study and the largest holds three.

SARS-CoV-2 in stool, the most requested combination, is ten groups:

| reference event | unit | model | studies |
|---|---|---|---|
| symptom onset | gc/mL | exponential | 3 |
| symptom onset | gc/mL | gamma | 2 |
| symptom onset | gc/dry gram | exponential, gamma, gamma_shifted | 1 each |
| enrollment | gc/mL | exponential, gamma, gamma_shifted | 1 each |
| confirmation date | gc/dry gram | exponential, gamma | 1 each |

So the user's real difficulty is not combination. It is that "SARS-CoV-2 in
stool" names a choice across three incommensurable axes — unit (gc/mL and
gc/dry gram are different physical quantities), reference event (different time
origins), and model (fitted to different observation sets, hence not
AIC-comparable, as `docs/modeling-methods.md` sets out) — and the library offers
no help making it and no way to see it.

## Decision: make the fragmentation navigable, do not paper over it

The tempting alternative is to attack the fragmentation: pool across reference
events with per-event offsets, and across units with declared conversions. That
would turn 71 singletons into far fewer, larger groups. It is deliberately **not**
done here. Every such pooled estimate would rest on a conversion factor or an
offset that the catalog does not contain and that would have to be defended per
study; the result would look better supported than it is. Cross-unit and
cross-event pooling can be argued on its own merits in a later spec.

What is added instead is discovery, a documented default, and one correctness
fix. No new statistical assumption, no change to any fitted number, and no change
to the existing ensemble compatibility rules.

## Reference-event classes

Seven reference events appear in the catalog, and they are not the same kind of
thing:

| class | events | fits | offset from infection |
|---|---|---|---|
| `exposure` | inoculation, vaccination | 18 | zero — the event *is* the exposure |
| `landmark` | symptom onset | 61 | the incubation period; defined, and in the literature |
| `administrative` | enrollment, confirmation date, hospital admission, treatment | 47 | none — reflects testing behaviour and health-system access |

This distinction is load-bearing for agent-based use. `simulate_shedding` takes
an `incubation_period` and shifts the time origin back to infection, and today it
accepts all seven events identically and stamps `time_origin="infection"` on the
result. Shifting an enrollment-anchored fit by five days therefore asserts
something untrue: when a subject enrolled in a study has no fixed relation to
when they were infected.

An unrecognised event classifies as `administrative`, the conservative default —
a new dataset must not crash the picker, but neither should it silently acquire a
defensible clock it has not earned.

## Ranking

One rule, applied in strict order:

1. **Reference-event class** — `exposure` and `landmark` before `administrative`.
2. **Model** — `gamma_shifted`, then `gamma`, then `exponential`.
3. **Evidence** — studies, then subjects, then measurements.
4. **The sorted key tuple**, so ties resolve deterministically and a given
   catalog always yields the same pick.

Rule 2 prefers the rise-capable models because for a wastewater model the
pre-symptomatic rise is the epidemiologically interesting part, and an
exponential asserts by construction that an agent sheds maximally on the day of
the reference event. A model appears in the catalog only if its gates passed, so
presence is already the identifiability signal; the ranking needs no separate
check.

Rule 1 outranks rule 2 because a defensible clock beats a better curve shape. A
three-study exponential on symptom onset is preferred to a one-study
`gamma_shifted` on enrollment: the latter resolves a rise, but measures it from a
date with no fixed relation to infection, so the shape it recovers cannot be
placed in an agent's timeline.

Rule 3 comes last of the substantive rules, which is the deliberate cost of this
design: the default will often pass over a better-supported fit for a
worse-supported one with a better clock or curve. `shedding_options` exists so
that trade is visible rather than hidden, and any key can be pinned to override
it.

## API

```python
sh.shedding_options(biomarker="SARS-CoV-2", specimen="stool")
```

Returns a `DataFrame`, one row per compatible group, sorted best first:
`reference_event`, `event_class`, `unit`, `model`, `n_studies`, `n_subjects`,
`n_measurements`, `rank`. Accepts any subset of the five keys, and a `catalog=`
argument defaulting to the shipped one. Rank 1 is by construction what
`shedding_for` returns.

```python
source = sh.shedding_for("SARS-CoV-2", "stool")           # or model=, unit=, reference_event=
traj = sh.simulate_shedding(source, n_individuals=1000,
                            times=np.arange(0, 30), incubation_period=5.0, seed=42)
```

`shedding_for` calls `shedding_options` and takes rank 1, so there is a single
ranking implementation rather than two that can drift apart. It returns a
`SheddingEnsemble` in every case, single-component when one study matched;
`make_ensemble` already guarantees a one-component ensemble consumes the
generator exactly as the underlying fit does, so callers get one type and one
code path regardless of how many studies backed the answer.

The choice is exposed as structured data on `source.selection` — what was picked,
what was passed over, and which rule decided it — rather than only printed, so it
can be asserted on in a test and recorded in a model run's provenance. Passing
keys that match nothing raises, listing what is available, following
`build_ensemble`'s existing pattern.

## The time-origin fix

`simulate_shedding` gains no new arguments. Its `attrs` change:

- `time_origin` is `"infection"` only when the event is a `landmark` and an
  incubation period was applied.
- Applying an incubation period to an `exposure` event warns — the event already
  is the exposure, so there is nothing to bridge — and records
  `f"{event}_shifted"`.
- Applying one to an `administrative` event warns that the origin is not
  infection, and records `f"{event}_shifted"`.
- With no incubation period the origin remains the event itself, as now.
- `attrs` gains `reference_event_class`.

Nothing is forbidden. An administrative fit remains simulable and shiftable; the
false claim merely stops being silent.

## Testing

- Every `reference_event` in the shipped catalog is classified deliberately —
  the same shape of guard as `test_shipped_catalog_covers_every_dataset`, so a
  new event is a decision rather than a silent fallback.
- The SARS-CoV-2 stool ranking is pinned explicitly, including that it picks the
  two-study `gamma` on symptom onset over the three-study `exponential`, and over
  every administrative group.
- `shedding_options(...).iloc[0]` and `shedding_for(...)` agree, on the shipped
  catalog and on synthetic catalogs — the property that keeps the two surfaces
  honest.
- Each of `model`, `unit` and `reference_event` pins correctly and ranks within
  the remainder.
- Ranking is deterministic across repeated calls and independent of catalog fit
  order.
- A single-component result simulates identically to the bare fit for a fixed
  seed.
- Warnings fire for `administrative` + incubation and `exposure` + incubation,
  and `time_origin` records the shifted form rather than `"infection"`.
- No match raises, and the message names the available combinations.

## Out of scope

Unit conversion. Cross-reference-event and cross-unit pooling. Any change to
existing fits, catalog contents, published parameters, or the compatibility keys
`make_ensemble` enforces. Automatic model selection by fit statistic, which
`docs/modeling-methods.md` argues against and which this design deliberately
replaces with a documented, overridable preference order.
