# Catalog fit plots — design

**Date:** 2026-07-27
**Status:** implemented 2026-07-29, with three amendments recorded below

## Problem

The shipped catalog holds 83 fitted shedding curves — 29 studies, 8 biomarkers,
10 specimens, both models where each is admissible. There is currently no way to
see them. `catalog.table` lists their parameters as numbers, and
`plot_simulated_shedding` draws one fit's simulated cohort, but nothing renders
the fitted curves themselves side by side, which is what a modeller choosing
between studies actually wants to look at.

The goal is one function that draws the median individual of every fit in the
catalog, grouped so that curves sharing a panel are genuinely comparable.

## The comparability constraint

Curves cannot be overlaid across studies that disagree on either axis.

The catalog spans five units — `gc/mL`, `gc/dry gram`, `gc/wet gram`, `gc/swab`,
`pfu/mL` — so the y axis is not shared across all of it. It also spans five
reference events — symptom onset, enrollment, confirmation date, vaccination,
hospital admission — so `t = 0` does not mean the same thing either. Overlaying
across those would produce a picture that reads as a comparison while comparing
nothing.

The grouping key is therefore `(biomarker, specimen, unit, reference_event)`.
That yields 31 groups, of which 7 hold more than one study:

| biomarker | specimen | unit | reference event | studies | fits |
|---|---|---|---|---|---|
| SARS-CoV-2 | stool | gc/mL | symptom onset | 3 | 6 |
| rotavirus vaccine | stool | gc/wet gram | vaccination | 2 | 5 |
| SARS-CoV-2 | nasopharyngeal_swab | gc/mL | symptom onset | 3 | 5 |
| SARS-CoV-2 | oropharyngeal_swab | gc/mL | symptom onset | 3 | 4 |
| SARS-CoV-2 | nasopharyngeal_swab | gc/swab | symptom onset | 2 | 3 |
| norovirus | stool | gc/wet gram | symptom onset | 3 | 3 |
| SARS-CoV-2 | sputum | gc/mL | symptom onset | 2 | 2 |

The remaining 24 groups hold a single study each (55 fits). They are still drawn:
a one-study panel is a legitimate view of that study, and dropping them would
silently hide two thirds of the catalog. The layout simply makes it obvious which
panels support a cross-study read and which do not.

A coarser grouping (biomarker × specimen, 17 panels) was considered and rejected:
it buys more overlap per panel by mixing units and time origins on shared axes,
which is precisely the misleading picture above.

## API

Added to `shedding_hub/viz.py`, exported from `shedding_hub`:

```python
def plot_catalog_fits(
    catalog: SheddingCatalog | list[SheddingFit],
    *,
    biomarker: str | None = None,
    specimen: str | None = None,
    unit: str | None = None,
    reference_event: str | None = None,
    dataset_ids: Sequence[str] | None = None,
    n_days: float | None = None,
    show_extrapolation: bool = True,
    ncols: int = 3,
    figsize_per_panel: tuple[float, float] = (4.5, 3.2),
) -> Figure
```

Accepting either a `SheddingCatalog` or a bare list of `SheddingFit` keeps it
usable on a fresh fit of private data, matching how `make_ensemble` accepts fits
from anywhere. The filter arguments reuse the names already used by
`plot_time_courses` and `plot_mean_trajectory`.

## Rendering

**One panel per group**, arranged in a grid `ncols` wide, titled with the group
key. Panels are ordered by descending study count, so the cross-study panels come
first.

**One colour per study** within a panel, cycling `TABLEAU_COLORS` as the rest of
`viz.py` does. **Solid line for the exponential model, dashed for gamma**, so a
study fitted with both contributes two lines of one colour and the models can be
compared at a glance. Where the gamma model was refused — 38 analytes, mostly for
showing no rise — there is simply one line, and the absence is itself
informative.

**The curve** is `log10_concentration(model, fit.median_params, t)`: the median
individual, the same object the notebook draws in red. Not the simulated cohort
median, which is a different quantity and would need sampling per panel.

**Each study's censoring limit** is drawn as a faint horizontal line in that
study's colour. Limits differ between studies, so one shared line would be wrong.

**Extrapolation is shown by opacity, not linestyle.** Each curve is drawn in two
segments, split at that study's `median_first_observed_day`: before it at
`alpha=0.35`, after it at full opacity. The faded portion is the stretch no data
constrains. This matters more than it sounds — 11 of the 22 gamma fits have their
fitted peak *earlier* than the median first observation, so for half of them the
entire rise phase is functional form rather than measurement. Woelfel stool is the
clearest case: peak at day 1.2, first observation at day 6. Setting
`show_extrapolation=False` draws every curve at full opacity.

The linestyle channel is spent on the model and the colour channel on the study,
which is why extrapolation gets opacity. It is the only channel left, and it
carries the right connotation: faded means less certain.

**X horizon.** A fit does not record its last observed day, so with `n_days=None`
the horizon is derived per panel as `max(peak_day + 5 × half_life_days)` over the
panel's fits, clamped to `[7, 60]` days. Five half-lives puts the curve about 1.5
log10 below its peak, which is far enough to show the decline and near enough to
avoid a long flat tail. The clamp guards both ends: a runaway `half_life_days`
would otherwise stretch the axis uselessly, and a very fast decay would give a
panel two days wide. Passing `n_days` overrides the derivation for every panel,
which is also how to put panels on a common axis.

Curves start at a small positive epsilon (0.05 days) rather than zero, because
the gamma model is undefined at `t = 0`: `c(t) = c0·t^b0·e^(−a0·t)` sends the
concentration to zero as `t → 0`, so `log10` diverges. The exponential model is
defined there and simply starts at 0.05 too, for one code path.

## Errors

Filters that match no fit raise `ValueError` naming the combinations that do
exist, following `build_ensemble`'s precedent — a silently empty figure is worse
than a refusal. An empty catalog raises likewise.

## Testing

Built with TDD, in `tests/test_viz.py` alongside the other plot tests, using
synthetic `SheddingFit` objects rather than the shipped catalog so the tests do
not move when the catalog is rebuilt.

- Fits differing only in `unit` are drawn in separate panels — the comparability
  rule, and the reason the grouping key is what it is.
- Fits differing only in `reference_event` are likewise separated.
- A study fitted with both models yields two lines of one colour with different
  linestyles.
- Each filter argument reduces the panel count as expected, and `dataset_ids`
  restricts to the named studies.
- With `show_extrapolation=True`, a fit whose `median_first_observed_day` falls
  inside the drawn range produces both a faded and a full-opacity segment; with
  it `False`, every segment is full opacity.
- `n_days` overrides the derived horizon; the derived horizon respects the clamp
  at both ends.
- No-match and empty-catalog inputs raise `ValueError`.
- Panels are ordered with the multi-study groups first.
- The figure returned has one axis per group and is closed in the pyplot state,
  matching the convention in `shedding_peak.py` and `shedding_simulate.py`.

## Amendments found during implementation

Three things the design above did not anticipate, each caught by drawing the
actual catalog rather than the synthetic fits the tests use.

**A study can contribute several analytes to one panel.** The encoding above
spends colour on the study and linestyle on the model, which assumes at most one
fit per (study, model) per panel. Eight (study, panel, model) combinations in the
catalog break that assumption, and `natarajan2022gastrointestinal` contributes
*fourteen* stool assays to a single panel — gene target × assay chemistry × lab,
with median peaks spanning 3.26 to 5.43 log10. Labelling them all
`natarajan2022gastrointestinal (exponential)` produced fourteen identical legend
keys for fourteen materially different curves, which is a legend that lies.

The analyte is therefore named in the label — `study analyte (model)` — but only
for studies contributing more than one analyte to that panel, so the common case
keeps a short label and both models of one analyte stay distinguished by
linestyle alone. Colour still means study, as designed: within-study spread reads
as a band of one colour, which is the honest picture of assay disagreement.

**A crowded legend buries its own panel.** Seventeen entries at 7pt in a
4.5×3.2 inch panel covered the curves entirely. The legend is capped at six
entries with a `+ N more` line, so a truncated legend never reads as the whole
panel.

**The derived horizon can waste most of the y axis.** Deriving the horizon as the
`max` over a panel's fits means a slow-decaying fit stretches the axis while a
fast one plunges over that span — to −9.1 log10 gc/mL in one panel, and −30 in a
synthetic worst case. Ten of the 31 panels gave more than a quarter of their
height to concentrations below every censoring limit in the panel; one gave 85%.

The y axis is therefore floored just below the panel's lowest censoring limit,
since nothing under it was measurable by any study there. The floor also clears
the weakest curve's own peak: two catalog fits (`hakki2022onset`
asymptomatic_cultivable, `kissler2021viral` AN_OPS_SARSCoV2_viral) have median
individuals peaking *below* their own limit, and flooring at the limit alone
would have dropped them off the bottom of the panel.

## Out of scope

**No CSV export.** `load_shedding_catalog().table.to_csv(path)` already writes
every parameter estimate the table carries. Adding a wrapper would be API surface
for a one-liner. `population_cov` is deliberately absent from the table — it is a
`k × k` matrix per fit and does not belong in a flat browsing surface; reach for
`fit.population_cov` when you need it.

**No refitting.** `scripts/build_shedding_catalog.py` already fits both models
against every analyte in `data/` and records a reason for each of the 207 refused
combinations. This design consumes that output and does not change it.

**No new notebook.** A cell may be added to the existing walkthrough showing the
function on one filtered group, but the walkthrough's subject is the estimation
pipeline for a single study and should stay that.
