# Benchmark analysis plan

**Status: DRAFT. Not yet confirmed and not tagged.** This file becomes the pre-specified plan when Yuke Wang confirms the items in section 1 and the commit is tagged `benchmark-plan-v1`. No output of the full benchmark run is compared with released data before that tag exists. Any change after the tag is a documented deviation.

Draft written 8 Oct 2026. Reference data: Shedding Hub release `v1.1.0` (commit `1fe0e92`, DOI `10.5281/zenodo.23088307`).

## 1. Items to confirm before tagging

Each has a proposed default. The pilot (nine study runs, used for debugging only) informed items 4 to 8.

| # | Item | Proposed default |
|---|---|---|
| 1 | Primary outcome | Measurement-level agreement of S0 on the 37 manual-era studies (section 3) |
| 2 | Tolerances | The values in section 6 |
| 3 | Which tolerance set a study gets | By `input_csv_origin` in the manifest: `figure_digitized` gets the figure tolerances, every other origin gets the reported-value tolerances |
| 4 | `gene_target` | Listed, not scored as an error (free text) |
| 5 | Attributes only the draft records | Listed, not scored as an error |
| 6 | Attributes only the reference records | Scored, as `missing attribute` |
| 7 | Findings that hit no discrepancy | Adjudicated, not counted as false alarms by default |
| 8 | Curator decisions that the sources cannot settle | An adjudication outcome of their own; left out of error counts and reported separately |
| 9 | Bootstrap | 2,000 replicates, resampling whole studies, seed 20261008 |

## 2. Definitions

- **S0:** the extraction agent's draft for a study, never edited by a person.
- **F:** the review agent's findings on S0.
- **S2:** the dataset in release `v1.1.0`. For `obara2008single` the reference is the archived hand-curated file; the manifest's `reference_yaml` column records this.
- **Truth:** S2 corrected by adjudication (section 8). Metrics are computed against truth.
- **Analysis set:** 144 studies in `benchmark/manifest.csv`: 37 manual-era and 107 AI-era (16 pipeline v1, 91 pipeline v2). Two datasets in the release are excluded: `jones2021estimating` (population-scale routine testing) and `cdc2024nhphrn` (no article).
- **Pipeline under test:** tag `pipeline-v2-benchmark` in the private agents repository, with `gpt-5.2-2025-12-11` for extraction and `claude-opus-4-8` for review. Prompts and settings are those of production.

## 3. Primary outcome

**Measurement-level agreement of S0 with truth on the 37 manual-era studies (Analysis A).** It is the share of truth measurements that S0 reproduces: the same participant, the same analyte, time within tolerance, and value within tolerance or the same qualitative result.

It is reported with a 95% confidence interval from a cluster bootstrap over studies, since measurements within a study are not independent.

**What this outcome does and does not measure.** The input CSV for each manual-era study was built for the benchmark from the curator's raw files, with the curator's non-detect rule and value conversions applied by code (`benchmark/manual_curation.py`), as production curators did for the AI era. The outcome therefore measures how well the agent structures a curator-prepared file. The draft and the reference share the curator's decisions and any error in the curator's own file, so the outcome is an upper bound on accuracy against the article.

## 4. Analyses

### A. Manual-era studies (37; primary)

The reference for these studies was curated by hand with no AI draft.

- **Analyte fields:** agreement on biomarker, specimen, unit, reference event, limit of detection and limit of quantification.
- **Measurements:** recall (truth measurements found in S0), precision (S0 measurements that exist in truth), value error for matched pairs (absolute log10 difference for concentrations, absolute difference in cycles for Ct), and agreement on non-detect status.
- **Counts:** agreement on the number of participants and of analytes.
- **Participant attributes:** agreement on attributes the reference records.
- **Strata, reported for every metric:** `input_csv_origin` (14 data repository, 14 figure-digitized, 9 table or supplement). Repository-sourced studies are reported separately, since structuring a clean file is an easier task.
- **Caveat to report:** 33 of the 37 studies are SARS-CoV-2.

### B. AI-era studies (107; secondary)

The same metrics as A, reported as **agreement**, not accuracy. The reference for these studies began as an earlier AI draft that people corrected, and 55 of the 102 production input CSVs were edited after they were first committed, some during review. Agreement is therefore likely to be optimistic.

- **By pipeline version:** 16 v1 and 91 v2, not pooled without saying so.
- **By production extraction model:** descriptive only (19 of the 107 references started from a `gpt-5-mini` draft).
- **By `input_csv_origin`:** 54 figure-digitized, 41 table or supplement, 5 data repository, 2 author-shared, 5 with no input CSV.

### Review agent as an error detector (A and B)

- **Error:** an adjudicated difference between S0 and truth that is scored (section 5).
- **Hit:** a finding whose path is the error's place in S0 or an ancestor of it, and whose check can speak about that type of error (section 7).
- **Sensitivity:** errors with at least one hitting finding, over all errors, by error type.
- **Precision:** findings that hit an error or are judged correct at adjudication, over all scored findings. Reported with and without low-severity findings.
- **Workload:** findings per study and per 100 measurements.
- Schema-gap findings and `model.observation` notes are counted and described, not scored.

### C. Curation throughput from git history

Studies and measurements per month by era, batch sizes, distinct contributors, and the rate of edits within 90 days of merge, excluding commits that touch more than 10 dataset folders. Reported as throughput, not person-time.

### E. Cost and time

Tokens, cost in US dollars and wall-clock time per study and per 100 measurements, from the run logs. Median and interquartile range.

### F. Downstream equivalence (manual-era studies)

Fit the package's curve families to S0 and to truth in one pinned environment, and compare the population summaries. Agreement only, with no biological interpretation.

### Reproducibility

Twenty studies (`replicates.csv`, seed 20261002, stratified by era and input type) are each run three times from the same frozen evidence. Reported: between-run agreement on fields, Jaccard similarity of matched measurement sets, the standard deviation of per-study agreement, and the overlap of finding locations.

## 5. Alignment and discrepancies

Implemented in `benchmark/align.py` and `benchmark/discrepancies.py`, and tested on perturbed copies of released datasets (`tests/test_benchmark_compare.py`).

1. **Analytes** are matched on biomarker, specimen, unit and reference event, never on the key. Two analytes may be paired when biomarker and specimen both agree, or when one of the two differs while unit and reference event agree.
2. **Participants** are matched with the Hungarian algorithm, maximizing shared (analyte, time, value) triples, then shared times, then attributes. Two participants may be paired if they share at least one triple, or if at least half of the measurements of the larger one agree on time.
3. **Measurements** within a matched participant and analyte are matched on time, then on the nearest value, each used once. Leftover measurements whose values agree are paired as the same reading at a different time.

**Types of discrepancy**

| Level | Types |
|---|---|
| Study | missing analyte, extra analyte, participant-count mismatch, no draft, unreadable draft |
| Analyte field | wrong biomarker, specimen, unit, reference_event, limit_of_detection, limit_of_quantification |
| Participant | missing, extra, wrong attribute value, attribute in the wrong field, missing attribute |
| Measurement | omission, extra, wrong time, wrong value, wrong non-detect status |

**Listed but not scored as errors:** a `gene_target` that differs (the field is free text), and an attribute that only the draft records (the reference often leaves attributes out).

A failed or unreadable extraction is a result. It counts every truth measurement of that study as not reproduced.

## 6. Tolerances

| Quantity | Figure-digitized input | Any other input |
|---|---|---|
| Time | within 0.5 day | exact |
| Concentration | within 0.1 log10 (about 26%) | relative difference of at most 0.5% |
| Cycle threshold | within 0.5 cycle | exact |
| Limit of detection or quantification | equal | equal |

"Exact" allows for floating-point rounding only.

## 7. Matching findings to discrepancies

Implemented in `benchmark/match_findings.py`.

| Check | Can hit |
|---|---|
| `analyte.value_vs_paper` | a wrong analyte field; a wrong value, non-detect status or time when the path points at a measurement |
| `attr.value_vs_paper`, `attr.misassigned_vs_paper` | attribute errors |
| `analyte.missing_vs_paper` | missing analyte |
| `participant.count_vs_paper` | participant-count mismatch, missing or extra participant |
| `sample.count_vs_paper`, `sample.types_vs_paper` | measurement omissions and extras |
| `value.censored_string`, `value.nonpositive`, `value.nan` | wrong value or non-detect status |
| `schema.*` | not scored (schema gaps) |
| `model.observation` | not scored (notes) |

A finding cannot hit a difference that is not scored.

## 8. Adjudication

- **Discrepancies:** every scored discrepancy in the 37 manual-era studies, or a random sample of 300 if there are more, plus a random sample of 200 from the AI era.
- **Each item shows** the study, where to look in the source, and the two candidate values labelled A and B in random order. It does not show which came from the agent, or any finding.
- **Outcomes:** A correct, B correct, both acceptable, both wrong (with the correct value), or curator decision that the sources cannot settle. The last is left out of the error counts and reported separately.
- **Unmatched findings:** every scored finding that hits no discrepancy in the manual-era studies, and a random sample of 200 from the AI era, are judged correct, incorrect or cannot tell against the source. The comparison cannot see an error that the draft and the reference share, so such a finding is not a false alarm by default.
- **Agreement:** a second adjudicator rechecks a 10% sample, and Cohen's kappa is reported.

## 9. Limitations to state in the paper

- The models may have seen the public datasets during training.
- The extraction prompt names one benchmark study (`kimse2020viral`).
- Manual-era inputs carry the curator's decisions, so agreement there is an upper bound (section 3).
- AI-era inputs are the final production files and may include corrections made during review (Analysis B).
- For three manual-era studies the curator's script does something the input does not reproduce (`hakki2022onset`, `kimse2020viral`, `lavezzo2020suppression`); `benchmark/manual_inputs/NOTES.csv` says what.
- The pipeline is not public. The plan, the inputs, the comparison code and the run outputs are.
