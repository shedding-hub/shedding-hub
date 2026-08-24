# Shedding Hub { .sh-visually-hidden }

<div class="sh-hero" markdown>
![Shedding Hub](assets/shedding-hub-logo.png#only-light)
![Shedding Hub](assets/shedding-hub-logo-light.png#only-dark)
</div>

<p class="sh-tagline" markdown>
A Python package to <strong>load</strong>, <strong>visualise</strong> and
<strong>model</strong> biomarker shedding data — and to <strong>simulate</strong>
shedding kinetics from parameters estimated on real studies.
</p>

Biomarkers — viral RNA, drug metabolites, antibody titres — are shed into
stool, urine, saliva and respiratory specimens on a time course that decides
what testing and surveillance can actually detect. How long someone stays
positive, which specimen to sample and when, how much virus reaches a sewer:
these are questions about shedding kinetics, and the answers are scattered
across thousands of papers, buried in tables, figures and supplementary
files, in units and reference events that no two studies agree on.

Shedding Hub collects them. Every dataset is curated from a peer-reviewed
study into validated YAML against a shared JSON schema, reviewed by at least
two people in the open on GitHub, and re-derived from its extraction script on
every commit, so the numbers cannot silently drift from the paper they came
from. The result is a FAIR repository — findable, accessible, interoperable,
reusable — and this package is how you use it from Python, without writing the
preprocessing yourself.

```bash
pip install shedding-hub
```

Python 3.10 or newer is required.

## What it does

=== "Load"

    Datasets load by their `[author][year][word]` identifier, pinned to a
    commit or a pull request when you need an analysis to stay reproducible.

    ```python
    import shedding_hub as sh

    data = sh.load_dataset("woelfel2020virological")
    ```

=== "Visualise"

    Time courses, heatmaps and per-participant panels, returned as matplotlib
    figures you can restyle or save.

    ```python
    fig = sh.plot_time_course(data, specimen="stool")
    ```

=== "Model"

    Censored maximum-likelihood fits of rise-and-decay curves, handling the
    non-detects properly instead of dropping them or substituting a limit.

    ```python
    fit = sh.fit_shedding_model(data, analyte="stool")
    ```

=== "Simulate"

    Draw shedding trajectories for synthetic individuals from those fitted
    parameters — the input an agent-based or wastewater model needs.

    ```python
    source = sh.shedding_for("SARS-CoV-2", "stool")
    sim = sh.simulate_shedding(source, n_individuals=500, times=range(0, 31))
    ```

## Fitted parameters, already computed

You do not have to fit anything to start simulating. The package ships a
catalog of **159 converged fits over 53 studies**, spanning 10 biomarkers and 18
specimen types, each one a censored ML fit with its sample size, censoring rate
and diagnostics attached.

```python
catalog = sh.load_shedding_catalog()
```

Where several studies measure the same biomarker in the same specimen —
usually with different units, reference events and populations —
[`shedding_for`](reference/selection.md) picks a defensible default and tells
you why, or pools them into a cross-study ensemble. You can always override it.

## What people use it for

The simulation path draws synthetic trajectories by Monte Carlo from the
fitted parameter distributions, carrying inter-individual variability rather
than tracing one average curve. That makes the fits usable as inputs to new
work, not just as descriptions of past studies:

- **Diagnostic testing strategy** — which specimen to sample, and how long after exposure a test can still find something.
- **Surveillance design** — case-based and wastewater-based, where the shed load per infected person sets what a catchment signal means.
- **Disease control measures** — how long isolation or quarantine has to run to cover the infectious period.
- **Transmission and wastewater models** — per-agent shedding trajectories, the input an agent-based model consumes directly.
- **Study planning** — scenario analysis, power calculations and sensitivity analyses against parameter distributions taken from the literature rather than from one convenient study.

## Where to go next

- **[Getting started](getting-started.md)** — load a dataset, summarise it, plot it.
- **[Tutorial](tutorial.md)** — simulate shedding for synthetic individuals, end to end.
- **[Modeling methods](modeling-methods.md)** — what the fitted estimates mean, and what they do not support.
- **[Reference](reference/datasets.md)** — every public function, with a worked example.

Read the modeling methods page before you publish anything built on these
fits. Shedding curves are estimated from small, heterogeneous studies with
substantial censoring, and the page is explicit about which comparisons the
parameters support and which they do not.

## Project

The data repository, the extraction scripts and this package all live in one
place, and the datasets are contributed by the research community.

- [github.com/shedding-hub/shedding-hub](https://github.com/shedding-hub/shedding-hub) — data, code, and the contribution workflow
- [pypi.org/project/shedding-hub](https://pypi.org/project/shedding-hub/) — releases
- [shedding-hub.github.io](https://shedding-hub.github.io/) — the project website and dataset browser
