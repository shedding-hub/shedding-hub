# Package Documentation Site Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace a hand-written package page that documents 19 of 42 public names — three of which raise `AttributeError` — with a generated, executed, per-function reference on Read the Docs.

**Architecture:** MkDocs Material renders `docs/` and mkdocstrings generates the reference from docstrings, so signatures and prose cannot disagree with the code. Every public name carries a worked example written as a real doctest, run by `pytest --doctest-modules` in CI, because the six examples that already exist were never executed and all six fail. Plotting functions additionally get a figure rendered at build time, since a doctest can only assert `Figure`.

**Tech Stack:** MkDocs Material, mkdocstrings[python], Read the Docs, pytest doctests, matplotlib. No new *install* dependencies — everything added is docs/dev tooling.

## Global Constraints

- No new entries in `pyproject.toml`'s `dependencies`. Docs tooling goes in `requirements.in` (pip-compile source) and a separate `docs/requirements.txt` for Read the Docs.
- `black --check .` must pass. Run `black .` before every commit.
- **Never invent doctest expected output.** Write the example, run it, paste what it actually printed, run it again. An invented value that happens to be wrong is worse than no example.
- **Never let a numpy scalar be a doctest's expected value.** Indexing a DataFrame
  returns `np.int64`/`np.float64`, whose repr changed in numpy 2: `4` under 1.x,
  `np.int64(4)` under 2.x. Four examples written against local numpy 1.26 failed
  CI on the pinned 2.3.5. Wrap the lookup in `int(...)` or `float(...)`, which
  reprs identically on both.
- **Do not print DataFrames in doctests.** Their repr depends on terminal width and pandas version. Follow the convention the README already uses: assert on `list(df.columns)`, `df.shape`, or a single indexed value.
- Examples must run **offline and deterministically**. Use `sh.load_shedding_catalog()` (ships inside the package) for the fitting/simulation/selection surface, and `sh.load_dataset(name, local='./data')` for the dataset surface.
- The shipped catalog holds 126 fits over 40 datasets; `sh.MODELS` is `('exponential', 'gamma', 'gamma_shifted')`. Do not restate other counts in docstrings — they go stale.
- `docs/superpowers/**` must never be published. It holds internal design records.
- Tests must not reload the shipped catalog per test; use the session-scoped `shipped_catalog` fixture in `tests/conftest.py`.
- `plot_shedding_duration(s)` and `plot_shedding_peak(s)` take the DataFrame returned by `calc_shedding_duration`/`calc_shedding_peak`, not a dataset. `calc_shedding_durations`/`calc_shedding_peaks` take dataset *ids* and fetch from GitHub, so examples must not call them — they would need the network. The plural *plotters* need `output='summary'`; the singular ones take the default `output='individual'`. All four call shapes are verified.

**Out of scope for this plan:** rewriting `shedding-hub.github.io/package.html`. That file lives in a *different repository* (`shedding-hub/shedding-hub.github.io`) which is not checked out here. It is a follow-up, tracked at the end.

---

### Task 1: Export the three orphaned plotting functions

**Files:**
- Modify: `shedding_hub/__init__.py`
- Test: `tests/test_viz.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `sh.plot_clearance_curve`, `sh.plot_detection_probability`, `sh.plot_value_distribution_by_time` at package level. `__all__` grows from 42 to 45. Later tasks document all 45.

`shedding_hub/viz.py` already defines all three, with 1,648–2,420 character docstrings. None is exported, so the website's `sh.plot_clearance_curve(...)` raises `AttributeError`. None has a single test reference anywhere in `tests/`, against 63 for `plot_time_course`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_viz.py`:

```python
import pytest

from matplotlib.figure import Figure


@pytest.mark.parametrize(
    "name",
    ["plot_clearance_curve", "plot_detection_probability",
     "plot_value_distribution_by_time"],
)
def test_orphaned_plots_are_exported(name):
    """
    All three are implemented and were documented on the project website, but
    never exported, so every call in that documentation raised AttributeError.
    """
    import shedding_hub as sh

    assert hasattr(sh, name), f"sh.{name} is documented but not exported"
    assert name in sh.__all__


def test_plot_clearance_curve_draws(woelfel_dataset):
    import shedding_hub as sh

    fig = sh.plot_clearance_curve(woelfel_dataset, specimen="sputum")
    assert isinstance(fig, Figure)
    assert len(fig.axes) >= 1


def test_plot_detection_probability_draws(woelfel_dataset):
    import shedding_hub as sh

    fig = sh.plot_detection_probability(woelfel_dataset, specimen="sputum")
    assert isinstance(fig, Figure)
    assert len(fig.axes) >= 1


def test_plot_value_distribution_by_time_draws(woelfel_dataset):
    import shedding_hub as sh

    fig = sh.plot_value_distribution_by_time(woelfel_dataset, specimen="sputum")
    assert isinstance(fig, Figure)
    assert len(fig.axes) >= 1
```

- [ ] **Step 2: Add the dataset fixture these tests need**

`tests/test_viz.py` may already load this dataset ad hoc. Add a shared fixture to `tests/conftest.py` so it is loaded once:

```python
@pytest.fixture(scope="session")
def woelfel_dataset():
    """A real dataset for plot smoke tests, read from the repo rather than the network."""
    from shedding_hub import load_dataset

    return load_dataset("woelfel2020virological", local="./data")
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `python -m pytest tests/test_viz.py -q -k "orphaned or clearance_curve_draws or detection_probability_draws or value_distribution_by_time_draws"`
Expected: FAIL — `AttributeError: module 'shedding_hub' has no attribute 'plot_clearance_curve'`

- [ ] **Step 4: Export the three functions**

In `shedding_hub/__init__.py`, find the existing viz import block:

```python
from .viz import (
    plot_time_course,
    plot_time_courses,
    plot_shedding_heatmap,
    plot_mean_trajectory,
    plot_catalog_fits,
    plot_fit_diagnostic,
)
```

Replace with:

```python
from .viz import (
    plot_time_course,
    plot_time_courses,
    plot_shedding_heatmap,
    plot_mean_trajectory,
    plot_catalog_fits,
    plot_fit_diagnostic,
    # Implemented since before 0.1.3 and documented on the project website, but
    # never exported until now, so every documented call raised AttributeError.
    plot_clearance_curve,
    plot_detection_probability,
    plot_value_distribution_by_time,
)
```

Then add the three names to `__all__`, beside the other `plot_` entries.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_viz.py -q -k "orphaned or clearance_curve_draws or detection_probability_draws or value_distribution_by_time_draws"`
Expected: PASS, 6 tests (3 parametrized + 3 draw tests).

- [ ] **Step 6: Confirm the public surface is now 45**

Run:

```bash
python -c "import shedding_hub as sh; print(len(sh.__all__))"
```
Expected: `45`

- [ ] **Step 7: Format and commit**

```bash
black shedding_hub/__init__.py tests/test_viz.py tests/conftest.py
git add shedding_hub/__init__.py tests/test_viz.py tests/conftest.py
git commit -m "fix: export the three plotting functions the website already documents"
```

---

### Task 2: Run module doctests in CI, and repair the six that fail

**Files:**
- Modify: `pyproject.toml` (add a `[tool.pytest.ini_options]` block)
- Modify: `shedding_hub/stats.py` (six docstrings)
- Modify: `.github/workflows/build.yaml`

**Interfaces:**
- Consumes: nothing.
- Produces: `pytest --doctest-modules shedding_hub/` passes and runs in CI. Later tasks add examples that this same command verifies.

Six names in `stats.py` carry an `Examples:` block. All six fail, because they were written as illustration — `>>> print(...)` with no expected output beneath it — and CI has only ever run doctests against `README.md`.

- [ ] **Step 1: Watch all six fail**

Run: `python -m pytest --doctest-modules shedding_hub/stats.py -q`
Expected: `6 failed`, each reporting `Expected nothing / Got: ...`.

- [ ] **Step 2: Repair the six examples**

For each of `calc_shedding_summary`, `calc_detection_summary`, `calc_clearance_summary`, `calc_value_summary`, `calc_dataset_summary` and `compare_datasets` in `shedding_hub/stats.py`, rewrite the `Examples:` block so every `>>>` line that produces output has that output written beneath it.

Do not print DataFrames. Use the README's convention. For example, `calc_dataset_summary`'s block becomes:

```
    Examples:
        >>> import shedding_hub as sh
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> summary = sh.calc_dataset_summary(data)
        >>> sorted(summary)  # doctest: +ELLIPSIS
        [...]
        >>> summary['n_participants']
        9
```

**Determine every expected value by running it, never by guessing.** For each docstring: write the example, run
`python -m pytest --doctest-modules shedding_hub/stats.py -q`, read the `Got:` block, paste that exact text as the expected output, and run again until it passes.

- [ ] **Step 3: Verify all six pass**

Run: `python -m pytest --doctest-modules shedding_hub/stats.py -q`
Expected: `6 passed`.

- [ ] **Step 4: Make doctests part of the default test run**

Add to `pyproject.toml`:

```toml
[tool.pytest.ini_options]
# Docstring examples are documentation users copy, so they are tested like any
# other claim. Six examples predating this setting had never been executed and
# all six were broken.
addopts = "--doctest-modules --doctest-glob='*.md'"
testpaths = ["tests", "shedding_hub"]
doctest_optionflags = "ELLIPSIS NORMALIZE_WHITESPACE"
```

- [ ] **Step 5: Run the whole suite**

Run: `python -m pytest -q`
Expected: PASS. The count rises above 506 as doctests are collected. If a doctest in a module other than `stats.py` now fails, fix that docstring the same way — by running it, not by guessing.

- [ ] **Step 6: Say so in CI**

In `.github/workflows/build.yaml`, the `Run doctests` step currently covers only the README. Rename and extend it so the intent is legible:

```yaml
      - name: Run README doctests
        run: python -m doctest -o ELLIPSIS -o NORMALIZE_WHITESPACE README.md
```

The module doctests now run inside the existing `Run tests` step via `addopts`, so no new step is needed — but confirm by checking that the pytest output in CI reports more tests than before.

- [ ] **Step 7: Format and commit**

```bash
black shedding_hub/stats.py
git add pyproject.toml shedding_hub/stats.py .github/workflows/build.yaml
git commit -m "test: execute docstring examples, and fix the six that never ran"
```

---

### Task 3: Examples for the dataset surface

**Files:**
- Modify: `shedding_hub/util.py` (5 names), `shedding_hub/shedding_duration.py` (4), `shedding_hub/shedding_peak.py` (4)

**Interfaces:**
- Consumes: the doctest runner from Task 2.
- Produces: worked examples on 13 public names. No API change.

Names to cover: `load_dataset`, `check_dataset`, `normalize_str`, `folded_str`, `literal_str`, `calc_shedding_duration`, `calc_shedding_durations`, `plot_shedding_duration`, `plot_shedding_durations`, `calc_shedding_peak`, `calc_shedding_peaks`, `plot_shedding_peak`, `plot_shedding_peaks`.

- [ ] **Step 1: Add an `Examples:` block to each of the 13 docstrings**

Follow the module's existing Google-style layout — `Args:`, `Returns:`, then `Examples:` last. Two patterns, by return type.

For a data-returning function, show the value:

```
    Examples:
        >>> import shedding_hub as sh
        >>> sh.normalize_str('  Sputum Sample ')
        'sputum sample'
```

For a plotting function, assert the type — the picture is supplied separately by Task 6.
**Note the shape**: the duration and peak plotters take the DataFrame produced by
the matching `calc_*` function, *not* a dataset. Passing a dataset raises.

```
    Examples:
        >>> import shedding_hub as sh
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> durations = sh.calc_shedding_duration(data)
        >>> fig = sh.plot_shedding_duration(durations)
        >>> type(fig).__name__
        'Figure'
```

`normalize_str`'s output above is illustrative of the *shape*, not a verified value — run it and use what it prints.

- [ ] **Step 2: Run them, pasting real output**

Run: `python -m pytest --doctest-modules shedding_hub/util.py shedding_hub/shedding_duration.py shedding_hub/shedding_peak.py -q`

Iterate until it passes: read each `Got:` block and paste it as the expected output. Never invent a value.

- [ ] **Step 3: Verify every name in these modules now has an example**

Run:

```bash
python - <<'PY'
import inspect, shedding_hub as sh
mods = ("util", "shedding_duration", "shedding_peak")
missing = [
    n for n in sh.__all__
    if getattr(getattr(sh, n), "__module__", "").split(".")[-1] in mods
    and ">>>" not in (inspect.getdoc(getattr(sh, n)) or "")
]
print("missing examples:", missing)
assert not missing, missing
PY
```
Expected: `missing examples: []`

- [ ] **Step 4: Format and commit**

```bash
black shedding_hub/util.py shedding_hub/shedding_duration.py shedding_hub/shedding_peak.py
git add shedding_hub/util.py shedding_hub/shedding_duration.py shedding_hub/shedding_peak.py
git commit -m "docs: worked examples for the dataset surface"
```

---

### Task 4: Examples for the visualisation surface

**Files:**
- Modify: `shedding_hub/viz.py` (9 names)

**Interfaces:**
- Consumes: Task 1's exports, Task 2's doctest runner.
- Produces: worked examples on all 9 `viz` names.

Names: `plot_time_course`, `plot_time_courses`, `plot_shedding_heatmap`, `plot_mean_trajectory`, `plot_catalog_fits`, `plot_fit_diagnostic`, `plot_clearance_curve`, `plot_detection_probability`, `plot_value_distribution_by_time`.

- [ ] **Step 1: Add an `Examples:` block to each of the 9 docstrings**

Every one returns a `Figure`, so all follow one pattern:

```
    Examples:
        >>> import shedding_hub as sh
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> fig = sh.plot_time_course(data, specimen='sputum')
        >>> type(fig).__name__
        'Figure'
```

`plot_catalog_fits` and `plot_fit_diagnostic` take a fit or catalog rather than a dataset, so theirs use the shipped catalog:

```
    Examples:
        >>> import shedding_hub as sh
        >>> catalog = sh.load_shedding_catalog()
        >>> fit = catalog.select(
        ...     dataset_id='woelfel2020virological', analyte='stool', model='gamma'
        ... )
        >>> data = sh.load_dataset('woelfel2020virological', local='./data')
        >>> fig = sh.plot_fit_diagnostic(fit, data)
        >>> type(fig).__name__
        'Figure'
```

- [ ] **Step 2: Run them**

Run: `python -m pytest --doctest-modules shedding_hub/viz.py -q`
Iterate until it passes, pasting real output.

- [ ] **Step 3: Verify coverage of this module**

```bash
python - <<'PY'
import inspect, shedding_hub as sh
missing = [
    n for n in sh.__all__
    if getattr(getattr(sh, n), "__module__", "").endswith("viz")
    and ">>>" not in (inspect.getdoc(getattr(sh, n)) or "")
]
print("missing examples:", missing)
assert not missing, missing
PY
```
Expected: `missing examples: []`

- [ ] **Step 4: Format and commit**

```bash
black shedding_hub/viz.py
git add shedding_hub/viz.py
git commit -m "docs: worked examples for the visualisation surface"
```

---

### Task 5: Examples for the modelling surface

**Files:**
- Modify: `shedding_hub/shedding_models.py`, `shedding_fit.py`, `shedding_catalog.py`, `shedding_ensemble.py`, `shedding_simulate.py`, `shedding_select.py`

**Interfaces:**
- Consumes: Task 2's doctest runner.
- Produces: worked examples on the remaining public names — the 0.2.0 surface plus the three constants.

Names: `MODELS`, `PARAM_NAMES`, `REFERENCE_EVENT_CLASSES`, `SheddingFit`, `SheddingDataError`, `fit_shedding_model`, `SheddingCatalog`, `fit_shedding_models`, `load_shedding_catalog`, `SheddingEnsemble`, `make_ensemble`, `simulate_shedding`, `plot_simulated_shedding`, `Selection`, `classify_reference_event`, `shedding_for`, `shedding_options`.

Everything here works off the shipped catalog, so no example needs the network or a dataset.

- [ ] **Step 1: Add the examples**

These four are verified — the outputs below were produced by running them, and may be pasted as-is:

```
    >>> import shedding_hub as sh
    >>> sh.classify_reference_event('enrollment')
    'administrative'
```

```
    >>> import shedding_hub as sh
    >>> sh.MODELS
    ('exponential', 'gamma', 'gamma_shifted')
```

```
    >>> import shedding_hub as sh
    >>> sorted(sh.REFERENCE_EVENT_CLASSES)[:3]
    ['confirmation date', 'enrollment', 'hospital admission']
```

```
    >>> import shedding_hub as sh
    >>> catalog = sh.load_shedding_catalog()
    >>> len(catalog.fits)
    126
```

For `shedding_options`, avoid printing the frame:

```
    >>> options = sh.shedding_options('SARS-CoV-2', 'stool', catalog=catalog)
    >>> options.shape
    (10, 11)
    >>> options.iloc[0]['model']
    'gamma'
```

For `simulate_shedding`, follow the README and assert on the columns, and start times at 1 — the gamma model is undefined at `t <= 0`:

```
    >>> source = sh.shedding_for('SARS-CoV-2', 'stool', catalog=catalog)
    >>> traj = sh.simulate_shedding(
    ...     source, n_individuals=10, times=np.arange(1, 8), seed=42
    ... )
    >>> list(traj.columns)
    ['individual_id', 'time', 'log10_value', 'value', 'detected', 'source_dataset_id']
```

Write the rest to the same pattern, running each to obtain its output.

- [ ] **Step 2: Run them**

Run: `python -m pytest --doctest-modules shedding_hub/ -q`
Iterate until it passes.

- [ ] **Step 3: Verify every one of the 45 names now has an example**

```bash
python - <<'PY'
import inspect, shedding_hub as sh
missing = [n for n in sh.__all__ if ">>>" not in (inspect.getdoc(getattr(sh, n)) or "")]
print(f"{len(sh.__all__) - len(missing)}/{len(sh.__all__)} documented")
print("missing:", missing)
assert not missing, missing
PY
```
Expected: `45/45 documented`, `missing: []`

- [ ] **Step 4: Format and commit**

```bash
black shedding_hub/
git add shedding_hub/
git commit -m "docs: worked examples for the modelling surface"
```

---

### Task 6: Render an example figure for every plotting function

**Files:**
- Create: `scripts/build_doc_figures.py`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: Task 1's exports.
- Produces: `python scripts/build_doc_figures.py` writes one PNG per `plot_*` name into `docs/images/`, named `<function>.png`. Task 7's reference pages embed those paths.

A doctest on a plotting function can assert `Figure` and nothing else, which does not help a reader choosing between `plot_shedding_heatmap` and `plot_mean_trajectory`.

- [ ] **Step 1: Write the generator**

Create `scripts/build_doc_figures.py`:

```python
"""
Render one example figure per plotting function, for the documentation.

Run by the docs build, never committed: an image committed once goes stale
silently the moment a plot changes, which is the failure the documentation
design exists to correct. Regenerating means a plot that breaks fails the build.
"""

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import shedding_hub as sh  # noqa: E402

OUTPUT = REPO_ROOT / "docs" / "images"


def _dataset():
    return sh.load_dataset("woelfel2020virological", local=str(REPO_ROOT / "data"))


def _fit(catalog):
    return catalog.select(
        dataset_id="woelfel2020virological", analyte="stool", model="gamma"
    )


def main() -> int:
    import numpy as np

    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = _dataset()
    catalog = sh.load_shedding_catalog()
    fit = _fit(catalog)
    source = sh.shedding_for("SARS-CoV-2", "stool", catalog=catalog)
    traj = sh.simulate_shedding(
        source, n_individuals=200, times=np.arange(1, 31), seed=42
    )

    figures = {
        "plot_time_course": lambda: sh.plot_time_course(data, specimen="sputum"),
        "plot_time_courses": lambda: sh.plot_time_courses([data], specimen="sputum"),
        "plot_shedding_heatmap": lambda: sh.plot_shedding_heatmap(
            data, specimen="sputum", value="concentration"
        ),
        "plot_mean_trajectory": lambda: sh.plot_mean_trajectory(
            data, specimen="sputum", value="concentration"
        ),
        "plot_clearance_curve": lambda: sh.plot_clearance_curve(
            data, specimen="sputum"
        ),
        "plot_detection_probability": lambda: sh.plot_detection_probability(
            data, specimen="sputum"
        ),
        "plot_value_distribution_by_time": lambda: sh.plot_value_distribution_by_time(
            data, specimen="sputum"
        ),
        # These four take the DataFrame the matching calc_* returns, not a
        # dataset. calc_shedding_durations/peaks take dataset *ids* and read
        # from GitHub, so they are driven from the local dataset instead.
        "plot_shedding_duration": lambda: sh.plot_shedding_duration(
            sh.calc_shedding_duration(data)
        ),
        # The plural variants need the *summary* frame -- they read
        # shedding_duration_mean / n_participant, which output='individual'
        # does not produce.
        "plot_shedding_durations": lambda: sh.plot_shedding_durations(
            sh.calc_shedding_duration(data, output="summary")
        ),
        "plot_shedding_peak": lambda: sh.plot_shedding_peak(
            sh.calc_shedding_peak(data)
        ),
        "plot_shedding_peaks": lambda: sh.plot_shedding_peaks(
            sh.calc_shedding_peak(data, output="summary")
        ),
        "plot_fit_diagnostic": lambda: sh.plot_fit_diagnostic(fit, data),
        "plot_catalog_fits": lambda: sh.plot_catalog_fits(catalog),
        "plot_simulated_shedding": lambda: sh.plot_simulated_shedding(
            traj, source=source
        ),
    }

    for name, build in figures.items():
        # No try/except: a plot that cannot be drawn must fail the docs build
        # rather than leave the previous image in place.
        figure = build()
        figure.savefig(OUTPUT / f"{name}.png", dpi=110, bbox_inches="tight")
        print(f"wrote {name}.png")

    print(f"{len(figures)} figure(s) in {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Run it**

Run: `python scripts/build_doc_figures.py`
Expected: 14 lines of `wrote <name>.png`, then `14 figure(s) in ...`.

If a call signature differs from the guess above, correct the lambda — do not wrap it in a `try`. The point of the script is that a broken plot is loud.

- [ ] **Step 3: Confirm one figure per exported plotting function**

```bash
python - <<'PY'
import pathlib, shedding_hub as sh
plots = {n for n in sh.__all__ if n.startswith("plot_")}
made = {p.stem for p in pathlib.Path("docs/images").glob("*.png")}
print("plotting functions:", len(plots), "| images:", len(made))
print("missing images:", sorted(plots - made))
assert not plots - made
PY
```
Expected: `plotting functions: 14 | images: 14`, `missing images: []`

- [ ] **Step 4: Ignore the generated images**

Append to `.gitignore`:

```
# Regenerated by scripts/build_doc_figures.py on every docs build. Committing
# them would let a stale chart outlive the code that drew it.
docs/images/
site/
```

- [ ] **Step 5: Commit**

```bash
black scripts/build_doc_figures.py
git add scripts/build_doc_figures.py .gitignore
git commit -m "docs: render an example figure for every plotting function"
```

---

### Task 7: The MkDocs site and its reference pages

**Files:**
- Create: `mkdocs.yml`, `docs/index.md`, `docs/getting-started.md`, `docs/tutorial.md`, `docs/reference/*.md`
- Test: `tests/test_docs.py`

**Interfaces:**
- Consumes: Tasks 1–6.
- Produces: `mkdocs build --strict` produces `site/`. Task 8 publishes it.

- [ ] **Step 1: Write the coverage test first**

Create `tests/test_docs.py`:

```python
"""
The reference must not drift from the package.

A generated reference cannot disagree with the code about a signature, but it
can silently omit a module nobody added to the nav -- the same rot that left
the project website documenting 19 of 42 names, three of which did not exist.
"""

import pathlib

import shedding_hub as sh

REFERENCE = pathlib.Path(__file__).parent.parent / "docs" / "reference"


def _documented_names() -> set:
    text = "\n".join(p.read_text(encoding="utf-8") for p in REFERENCE.glob("*.md"))
    return {name for name in sh.__all__ if f"shedding_hub.{name}" in text}


def test_every_public_name_appears_in_the_reference():
    missing = sorted(set(sh.__all__) - _documented_names())
    assert not missing, (
        f"{len(missing)} public name(s) are exported but absent from "
        f"docs/reference/: {missing}. Add them to the matching module page."
    )


def test_every_callable_has_a_runnable_example():
    """
    Constants are exempt, and cannot not be.

    ``inspect.getdoc`` on a module-level constant returns its *container type's*
    docstring -- ``tuple.__doc__`` for MODELS, ``dict.__doc__`` for PARAM_NAMES
    and REFERENCE_EVENT_CLASSES -- so no example can ever be attached to the
    name itself. Those three carry executed examples in their defining module's
    docstring instead, and test_every_public_name_appears_in_the_reference still
    holds all 45 to appearing in the reference.
    """
    import inspect

    needs_example = [
        n
        for n in sh.__all__
        if callable(getattr(sh, n)) or inspect.isclass(getattr(sh, n))
    ]
    missing = sorted(
        n for n in needs_example if ">>>" not in (inspect.getdoc(getattr(sh, n)) or "")
    )
    assert not missing, (
        f"{len(missing)} callable(s) have no worked example: {missing}. "
        "Add an Examples: block; it is executed by --doctest-modules."
    )
```

- [ ] **Step 2: Run it to watch it fail**

Run: `python -m pytest tests/test_docs.py -q`
Expected: FAIL — `docs/reference/` does not exist yet, so all 45 names are missing.

- [ ] **Step 3: Add the docs tooling to the dependency sources**

Append to `requirements.in`:

```
mkdocs-material
mkdocstrings[python]
```

Create `docs/requirements.txt` for Read the Docs, which installs only what the docs build needs:

```
mkdocs-material
mkdocstrings[python]
-e .
```

- [ ] **Step 4: Write `mkdocs.yml`**

```yaml
site_name: Shedding Hub
site_description: Data and statistical models for biomarker shedding.
site_url: https://shedding-hub.readthedocs.io/
repo_url: https://github.com/shedding-hub/shedding-hub
edit_uri: blob/main/docs/

theme:
  name: material
  features:
    - navigation.sections
    - content.code.copy
    - search.suggest
  palette:
    - media: "(prefers-color-scheme: light)"
      scheme: default
      toggle: {icon: material/weather-night, name: Dark mode}
    - media: "(prefers-color-scheme: dark)"
      scheme: slate
      toggle: {icon: material/weather-sunny, name: Light mode}

# docs/ holds the project's internal design records. A default build would
# publish every one of them as a page.
exclude_docs: |
  superpowers/

plugins:
  - search
  - mkdocstrings:
      handlers:
        python:
          options:
            docstring_style: google
            show_source: false
            show_root_heading: true
            heading_level: 3
            members_order: source

markdown_extensions:
  - admonition
  - pymdownx.highlight
  - pymdownx.superfences
  - pymdownx.snippets:
      base_path: ["."]

nav:
  - Home: index.md
  - Getting started: getting-started.md
  - Tutorial: tutorial.md
  - Modeling methods: modeling-methods.md
  - Reference:
      - Datasets: reference/datasets.md
      - Statistics: reference/statistics.md
      - Visualisation: reference/visualisation.md
      - Fitting: reference/fitting.md
      - Catalog and ensembles: reference/catalog.md
      - Simulation: reference/simulation.md
      - Choosing a source: reference/selection.md
```

- [ ] **Step 5: Write the narrative pages**

`docs/index.md` — a short landing page. Include the development-version banner the design calls for:

```markdown
# Shedding Hub

Data and statistical models for biomarker shedding — viral RNA, drug
metabolites — in human specimens, for wastewater-based epidemiology.

!!! note "This documents the development version"
    These pages are built from `main`. The current release on PyPI is 0.2.0;
    behaviour described here may not be in it yet.

```bash
pip install shedding-hub
```

- [Getting started](getting-started.md) — load a dataset, summarise it, plot it.
- [Tutorial](tutorial.md) — simulate shedding for synthetic individuals.
- [Modeling methods](modeling-methods.md) — what the fitted estimates mean, and what they do not support.
```

`docs/getting-started.md` — install, then the three things a new user does, each a fenced `python` block: `load_dataset`, `calc_shedding_summary`, `plot_time_course`. Keep it under a screen; the reference carries the detail.

`docs/tutorial.md` — include the notebook rather than copying it, so there is one source:

```markdown
# Tutorial: simulating shedding

--8<-- "examples/simulating-shedding.md"
```

- [ ] **Step 6: Write the seven reference pages**

Each is a thin mkdocstrings stub. `docs/reference/selection.md`:

```markdown
# Choosing a source

::: shedding_hub.shedding_options
::: shedding_hub.shedding_for
::: shedding_hub.Selection
::: shedding_hub.classify_reference_event
::: shedding_hub.REFERENCE_EVENT_CLASSES
```

`docs/reference/visualisation.md` embeds the generated figure under each entry:

```markdown
# Visualisation

::: shedding_hub.plot_time_course

![plot_time_course](../images/plot_time_course.png)

::: shedding_hub.plot_time_courses

![plot_time_courses](../images/plot_time_courses.png)
```

…continuing for all 14 plotting functions. Distribute the remaining names across `datasets.md` (`load_dataset`, `check_dataset`, `normalize_str`, `folded_str`, `literal_str`), `statistics.md` (the six `calc_*`/`compare_datasets` names plus the duration and peak `calc_*` functions), `fitting.md` (`fit_shedding_model`, `SheddingFit`, `SheddingDataError`, `MODELS`, `PARAM_NAMES`), `catalog.md` (`load_shedding_catalog`, `SheddingCatalog`, `fit_shedding_models`, `make_ensemble`, `SheddingEnsemble`) and `simulation.md` (`simulate_shedding`, `plot_simulated_shedding`).

- [ ] **Step 7: Run the coverage test**

Run: `python -m pytest tests/test_docs.py -q`
Expected: PASS. If it names missing entries, add them to the matching page — that is the test doing its job.

- [ ] **Step 8: Build the site**

```bash
python scripts/build_doc_figures.py
python -m mkdocs build --strict
```
Expected: builds into `site/` with no warnings. `--strict` turns broken links and missing nav targets into failures.

- [ ] **Step 9: Check the tutorial rendered rather than stubbed**

```bash
python - <<'PY'
import pathlib
html = pathlib.Path("site/tutorial/index.html").read_text(encoding="utf-8")
assert "shedding_for" in html, "the notebook include did not render"
assert "jupytext" not in html.split("<body")[1][:2000], "front matter leaked into the page"
print("tutorial rendered, front matter stripped")
PY
```

If the jupytext front matter renders, strip it with a `pymdownx.snippets` section marker in `examples/simulating-shedding.md` and reference that section from `docs/tutorial.md`.

- [ ] **Step 10: Commit**

```bash
black tests/test_docs.py
git add mkdocs.yml requirements.in docs/ tests/test_docs.py
git commit -m "docs: a generated reference that a test keeps honest"
```

---

### Task 8: Publish on Read the Docs, and build the docs in CI

**Files:**
- Create: `.readthedocs.yaml`
- Modify: `.github/workflows/build.yaml`

**Interfaces:**
- Consumes: Task 7's `mkdocs.yml`.
- Produces: a published site, and a CI gate that fails on a broken docs build.

- [ ] **Step 1: Write `.readthedocs.yaml`**

```yaml
version: 2

build:
  os: ubuntu-24.04
  tools:
    python: "3.11"
  jobs:
    pre_build:
      # The reference embeds one figure per plotting function. They are
      # generated rather than committed so a stale chart cannot outlive the
      # code that drew it.
      - python scripts/build_doc_figures.py

mkdocs:
  configuration: mkdocs.yml
  fail_on_warning: true

python:
  install:
    - requirements: docs/requirements.txt
```

- [ ] **Step 2: Add a CI step so a broken docs build fails the PR**

In `.github/workflows/build.yaml`, after the `Assert the wheel ships only the package` step:

```yaml
      - name: Build the documentation
        run: |
          pip install mkdocs-material "mkdocstrings[python]"
          python scripts/build_doc_figures.py
          python -m mkdocs build --strict
```

- [ ] **Step 3: Verify the workflow parses**

```bash
python -c "import yaml; d=yaml.safe_load(open('.github/workflows/build.yaml')); print([s.get('name') for s in d['jobs']['build']['steps']])"
```
Expected: the list ends with `Build the documentation`.

- [ ] **Step 4: Run the full suite and the docs build locally**

```bash
python -m pytest -q
python scripts/build_doc_figures.py && python -m mkdocs build --strict
black --check .
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add .readthedocs.yaml .github/workflows/build.yaml
git commit -m "docs: publish on Read the Docs, and gate the docs build in CI"
```

- [ ] **Step 6: Connect the repository on Read the Docs**

This is a manual step in the Read the Docs web interface and cannot be scripted here: import the `shedding-hub/shedding-hub` project, confirm it detects `.readthedocs.yaml`, and trigger the first build. Report the resulting URL — expected to be `https://shedding-hub.readthedocs.io/`.

---

## Follow-up, not in this plan

`shedding-hub.github.io/package.html` still carries the hand-maintained function
list, including the three calls that raised `AttributeError` before Task 1. That
file lives in the **`shedding-hub/shedding-hub.github.io` repository**, which is
not checked out here, so it needs its own change: reduce the page to install
instructions, a short quickstart, and a link to the Read the Docs site, deleting
the function list rather than migrating it.

## Self-Review Notes

**Spec coverage.** Read the Docs hosting → Task 8. MkDocs Material +
mkdocstrings → Task 7. Three orphans exported with tests first → Task 1.
Latest-only versioning with a banner → Task 7 Step 5 (`docs/index.md`) and Task
8. Worked examples for every name → Tasks 3–5, with the six broken ones repaired
in Task 2. Doctests executed in CI → Task 2. Rendered figures → Task 6.
`exclude_docs` for the internal records → Task 7 Step 4. Single source for the
tutorial → Task 7 Step 5. The anti-drift test → Task 7 Step 1. `package.html` is
explicitly deferred above, since it is in another repository.

**On expected output.** Tasks 3–5 give the example *code* and four fully
verified outputs, but instruct the implementer to obtain the rest by running.
This is deliberate, not a placeholder: expected output cannot be known without
executing, and a plan that invented values would produce exactly the fabricated
documentation this whole effort exists to eliminate.

**Naming consistency.** `scripts/build_doc_figures.py`, `docs/images/<name>.png`,
`docs/reference/*.md` and `tests/test_docs.py` are referred to identically across
Tasks 6, 7 and 8.

**Known risk.** Task 6's figure lambdas guess some call signatures from the
existing README and tests. Step 2 says to correct any that differ rather than
wrap them in `try`, so a wrong guess surfaces as a build failure rather than a
missing image.
