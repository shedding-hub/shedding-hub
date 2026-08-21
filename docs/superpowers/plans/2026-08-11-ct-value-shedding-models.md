# Cycle-threshold shedding models — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `fit_shedding_model` and `plot_fit_diagnostic` work on cycle-threshold analytes — 27 studies, 69 analytes, 12,537 measurements currently refused outright — without changing anything about concentration fits or the published catalog.

**Architecture:** Ct is affine in log10 concentration (`Ct = α − β·log10 C`), so the existing gamma/exponential models already describe it. The response is transformed once, at the single point in `prepare_observations` where a reported value becomes a model value: `depth = CT_REFERENCE − Ct`, with `CT_REFERENCE = 40.0` fixed for every analyte. Everything downstream — the left-censored likelihood, the optimizer, the population summary — is untouched.

**Tech Stack:** Python 3.10+, numpy, pandas, scipy, matplotlib, pytest.

Design: `docs/superpowers/specs/2026-08-11-ct-value-shedding-models-design.md`

## Global Constraints

- Python floor is **3.10** (`requires-python = ">=3.10"` in `pyproject.toml`). PEP 604 unions in annotations are fine.
- CI runs `black --check .` — run `black .` before every commit.
- pytest runs with `--doctest-modules` over `tests` and `shedding_hub`, so **every docstring example is executed as a test**.
- Doctests must not print numpy scalars directly — local and CI numpy versions repr them differently. Wrap in `int()` / `float()`.
- `CT_REFERENCE` is exactly `40.0`. It is a published convention: once a fit is serialized carrying it, changing it silently reinterprets stored heights.
- **New `SheddingFit` fields must have defaults.** Catalogs serialized before a field existed must stay loadable — precedent is `n_degenerate_subjects`.
- The published catalog (`shedding_hub/data/shedding_catalog.yaml`) must not gain or lose fits in this work, and no fitted number in it may change. It currently records **207** skips with `reason: ct_units` and that count must stay 207. The skip *message* text does change (Task 4) — that is the only permitted difference, and Task 12 verifies it is the only one.

---

### Task 1: `CT_REFERENCE` and the response transform

**Files:**
- Modify: `shedding_hub/shedding_fit.py` (constants near `NEGATIVE_VALUE`, line ~22)
- Test: `tests/test_shedding_fit.py`

**Interfaces:**
- Consumes: nothing
- Produces: `CT_REFERENCE: float`, `_to_response(value: float, value_type: str) -> float`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_shedding_fit.py`, and add `CT_REFERENCE, _to_response` to the existing `from shedding_hub.shedding_fit import (...)` block at the top:

```python
def test_to_response_maps_concentration_to_log10():
    assert _to_response(1e6, "concentration") == 6.0


def test_to_response_maps_ct_to_cycles_below_the_reference():
    # Ct 31 is the repository median; 40 - 31 = 9 cycles below the reference.
    assert _to_response(31.0, "ct") == 9.0


def test_to_response_is_decreasing_in_ct():
    # The sign flip is the whole point: less virus means a HIGHER Ct, so a
    # higher Ct must map to a LOWER response. Drop the negation and every
    # fitted curve inverts while still converging happily.
    assert _to_response(20.0, "ct") > _to_response(35.0, "ct")


def test_ct_reference_is_forty():
    assert CT_REFERENCE == 40.0
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/test_shedding_fit.py -k to_response -v`
Expected: FAIL at import — `ImportError: cannot import name '_to_response'`

- [ ] **Step 3: Implement**

In `shedding_hub/shedding_fit.py`, after `NEGATIVE_VALUE = "negative"`:

```python
# Ct is affine in log10 concentration (Ct = alpha - beta * log10 C), so the
# shedding models describe it unchanged once it is negated -- Ct falls as
# shedding rises -- and offset to keep fitted levels positive.
#
# The offset is one constant for every analyte, not each study's own cutoff.
# Recorded cutoffs run 37 to 41, so anchoring per study would put two studies
# measuring identical samples up to 4 cycles apart on height alone. 40 sits
# above the observed Ct median of 31 and above 95% of all readings, so fitted
# peak heights -- which occur at LOW Ct -- stay comfortably positive.
CT_REFERENCE = 40.0


def _to_response(value: float, value_type: str) -> float:
    """
    Map a reported measurement onto the scale the models are fitted on.

    Concentrations are fitted on log10. Cycle thresholds are fitted as cycles
    below ``CT_REFERENCE``, which is decreasing in Ct and therefore increasing
    in viral load, exactly like a log10 concentration.
    """
    if value_type == "ct":
        return CT_REFERENCE - float(value)
    return math.log10(float(value))
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/test_shedding_fit.py -k "to_response or ct_reference" -v`
Expected: PASS, 4 tests

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git add shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git commit -m "feat: add the Ct response transform"
```

---

### Task 2: Resolve the censoring limit on the Ct scale

**Files:**
- Modify: `shedding_hub/shedding_fit.py:143-178` (`_resolve_censoring_limit`)
- Test: `tests/test_shedding_fit.py`

**Interfaces:**
- Consumes: `_to_response` from Task 1
- Produces: `_resolve_censoring_limit(analyte_spec: dict, observed: np.ndarray, value_type: str = "concentration") -> float`

- [ ] **Step 1: Write the failing tests**

Add `_resolve_censoring_limit` to the import block, then:

```python
def test_censoring_limit_on_the_ct_scale():
    # An assay running to Ct 41 detects one cycle PAST the reference, so the
    # limit is negative. Nothing in the likelihood cares.
    spec = {"limit_of_detection": 41, "limit_of_quantification": "unknown"}
    assert _resolve_censoring_limit(spec, np.array([9.0, 5.0]), "ct") == -1.0


def test_censoring_limit_on_the_ct_scale_for_a_stricter_assay():
    spec = {"limit_of_detection": 37, "limit_of_quantification": "unknown"}
    assert _resolve_censoring_limit(spec, np.array([9.0]), "ct") == 3.0


def test_censoring_limit_for_concentration_is_unchanged():
    spec = {"limit_of_quantification": 100, "limit_of_detection": "unknown"}
    assert _resolve_censoring_limit(spec, np.array([6.0]), "concentration") == 2.0


def test_censoring_limit_falls_back_below_the_smallest_observed_ct():
    # No declared limit: fall back below the smallest response, which on the Ct
    # scale means just past the HIGHEST observed Ct. The fallback arithmetic is
    # identical in both spaces because `observed` is already transformed.
    spec = {"limit_of_detection": "unknown", "limit_of_quantification": "unknown"}
    with pytest.warns(UserWarning, match="cycles below"):
        limit = _resolve_censoring_limit(spec, np.array([9.0, 2.0]), "ct")
    assert limit < 2.0
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/test_shedding_fit.py -k censoring_limit -v`
Expected: FAIL — `_resolve_censoring_limit() takes 2 positional arguments but 3 were given`

- [ ] **Step 3: Implement**

Replace the body of `_resolve_censoring_limit`. Keep the existing docstring and append the paragraph shown:

```python
def _resolve_censoring_limit(
    analyte_spec: dict,
    observed: np.ndarray,
    value_type: str = "concentration",
) -> float:
    """
    ... keep the existing docstring text unchanged, then append: ...

    For a cycle-threshold analyte the declared limit is itself a Ct, so it is
    transformed the same way the observations are and the resolved limit is
    ``CT_REFERENCE - cutoff`` — zero at a cutoff of 40, negative where an assay
    runs to 41, positive where it stops at 37. The fallback branch needs no
    special case: ``observed`` has already been transformed, so "just below the
    smallest response" means the same thing on either scale.
    """
    for key in ("limit_of_quantification", "limit_of_detection"):
        limit = _numeric_limit(analyte_spec.get(key))
        if limit is not None:
            return _to_response(limit, value_type)

    smallest = float(observed.min())
    fallback = smallest - CENSORING_MARGIN
    scale = "cycles below reference" if value_type == "ct" else "log10"
    warnings.warn(
        "Falling back to a censoring limit of "
        f"{fallback:.4g} ({scale}) because no limit of quantification or "
        "detection is declared for this analyte.",
        UserWarning,
        stacklevel=2,
    )
    return fallback
```

Rename the parameter from `observed_log10` to `observed` in the signature only; update the one call site at `shedding_hub/shedding_fit.py:455` in Task 5.

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_shedding_fit.py -k censoring_limit -v`
Expected: PASS, 4 tests

Run: `pytest tests/test_shedding_fit.py -v`
Expected: PASS — the default argument keeps every existing caller working

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git add shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git commit -m "feat: resolve censoring limits on the Ct scale"
```

---

### Task 3: `Observations` carries its value type

**Files:**
- Modify: `shedding_hub/shedding_fit.py:106-126` (`Observations`), `:181-189` (`_record_dropped`)
- Test: `tests/test_shedding_fit.py`

**Interfaces:**
- Consumes: `_to_response` from Task 1
- Produces: `Observations.value_type: str` (default `"concentration"`); `_record_dropped(measurement, time, times, values, value_type="concentration")`

- [ ] **Step 1: Write the failing test**

```python
def test_dropped_ct_readings_are_recorded_on_the_response_scale():
    # Dropped points are drawn on the diagnostic plot, so they must share the
    # scale of the points that were kept or they land in the wrong place.
    times: list[float] = []
    values: list[float] = []
    _record_dropped({"value": 28.0}, -1.0, times, values, "ct")
    assert values == [12.0]


def test_observations_default_to_concentration():
    assert Observations(
        subject_index=np.zeros(1, int),
        times=np.zeros(1),
        values=np.zeros(1),
        censored=np.zeros(1, bool),
        censoring_limit=0.0,
    ).value_type == "concentration"
```

Add `Observations` and `_record_dropped` to the import block.

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_shedding_fit.py -k "dropped_ct or observations_default" -v`
Expected: FAIL — `_record_dropped() takes 4 positional arguments but 5 were given`

- [ ] **Step 3: Implement**

Add the field to `Observations`, after `censoring_limit`:

```python
    censoring_limit: float
    # Which scale ``values`` and ``censoring_limit`` are on: log10 concentration,
    # or cycles below CT_REFERENCE. Carried on the observations rather than
    # re-derived downstream, so a plot or a fit cannot disagree with the fitter
    # about what its own numbers mean.
    value_type: str = "concentration"
```

Then `_record_dropped`:

```python
def _record_dropped(
    measurement: dict,
    time: float,
    times: list,
    values: list,
    value_type: str = "concentration",
) -> None:
    """Note a discarded reading, if it can be placed on a plot at all."""
    value = measurement.get("value")
    if isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0:
        times.append(time)
        values.append(_to_response(float(value), value_type))
    elif value == NEGATIVE_VALUE:
        times.append(time)
        values.append(float("nan"))
```

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_shedding_fit.py -v`
Expected: PASS — `value_type` has a default, so existing construction sites are unaffected

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git add shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git commit -m "feat: record the value type on Observations"
```

---

### Task 4: Keep Ct analytes out of the catalog — **before** the raise is removed

This task must land before Task 5. `fit_shedding_models` currently relies on `prepare_observations` raising to skip Ct analytes; remove that raise first and the published catalog silently gains up to 207 fits whose heights are not commensurable with the rest.

**Files:**
- Modify: `shedding_hub/shedding_catalog.py:273-330` (`fit_shedding_models` signature and analyte loop)
- Test: `tests/test_shedding_catalog.py`

**Interfaces:**
- Consumes: `_is_ct_unit` from `shedding_hub/shedding_fit.py:129`
- Produces: `fit_shedding_models(..., value_types: tuple[str, ...] = ("concentration",))`

- [ ] **Step 1: Write the failing test**

In `tests/test_shedding_catalog.py`:

```python
def test_catalog_skips_ct_analytes_by_default(ct_dataset):
    catalog = fit_shedding_models([ct_dataset], models=("gamma",))
    assert len(catalog.fits) == 0
    assert (catalog.skipped["reason"] == "ct_units").all()


def test_catalog_fits_ct_analytes_when_asked(ct_dataset):
    catalog = fit_shedding_models(
        [ct_dataset], models=("gamma",), value_types=("concentration", "ct")
    )
    assert len(catalog.fits) == 1
    assert catalog.fits[0].value_type == "ct"
```

`test_catalog_fits_ct_analytes_when_asked` depends on Tasks 5 and 6; mark it `@pytest.mark.xfail(reason="enabled by Task 6", strict=False)` here and remove the marker in Task 6.

Add the `ct_dataset` fixture to `tests/conftest.py` so both test modules share it:

```python
@pytest.fixture
def ct_dataset():
    """One analyte in cycle threshold, three subjects that rise then fall."""
    curve = [(1, 32.0), (3, 24.0), (5, 22.0), (8, 27.0), (12, 33.0), (16, "negative")]
    return {
        "dataset_id": "ct_study",
        "analytes": {
            "swab": {
                "specimen": "nasopharyngeal_swab",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "cycle threshold",
                "limit_of_detection": 40,
                "limit_of_quantification": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "swab", "time": t, "value": v}
                    for t, v in [(t, v if isinstance(v, str) else v + shift)
                                 for t, v in curve]
                ]
            }
            for shift in (0.0, 1.5, -1.0)
        ],
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_shedding_catalog.py -k ct_analytes -v`
Expected: `test_catalog_skips_ct_analytes_by_default` FAILS with `TypeError: fit_shedding_models() got an unexpected keyword argument` only after Step 3 adds the parameter; before that it fails because the fixture is unknown. Either failure is the expected red.

- [ ] **Step 3: Implement**

Add the parameter to `fit_shedding_models` (alongside `min_observations`, `min_time`, `max_peak_above_observed` at `shedding_hub/shedding_catalog.py:277-279`):

```python
    value_types: tuple[str, ...] = ("concentration",),
```

Document it in the docstring's Args block:

```
        value_types: Which measurement scales may enter the catalog. Defaults
            to concentration only. Cycle-threshold fits are individually valid
            -- ``fit_shedding_model`` produces them -- but their heights are
            cycles below ``CT_REFERENCE`` rather than log10 concentrations, so
            an ensemble that averaged the two would be averaging incommensurable
            quantities. Opt in with ``("concentration", "ct")`` once that is
            resolved.
```

Then, inside `for dataset in datasets:` replace `for analyte in dataset.get("analytes", {}):` with:

```python
        for analyte, analyte_spec in dataset.get("analytes", {}).items():
            # Skipped here rather than left to prepare_observations, which now
            # accepts Ct analytes. Keeping the decision in the catalog builder is
            # what lets a caller fit one directly while the published catalog
            # stays concentration-only.
            if _is_ct_unit(analyte_spec.get("unit")) and "ct" not in value_types:
                for model in models:
                    skipped.append(
                        {
                            "dataset_id": dataset_id,
                            "analyte": analyte,
                            "model": model,
                            "reason": "ct_units",
                            "message": (
                                f"Analyte {analyte!r} is reported in "
                                f"{analyte_spec.get('unit')!r}. Cycle-threshold "
                                "fits are supported but excluded from the "
                                "catalog, whose heights are log10 "
                                "concentrations."
                            ),
                        }
                    )
                continue
            for model in models:
```

Add `_is_ct_unit` to the existing `from .shedding_fit import (...)` block at `shedding_hub/shedding_catalog.py:19`.

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_shedding_catalog.py -v`
Expected: PASS, with `test_catalog_fits_ct_analytes_when_asked` xfailing

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/shedding_catalog.py tests/test_shedding_catalog.py tests/conftest.py
git add shedding_hub/shedding_catalog.py tests/test_shedding_catalog.py tests/conftest.py
git commit -m "feat: gate catalog membership on value type"
```

---

### Task 5: `prepare_observations` accepts Ct analytes

**Files:**
- Modify: `shedding_hub/shedding_fit.py:257-264` (delete the raise), `:300`, `:327`, `:455`, `:457-465`, and the docstrings at `:68` and `:229-233`
- Test: `tests/test_shedding_fit.py` (including replacing the test at `:549`)

**Interfaces:**
- Consumes: `_to_response`, `_resolve_censoring_limit`, `Observations.value_type`, `_record_dropped` from Tasks 1–3
- Produces: `prepare_observations` returning `Observations` with `value_type="ct"` for cycle-threshold analytes

- [ ] **Step 1: Write the failing tests**

Replace the existing test at `tests/test_shedding_fit.py:549` that asserts `excinfo.value.reason == "ct_units"` — delete it, and add:

```python
def test_prepare_observations_accepts_ct_analytes(ct_dataset):
    obs = prepare_observations(ct_dataset, "swab", "gamma")
    assert obs.value_type == "ct"
    assert obs.n_subjects == 3


def test_ct_values_are_cycles_below_the_reference(ct_dataset):
    obs = prepare_observations(ct_dataset, "swab", "gamma")
    # First subject, day 1, Ct 32.0 -> 40 - 32 = 8.0.
    assert obs.values[0] == pytest.approx(8.0)


def test_ct_response_peaks_where_ct_is_lowest(ct_dataset):
    # The sign-flip guard at the level that matters. The fixture's lowest Ct is
    # 22.0 at day 5; that must be the LARGEST response, not the smallest.
    obs = prepare_observations(ct_dataset, "swab", "gamma")
    first = obs.subject_index == 0
    assert obs.times[first][np.argmax(obs.values[first])] == 5


def test_ct_censoring_limit_comes_from_the_declared_cutoff(ct_dataset):
    obs = prepare_observations(ct_dataset, "swab", "gamma")
    # Cutoff 40 == CT_REFERENCE, so the limit is exactly zero.
    assert obs.censoring_limit == pytest.approx(0.0)


def test_concentration_analytes_are_untouched(simple_dataset):
    obs = prepare_observations(simple_dataset, "stool", "gamma")
    assert obs.value_type == "concentration"
    assert obs.values[0] == pytest.approx(6.0)
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/test_shedding_fit.py -k ct_ -v`
Expected: FAIL — `SheddingDataError: ... Cycle-threshold values are inversely related ...` with `reason == "ct_units"`

- [ ] **Step 3: Implement**

Delete the raise at `shedding_hub/shedding_fit.py:257-264` entirely and put the value type in its place:

```python
    # Cycle thresholds are affine in log10 concentration, so both models
    # describe them once the response is transformed. See ``_to_response``.
    value_type = "ct" if _is_ct_unit(analyte_spec.get("unit")) else "concentration"
```

At `:300`, pass the value type when recording a dropped reading:

```python
                _record_dropped(
                    measurement, time, dropped_times, dropped_values, value_type
                )
```

At `:327`, transform the kept reading:

```python
            values.append(_to_response(float(value), value_type))
```

At `:455`:

```python
    censoring_limit = _resolve_censoring_limit(analyte_spec, observed, value_type)
```

And in the `return Observations(...)` block add:

```python
        value_type=value_type,
```

Update the two docstrings that advertise the removed reason. At `:68`, drop `ct_units` from the list of reasons. At `:229-233`, drop `"uses cycle-threshold units"` from the `Raises:` text and remove `ct_units` from the reason list, then add to the `Returns:` text:

```
        An ``Observations`` instance with subject indices renumbered
        contiguously. For a cycle-threshold analyte, ``values`` are cycles
        below ``CT_REFERENCE`` rather than log10 concentrations, and
        ``value_type`` says which.
```

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_shedding_fit.py -v`
Expected: PASS

Run: `pytest tests/ shedding_hub/ -q`
Expected: PASS — Task 4's catalog gate keeps `tests/test_shedding_catalog.py` green

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git add shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git commit -m "feat: fit cycle-threshold analytes"
```

---

### Task 6: `SheddingFit` records its scale and what may be compared

**Files:**
- Modify: `shedding_hub/shedding_fit.py:705-760` (dataclass fields), `:1589+` (construction), and add the comparability constant near `CT_REFERENCE`
- Test: `tests/test_shedding_fit.py`

**Interfaces:**
- Consumes: `Observations.value_type` from Task 3
- Produces: `SheddingFit.value_type`, `.ct_reference`, `.ct_cutoff`, `.comparable_with(other) -> tuple[str, ...]`; module constant `VALUE_TYPE_INVARIANT_PARAMETERS`

- [ ] **Step 1: Write the failing tests**

```python
def test_ct_fit_records_its_scale(ct_dataset):
    fit = fit_shedding_model(ct_dataset, analyte="swab", model="gamma")
    assert fit.value_type == "ct"
    assert fit.ct_reference == 40.0
    assert fit.ct_cutoff == 40.0


def test_concentration_fit_has_no_ct_metadata(simple_dataset):
    fit = fit_shedding_model(simple_dataset, analyte="stool", model="gamma")
    assert fit.value_type == "concentration"
    assert fit.ct_reference is None


def test_only_temporal_parameters_compare_across_value_types(
    ct_dataset, simple_dataset
):
    ct = fit_shedding_model(ct_dataset, analyte="swab", model="gamma")
    conc = fit_shedding_model(simple_dataset, analyte="stool", model="gamma")
    assert ct.comparable_with(conc) == VALUE_TYPE_INVARIANT_PARAMETERS
    assert "peak_day" in ct.comparable_with(conc)
    assert "half_life_days" not in ct.comparable_with(conc)


def test_everything_compares_within_a_value_type(ct_dataset):
    fit = fit_shedding_model(ct_dataset, analyte="swab", model="gamma")
    assert "half_life_days" in fit.comparable_with(fit)
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/test_shedding_fit.py -k "records_its_scale or compare" -v`
Expected: FAIL — `AttributeError: 'SheddingFit' object has no attribute 'value_type'`

- [ ] **Step 3: Implement**

Near `CT_REFERENCE`:

```python
# Peak time, onset and rise duration are ratios of b0 to a0. Fitting Ct returns
# a0 and b0 both multiplied by the unknown standard-curve slope, so the slope
# cancels in the ratio and these three transfer between value types exactly,
# with no assumption about PCR efficiency. Everything else carries either the
# slope (a0, half-life) or the assay intercept (height), and studies report
# neither.
VALUE_TYPE_INVARIANT_PARAMETERS = ("peak_day", "t0", "rise_days")

ALL_COMPARABLE_PARAMETERS = VALUE_TYPE_INVARIANT_PARAMETERS + (
    "a0",
    "half_life_days",
    "peak_height",
)
```

Add to `SheddingFit`, at the end of the field list so defaults stay valid:

```python
    # Which scale this fit's heights live on. Defaults to concentration so a
    # catalog serialized before Ct support stays loadable and reads correctly.
    value_type: str = "concentration"
    # The reference heights are measured below, and the analyte's own detection
    # cutoff. Recorded rather than assumed: a fit serialized under one
    # convention must not be silently reinterpreted under another, and a height
    # reads back as a minimum Ct via ``ct_reference - height``. Both None for
    # concentration fits.
    ct_reference: float | None = None
    ct_cutoff: float | None = None
```

Add the method to `SheddingFit`:

```python
    def comparable_with(self, other: "SheddingFit") -> tuple[str, ...]:
        """
        Parameters of this fit that may be compared with ``other``'s.

        Within one value type everything compares. Across value types only the
        temporal parameters do — see ``VALUE_TYPE_INVARIANT_PARAMETERS``.

        Examples:
            >>> import shedding_hub as sh
            >>> data = sh.load_dataset('woelfel2020virological', local='./data')
            >>> fit = sh.fit_shedding_model(data, analyte='stool', model='gamma')
            >>> 'half_life_days' in fit.comparable_with(fit)
            True
        """
        if self.value_type == other.value_type:
            return ALL_COMPARABLE_PARAMETERS
        return VALUE_TYPE_INVARIANT_PARAMETERS
```

At the `return SheddingFit(` construction (`shedding_hub/shedding_fit.py:1589`), add:

```python
        value_type=observations.value_type,
        ct_reference=CT_REFERENCE if observations.value_type == "ct" else None,
        ct_cutoff=(
            _numeric_limit(analyte_spec.get("limit_of_detection"))
            or _numeric_limit(analyte_spec.get("limit_of_quantification"))
            if observations.value_type == "ct"
            else None
        ),
```

Export the new names from `shedding_hub/__init__.py`: add `CT_REFERENCE` and `VALUE_TYPE_INVARIANT_PARAMETERS` to the `from .shedding_fit import (...)` block and to `__all__`.

Finally, remove the `@pytest.mark.xfail` marker added in Task 4 from `test_catalog_fits_ct_analytes_when_asked`.

- [ ] **Step 4: Run the tests**

Run: `pytest tests/ shedding_hub/ -q`
Expected: PASS, including the previously-xfailing catalog test

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/ tests/
git add shedding_hub/ tests/
git commit -m "feat: record value type and comparability on SheddingFit"
```

---

### Task 7: Name the population height coordinate for its scale

**Files:**
- Modify: `shedding_hub/shedding_models.py:43-51` (add the helper after `POPULATION_COORDS`)
- Test: `tests/test_shedding_models.py`

**Interfaces:**
- Consumes: `POPULATION_COORDS`
- Produces: `population_coord_names(model: str, value_type: str = "concentration") -> tuple[str, ...]`

- [ ] **Step 1: Write the failing tests**

```python
def test_population_coord_names_default_to_concentration():
    assert population_coord_names("gamma") == POPULATION_COORDS["gamma"]


def test_population_coord_names_rename_the_height_for_ct():
    assert population_coord_names("gamma", "ct") == (
        "log_a0",
        "log_peak_day",
        "peak_cycles",
    )


def test_population_coord_names_leave_temporal_coordinates_alone():
    # t0 is a time on either scale and must not be renamed.
    assert population_coord_names("gamma_shifted", "ct")[-1] == "t0"
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_shedding_models.py -k population_coord_names -v`
Expected: FAIL — `ImportError: cannot import name 'population_coord_names'`

- [ ] **Step 3: Implement**

```python
def population_coord_names(
    model: str, value_type: str = "concentration"
) -> tuple[str, ...]:
    """
    Coordinate names for a model's population summary on a given value scale.

    Only the height coordinate differs. On the Ct scale it is cycles below
    ``CT_REFERENCE``, not a log10 concentration, and naming it ``peak_log10``
    there would invite exactly the cross-scale comparison that is invalid. The
    temporal coordinates keep their names because they mean the same thing on
    either scale.

    Examples:
        >>> from shedding_hub.shedding_models import population_coord_names
        >>> population_coord_names('gamma', 'ct')
        ('log_a0', 'log_peak_day', 'peak_cycles')
    """
    validate_model(model)
    names = POPULATION_COORDS[model]
    if value_type != "ct":
        return names
    return tuple("peak_cycles" if name == "peak_log10" else name for name in names)
```

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_shedding_models.py shedding_hub/shedding_models.py -v`
Expected: PASS

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/shedding_models.py tests/test_shedding_models.py
git add shedding_hub/shedding_models.py tests/test_shedding_models.py
git commit -m "feat: name the population height coordinate for its scale"
```

---

### Task 8: Parameter recovery on synthetic Ct data

Proves the pipeline recovers what it should, and pins the peak-time invariance claim the whole design rests on.

**Files:**
- Create: `tests/test_shedding_fit_ct.py`

**Interfaces:**
- Consumes: everything from Tasks 1–6

- [ ] **Step 1: Write the failing test**

```python
"""Recovery and invariance checks for cycle-threshold fitting."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from shedding_hub import fit_shedding_model
from shedding_hub.shedding_fit import CT_REFERENCE
from shedding_hub.shedding_models import log10_concentration


def _synthetic(value_type, a0=0.25, b0=1.5, c0=12.0, n_subjects=25, seed=0):
    """
    One cohort, emitted either as concentrations or as the Ct values that the
    SAME curve would have produced through an assay with slope 3.5 and
    intercept 38. Peak time is b0 / a0 = 6.0 days for both.
    """
    rng = np.random.default_rng(seed)
    times = np.array([1.0, 2.0, 4.0, 6.0, 8.0, 11.0, 15.0, 20.0])
    participants = []
    for _ in range(n_subjects):
        params = np.array([[a0, b0, c0 + rng.normal(0.0, 0.8)]])
        y = log10_concentration("gamma", params, times)[0]
        y = y + rng.normal(0.0, 0.2, y.size)
        if value_type == "ct":
            values = 38.0 - 3.5 * y
        else:
            values = 10.0**y
        participants.append(
            {
                "measurements": [
                    {"analyte": "a", "time": float(t), "value": float(v)}
                    for t, v in zip(times, values)
                ]
            }
        )
    unit = "cycle threshold" if value_type == "ct" else "gc/mL"
    return {
        "dataset_id": f"synthetic_{value_type}",
        "analytes": {
            "a": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": unit,
                "limit_of_detection": 40 if value_type == "ct" else 10,
                "limit_of_quantification": "unknown",
            }
        },
        "participants": participants,
    }


def test_ct_fit_recovers_the_true_peak_time():
    fit = fit_shedding_model(_synthetic("ct"), analyte="a", model="gamma")
    assert fit.peak_day == pytest.approx(6.0, rel=0.25)


def test_ct_and_concentration_agree_on_peak_time():
    # The central claim: the standard-curve slope cancels in b0/a0, so both
    # scales must land on the same peak day even though the assay applied a
    # slope of 3.5 and an intercept of 38 to one of them.
    ct = fit_shedding_model(_synthetic("ct"), analyte="a", model="gamma")
    conc = fit_shedding_model(_synthetic("concentration"), analyte="a", model="gamma")
    assert ct.peak_day == pytest.approx(conc.peak_day, rel=0.15)


def test_ct_decay_rate_carries_the_assay_slope():
    # The other half of the claim, stated as a fact rather than a hope: a0 is
    # NOT invariant. It comes back scaled by the assay slope of 3.5, which is
    # exactly why half_life_days is excluded from comparable_with.
    ct = fit_shedding_model(_synthetic("ct"), analyte="a", model="gamma")
    conc = fit_shedding_model(_synthetic("concentration"), analyte="a", model="gamma")
    a0_ct = float(ct.subject_params["a0"].median())
    a0_conc = float(conc.subject_params["a0"].median())
    assert a0_ct / a0_conc == pytest.approx(3.5, rel=0.35)


def test_a_ct_fit_is_a_peak_not_a_trough():
    # If the sign flip were dropped anywhere in the chain, the optimizer would
    # still converge -- on an inverted curve. Assert the fitted median rises
    # from day 1 to the peak and falls after it.
    fit = fit_shedding_model(_synthetic("ct"), analyte="a", model="gamma")
    params = fit.median_params()
    heights = log10_concentration("gamma", params, np.array([1.0, 6.0, 20.0]))[0]
    assert heights[1] > heights[0]
    assert heights[1] > heights[2]
```

`SheddingFit.median_params()` is defined at `shedding_hub/shedding_fit.py:760` and returns a 1-D array of length `k`. `log10_concentration` calls `np.atleast_2d` on its `params`, so passing that array straight through gives a `(1, m)` result and `[0]` indexes the single row — as written above.

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_shedding_fit_ct.py -v`
Expected: FAIL if any earlier task is incomplete; this is the integration gate

- [ ] **Step 3: No implementation**

This task adds no production code. If a test fails, the defect is in Tasks 1–6 — fix it there rather than loosening the tolerance.

- [ ] **Step 4: Run the whole suite**

Run: `pytest tests/ shedding_hub/ -q`
Expected: PASS

- [ ] **Step 5: Format and commit**

```bash
black tests/test_shedding_fit_ct.py
git add tests/test_shedding_fit_ct.py
git commit -m "test: recover Ct parameters and pin peak-time invariance"
```

---

### Task 9: Diagnostic plot on the Ct scale

**Note — this corrects the spec.** The design says the plot should invert the Ct axis. It should not: the fitted response is already `40 − Ct`, which increases with viral load, so plotting it puts low Ct at the top with no inversion at all. Inverting on top of that would flip it back upside down. What the axis needs is a **tick formatter** showing real Ct numbers, and a label saying so.

**Files:**
- Modify: `shedding_hub/viz.py:3429-3430` (axis labels in `plot_fit_diagnostic`)
- Test: `tests/test_viz.py`

**Interfaces:**
- Consumes: `SheddingFit.value_type`, `CT_REFERENCE`

- [ ] **Step 1: Write the failing tests**

```python
def test_fit_diagnostic_labels_the_ct_axis(ct_dataset):
    fit = fit_shedding_model(ct_dataset, analyte="swab", model="gamma")
    fig = plot_fit_diagnostic(fit, ct_dataset)
    assert fig.axes[0].get_ylabel() == "Ct (cycle threshold)"


def test_fit_diagnostic_ct_ticks_read_as_ct_not_as_depth():
    # A response of 8.0 is Ct 32. The tick must say 32.
    fit = fit_shedding_model(ct_dataset(), analyte="swab", model="gamma")
    fig = plot_fit_diagnostic(fit, ct_dataset())
    formatter = fig.axes[0].yaxis.get_major_formatter()
    assert formatter(8.0, None) == "32"


def test_fit_diagnostic_is_not_inverted_for_ct(ct_dataset):
    # The response already increases with viral load. Inverting would put the
    # peak at the bottom.
    fit = fit_shedding_model(ct_dataset, analyte="swab", model="gamma")
    fig = plot_fit_diagnostic(fit, ct_dataset)
    bottom, top = fig.axes[0].get_ylim()
    assert bottom < top


def test_fit_diagnostic_concentration_label_unchanged(woelfel_dataset):
    fit = fit_shedding_model(woelfel_dataset, analyte="stool", model="gamma")
    fig = plot_fit_diagnostic(fit, woelfel_dataset)
    assert "log10 concentration" in fig.axes[0].get_ylabel()
```

Use whatever concentration fixture `tests/test_viz.py` already has in place of `woelfel_dataset`.

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_viz.py -k ct -v`
Expected: FAIL — ylabel is `log10 concentration (cycle threshold)`

- [ ] **Step 3: Implement**

Replace `shedding_hub/viz.py:3429-3430`:

```python
    ax.set_xlabel(f"Days after {fit.reference_event}")
    if getattr(fit, "value_type", "concentration") == "ct":
        # Points are plotted as cycles below CT_REFERENCE, which rises with
        # viral load, so the curve already reads as a peak and the axis must
        # NOT be inverted. Only the tick labels need converting back to the Ct
        # the study actually reported.
        ax.set_ylabel("Ct (cycle threshold)")
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda value, _: f"{CT_REFERENCE - value:.0f}")
        )
    else:
        ax.set_ylabel(
            f"log10 concentration ({unit})" if unit else "log10 concentration"
        )
```

Two import changes at the top of `shedding_hub/viz.py`. Add a new line after `from matplotlib.figure import Figure` (line 6):

```python
from matplotlib.ticker import FuncFormatter
```

and extend the existing import at line 10 — do not add a second statement:

```python
from .shedding_fit import CT_REFERENCE, prepare_observations
```

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_viz.py -q`
Expected: PASS

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/viz.py tests/test_viz.py
git add shedding_hub/viz.py tests/test_viz.py
git commit -m "feat: plot Ct fits with Ct tick labels"
```

---

### Task 10: Ensembles refuse to mix value types

**Files:**
- Modify: `shedding_hub/shedding_ensemble.py` (`make_ensemble`)
- Test: `tests/test_shedding_ensemble.py`

**Interfaces:**
- Consumes: `SheddingFit.value_type`

- [ ] **Step 1: Write the failing test**

```python
def test_ensemble_refuses_to_mix_value_types(ct_dataset, simple_dataset):
    ct = fit_shedding_model(ct_dataset, analyte="swab", model="gamma")
    conc = fit_shedding_model(simple_dataset, analyte="stool", model="gamma")
    with pytest.raises(ValueError, match="value type"):
        make_ensemble([ct, conc])


def test_ensemble_accepts_ct_fits_on_their_own(ct_dataset):
    ct = fit_shedding_model(ct_dataset, analyte="swab", model="gamma")
    assert make_ensemble([ct, ct]) is not None
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_shedding_ensemble.py -k value_types -v`
Expected: FAIL — `make_ensemble` averages them without complaint

- [ ] **Step 3: Implement**

Near the top of `make_ensemble`, after its existing validation:

```python
    value_types = {getattr(fit, "value_type", "concentration") for fit in fits}
    if len(value_types) > 1:
        raise ValueError(
            "Cannot build an ensemble from fits of more than one value type "
            f"({sorted(value_types)}). A Ct fit's height is cycles below "
            "CT_REFERENCE and a concentration fit's is log10 concentration; "
            "averaging them would average incommensurable quantities. Peak "
            "times are comparable across value types -- see "
            "SheddingFit.comparable_with -- but the population summary is not."
        )
```

- [ ] **Step 4: Run the tests**

Run: `pytest tests/ shedding_hub/ -q`
Expected: PASS

- [ ] **Step 5: Format and commit**

```bash
black shedding_hub/shedding_ensemble.py tests/test_shedding_ensemble.py
git add shedding_hub/shedding_ensemble.py tests/test_shedding_ensemble.py
git commit -m "feat: refuse ensembles that mix value types"
```

---

### Task 11: Validate against a study reporting both scales

`kissler2021viral` measures Ct and concentration on the same 68 subjects at the same 2,406 timepoints. This turns the invariance claim from a simulation result into an empirical one.

**Files:**
- Modify: `tests/test_shedding_fit_ct.py`

- [ ] **Step 1: Write the test**

```python
def test_kissler_peak_times_agree_across_value_types():
    """
    The same 68 subjects, measured both ways at the same 2,406 timepoints.
    Peak time is a ratio of b0 to a0, so the assay's standard-curve slope
    cancels and the two fits must land on the same day.
    """
    import shedding_hub as sh

    data = sh.load_dataset("kissler2021viral", local="./data")
    ct = sh.fit_shedding_model(
        data, analyte="AN_OPS_SARSCoV2_ct", model="gamma"
    )
    conc = sh.fit_shedding_model(
        data, analyte="AN_OPS_SARSCoV2_viral", model="gamma"
    )
    # The concentration fit gives peak_day = 0.77 over 51 retained subjects
    # (measured 2026-08-11). A purely relative tolerance on a number that small
    # would be brittle, so allow a day either way: "the two scales agree on when
    # shedding peaks, to within a day" is the claim worth defending.
    assert ct.peak_day == pytest.approx(conc.peak_day, rel=0.25, abs=1.0)


def test_kissler_ct_and_concentration_are_affinely_related():
    """
    Ct = alpha - beta * log10 C is the assumption the whole design rests on.
    Check it directly on the matched pairs, and record the empirical slope.
    """
    import shedding_hub as sh

    data = sh.load_dataset("kissler2021viral", local="./data")
    pairs = []
    for participant in data["participants"]:
        by_time = {}
        for m in participant.get("measurements") or []:
            if not isinstance(m.get("value"), (int, float)):
                continue
            by_time.setdefault(m["time"], {})[m["analyte"]] = m["value"]
        for readings in by_time.values():
            if {"AN_OPS_SARSCoV2_ct", "AN_OPS_SARSCoV2_viral"} <= set(readings):
                pairs.append(
                    (
                        np.log10(readings["AN_OPS_SARSCoV2_viral"]),
                        readings["AN_OPS_SARSCoV2_ct"],
                    )
                )
    assert len(pairs) > 500
    log10_conc, ct_values = np.array(pairs).T
    slope, _ = np.polyfit(log10_conc, ct_values, 1)
    # Negative because Ct falls as concentration rises, and near the -3.32 of a
    # perfectly efficient assay.
    assert -5.0 < slope < -2.0
    corr = float(np.corrcoef(log10_conc, ct_values)[0, 1])
    assert corr < -0.9
```

Verified 2026-08-11: the dataset's two analytes are exactly `AN_OPS_SARSCoV2_ct` and `AN_OPS_SARSCoV2_viral`, and the concentration analyte fits under `gamma` (peak_day 0.77, 51 of 68 subjects retained), `gamma_shifted` (−0.52, 60 subjects) and `exponential`. No gate blocks it. If the *Ct* analyte turns out to trip `no_rise_observed`, that is a finding worth reporting rather than working around — it would mean the two scales disagree about whether a rise was observed at all, which the design says cannot happen.

- [ ] **Step 2: Run the tests**

Run: `pytest tests/test_shedding_fit_ct.py -v`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
black tests/test_shedding_fit_ct.py
git add tests/test_shedding_fit_ct.py
git commit -m "test: validate Ct invariance on kissler2021viral"
```

---

### Task 12: Documentation and catalog integrity

**Files:**
- Modify: `docs/modeling-methods.md:140`
- Verify: `shedding_hub/data/shedding_catalog.yaml` unchanged

- [ ] **Step 1: Update the skip-reason table**

Replace the `ct_units` row at `docs/modeling-methods.md:140`:

```markdown
| `ct_units` | 207 | Cycle thresholds. Fittable directly with `fit_shedding_model`, but excluded from the catalog: their heights are cycles below `CT_REFERENCE`, not log10 concentrations, so an ensemble cannot average them together |
```

Add after the table:

```markdown
### Cycle-threshold fits

`Ct = α − β·log10 C` is affine, so a gamma curve in concentration is a gamma
curve in Ct. Cycle-threshold analytes are fitted on `CT_REFERENCE − Ct` —
cycles below a fixed reference of 40, which rises with viral load exactly as a
log10 concentration does.

Because `a0` and `b0` both come back multiplied by the unknown standard-curve
slope, that slope cancels in `b0/a0`. **Peak time, onset and rise duration
therefore compare directly with concentration fits, with no assumption about
PCR efficiency.** Decay rate, half-life and peak height do not:
`SheddingFit.comparable_with` says which is which for any pair of fits.
```

- [ ] **Step 2: Verify the catalog gained and lost no fits**

Task 4 deliberately writes a new, accurate skip message — the old one claimed "neither shedding model applies," which is no longer true — so the rebuilt catalog **will** differ. The diff must be confined to those message strings. Nothing else may move.

```bash
python scripts/build_shedding_catalog.py
```

```bash
python - <<'PY'
import subprocess, yaml

before = yaml.safe_load(
    subprocess.run(
        ["git", "show", "HEAD:shedding_hub/data/shedding_catalog.yaml"],
        capture_output=True, text=True, check=True,
    ).stdout
)
after = yaml.safe_load(open("shedding_hub/data/shedding_catalog.yaml", encoding="utf-8"))

def key(row):
    return (row["dataset_id"], row["analyte"], row["model"], row["reason"])

b_skip = sorted(map(key, before["skipped"]))
a_skip = sorted(map(key, after["skipped"]))
assert b_skip == a_skip, "the set of skipped analytes changed, not just the message"

assert len(before["fits"]) == len(after["fits"]), (
    f"fit count changed: {len(before['fits'])} -> {len(after['fits'])}"
)
assert before["fits"] == after["fits"], "a fit's numbers changed"

n_ct = sum(1 for row in after["skipped"] if row["reason"] == "ct_units")
assert n_ct == 207, f"expected 207 ct_units skips, got {n_ct}"
print(f"catalog OK: {len(after['fits'])} fits unchanged, {n_ct} ct_units skips")
PY
```

Expected: `catalog OK: ... fits unchanged, 207 ct_units skips`. If the fit count moved, Task 4's gate is not holding — that is a stop-and-fix, not a number to update.

Then confirm the textual diff really is message-only:

```bash
git diff --unified=0 shedding_hub/data/shedding_catalog.yaml | grep -E "^[-+]" | grep -v "^[-+][-+]" | grep -vc "message:"
```

Expected: `0`

- [ ] **Step 3: Run everything**

```bash
black --check .
pytest -q
python -m doctest -o ELLIPSIS -o NORMALIZE_WHITESPACE README.md
python -m mkdocs build --strict
```

Expected: all pass

- [ ] **Step 4: Commit**

```bash
git add docs/modeling-methods.md shedding_hub/data/shedding_catalog.yaml
git commit -m "docs: describe cycle-threshold fitting and its comparability limits"
```

---

## Self-review

**Spec coverage.** Response variable and fixed anchor → Tasks 1, 5. Censoring → Task 2. Cutoff resolution incl. the 15 analytes with no numeric limit → Task 2 (fallback branch, tested). Comparability tiering → Task 6. Population coordinate naming → Task 7. Visualisation → Task 9. Ensemble guard → Task 10. Verification plan items 1–4 → Tasks 8 and 11. The spec's `prepare_observations` change → Task 5.

**Spec deviation, deliberate:** the spec's "invert the Ct axis" is wrong given the response is already `40 − Ct`; Task 9 uses a tick formatter instead and says why. The spec should be amended to match once this lands.

**Not in the spec, added here:** Task 4 (catalog gate) and Task 12's catalog-integrity check. Removing the raise without them would silently add up to 207 fits to the published catalog. Task 4 must precede Task 5.

**Type consistency.** `value_type: str` with values `"concentration"` / `"ct"` throughout — `Observations.value_type` (Task 3) → `SheddingFit.value_type` (Task 6) → `plot_fit_diagnostic` (Task 9) → `make_ensemble` (Task 10). `_to_response(value, value_type)` has one signature, used in Tasks 1, 2, 3, 5. `CT_REFERENCE` is defined once in `shedding_fit.py` and imported by `viz.py`.

**Assumptions checked against the code while writing this plan**, so the implementer need not re-verify them: `SheddingFit.median_params()` exists (`shedding_fit.py:760`); `viz.py:10` already imports from `.shedding_fit` and has no `matplotlib.ticker` import; the `kissler2021viral` analyte keys are `AN_OPS_SARSCoV2_ct` and `AN_OPS_SARSCoV2_viral`; `_is_ct_unit` is at `shedding_fit.py:129`; the catalog's skip handler is at `shedding_catalog.py:354`.

Also checked: `kissler2021viral`'s concentration analyte clears every gate under all three models, so Task 11 has a working baseline to compare against.

**Placeholder scan:** none. Every step carries the code it needs, and the two tests that depend on later tasks (`test_catalog_fits_ct_analytes_when_asked`) are explicitly marked xfail in Task 4 and unmarked in Task 6 rather than left dangling.
