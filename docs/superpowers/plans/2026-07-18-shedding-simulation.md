# Shedding Simulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a modeler pick fitted shedding estimates from a table — one study, a chosen subset, or an ensemble across studies — and simulate synthetic infected individuals' shedding trajectories from them.

**Architecture:** Two curve models (exponential decay, gamma) are fitted per analyte by joint maximum likelihood with a left-censored normal likelihood, matching the Rstan tutorial at `shedding-hub.github.io/tutorials/Bayesian-workflow-Rstan.html`. Per-subject log-parameters are summarized into a multivariate normal population distribution; drawing from it produces new individuals. Fits across the repository are precomputed into a shipped YAML catalog.

**Tech Stack:** Python, numpy, scipy (`optimize.minimize`, `stats.norm`), pandas, matplotlib, pyaml, pytest.

**Spec:** `docs/superpowers/specs/2026-07-18-shedding-simulation-design.md`

## Global Constraints

- Formatting: `black` must pass (`black --check .`). CI enforces it. Run `black .` before every commit.
- Python version: the package targets modern Python; use `X | None` union syntax as the existing modules do.
- Model equations, verbatim from the spec (log10 scale, `LN10 = ln(10)`):
  - exponential: `log10 c(t) = (c0 - a0 * t) / LN10`
  - gamma: `log10 c(t) = (c0 + b0 * ln(t) - a0 * t) / LN10`
- `c0` is on the natural-log scale; `a0`, `b0`, `c0` are all strictly positive and modelled as `theta = log(params) ~ MVN(mu, Sigma)`.
- Censored likelihood term is `norm.logcdf((L - mu) / sigma)` where `L` is the log10 censoring limit — the analogue of Stan's `normal_lcdf`.
- `unit` is a hard constraint: never combine or compare estimates across different units.
- Analytes with unit `cycle threshold` are never fitted.
- Existing test convention: `tests/*.py` set `matplotlib.use("Agg")` before other imports and use `pytest` fixtures returning plain dicts.

**File structure** (the spec allowed splitting `simulate.py`; it is split here because the feature spans four distinct responsibilities):

| File | Responsibility |
| --- | --- |
| `shedding_hub/shedding_models.py` | Pure curve math: evaluate models, derived quantities. No I/O. |
| `shedding_hub/shedding_fit.py` | Dataset → observations → censored MLE → `SheddingFit`. |
| `shedding_hub/shedding_catalog.py` | `SheddingCatalog`, `SheddingEnsemble`, `make_ensemble`, `fit_shedding_models`, `load_shedding_catalog`. |
| `shedding_hub/shedding_simulate.py` | `simulate_shedding`, `plot_simulated_shedding`. |

---

### Task 1: Model curve math and dependencies

**Files:**
- Create: `shedding_hub/shedding_models.py`
- Modify: `pyproject.toml`
- Test: `tests/test_shedding_models.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `MODELS: tuple[str, ...]` = `("exponential", "gamma")`
  - `PARAM_NAMES: dict[str, tuple[str, ...]]` = `{"exponential": ("a0", "c0"), "gamma": ("a0", "b0", "c0")}`
  - `LN10: float`
  - `log10_concentration(model: str, params: np.ndarray, times: np.ndarray) -> np.ndarray` — `params` shape `(n, k)`, `times` shape `(m,)`, returns `(n, m)`
  - `log10_concentration_rowwise(model: str, params: np.ndarray, times: np.ndarray) -> np.ndarray` — `params` `(n, k)`, `times` `(n, m)`, returns `(n, m)`; each individual gets its own time row
  - `log10_concentration_pointwise(model: str, params: np.ndarray, times: np.ndarray) -> np.ndarray` — `params` `(n_obs, k)`, `times` `(n_obs,)`, returns `(n_obs,)`
  - `peak_day(model: str, params: np.ndarray) -> np.ndarray`
  - `half_life_days(model: str, params: np.ndarray) -> np.ndarray`
  - `validate_model(model: str) -> None` — raises `ValueError` for unknown model

- [ ] **Step 1: Add dependencies to `pyproject.toml`**

Replace the `dependencies` list:

```toml
dependencies = [
    "pyaml",
    "requests",
    "pandas",
    "matplotlib",
    "numpy",
    "scipy"
]
```

`numpy` is currently only transitive via pandas/matplotlib but is imported directly by this feature, so it must be declared.

- [ ] **Step 2: Write the failing tests**

Create `tests/test_shedding_models.py`:

```python
import numpy as np
import pytest

from shedding_hub.shedding_models import (
    LN10,
    MODELS,
    PARAM_NAMES,
    half_life_days,
    log10_concentration,
    log10_concentration_pointwise,
    log10_concentration_rowwise,
    peak_day,
    validate_model,
)


def test_models_and_param_names():
    assert MODELS == ("exponential", "gamma")
    assert PARAM_NAMES["exponential"] == ("a0", "c0")
    assert PARAM_NAMES["gamma"] == ("a0", "b0", "c0")


def test_exponential_matches_closed_form():
    params = np.array([[0.5, 20.0]])  # a0, c0
    times = np.array([0.0, 2.0, 10.0])
    got = log10_concentration("exponential", params, times)
    expected = (20.0 - 0.5 * times) / LN10
    np.testing.assert_allclose(got[0], expected)


def test_gamma_matches_closed_form():
    params = np.array([[0.5, 2.0, 10.0]])  # a0, b0, c0
    times = np.array([1.0, 4.0])
    got = log10_concentration("gamma", params, times)
    expected = (10.0 + 2.0 * np.log(times) - 0.5 * times) / LN10
    np.testing.assert_allclose(got[0], expected)


def test_gamma_is_nan_for_non_positive_times():
    params = np.array([[0.5, 2.0, 10.0]])
    got = log10_concentration("gamma", params, np.array([-1.0, 0.0, 1.0]))
    assert np.isnan(got[0, 0])
    assert np.isnan(got[0, 1])
    assert np.isfinite(got[0, 2])


def test_exponential_is_finite_for_negative_times():
    params = np.array([[0.5, 20.0]])
    got = log10_concentration("exponential", params, np.array([-3.0]))
    assert np.isfinite(got[0, 0])


def test_gamma_peaks_at_b_over_a():
    params = np.array([[0.5, 2.0, 10.0]])
    assert peak_day("gamma", params)[0] == pytest.approx(4.0)
    # confirm the curve really is maximal there
    times = np.linspace(0.1, 20, 500)
    values = log10_concentration("gamma", params, times)[0]
    assert times[np.argmax(values)] == pytest.approx(4.0, abs=0.1)


def test_exponential_peak_day_is_zero():
    params = np.array([[0.5, 20.0]])
    assert peak_day("exponential", params)[0] == 0.0


def test_half_life():
    params = np.array([[np.log(2.0), 20.0]])
    assert half_life_days("exponential", params)[0] == pytest.approx(1.0)


def test_rowwise_gives_each_individual_its_own_times():
    params = np.array([[0.5, 20.0], [1.0, 20.0]])
    times = np.array([[0.0, 1.0], [0.0, 1.0]])
    got = log10_concentration_rowwise("exponential", params, times)
    assert got.shape == (2, 2)
    np.testing.assert_allclose(got[0, 1], (20.0 - 0.5) / LN10)
    np.testing.assert_allclose(got[1, 1], (20.0 - 1.0) / LN10)


def test_pointwise_pairs_each_observation_with_its_params():
    params = np.array([[0.5, 20.0], [1.0, 20.0]])
    times = np.array([2.0, 2.0])
    got = log10_concentration_pointwise("exponential", params, times)
    assert got.shape == (2,)
    np.testing.assert_allclose(got[0], (20.0 - 1.0) / LN10)
    np.testing.assert_allclose(got[1], (20.0 - 2.0) / LN10)


def test_unknown_model_raises():
    with pytest.raises(ValueError, match="Unknown model"):
        validate_model("weibull")
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/test_shedding_models.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'shedding_hub.shedding_models'`

- [ ] **Step 4: Implement the module**

Create `shedding_hub/shedding_models.py`:

```python
"""
Parametric shedding-curve models.

Two models are supported, both taken from the Shedding Hub Rstan tutorial and
expressed on the log10 scale, which is the scale likelihoods are evaluated on:

- ``exponential``: ``c(t) = c0 * exp(-a0 * t)``, a pure decay appropriate when
  sampling begins at or after peak shedding.
- ``gamma``: ``c(t) = c0 * t**b0 * exp(-a0 * t)``, which rises and falls, peaking
  at ``t = b0 / a0``.

``c0`` is on the natural-log scale, so the log10 concentration at ``t = 0`` is
``c0 / ln(10)``. All parameters are strictly positive.

This module is pure math: no dataset handling, no I/O.
"""

import numpy as np

MODELS = ("exponential", "gamma")

PARAM_NAMES = {
    "exponential": ("a0", "c0"),
    "gamma": ("a0", "b0", "c0"),
}

LN10 = float(np.log(10.0))


def validate_model(model: str) -> None:
    """Raise ``ValueError`` unless ``model`` is a supported model name."""
    if model not in MODELS:
        raise ValueError(f"Unknown model {model!r}. Choose one of {list(MODELS)}.")


def _safe_log(times: np.ndarray) -> np.ndarray:
    """Natural log of ``times``, NaN where ``times <= 0``."""
    positive = times > 0
    return np.where(positive, np.log(np.where(positive, times, 1.0)), np.nan)


def log10_concentration(
    model: str, params: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """
    Evaluate the model for every combination of individual and time.

    Args:
        model: ``"exponential"`` or ``"gamma"``.
        params: Natural-scale parameters, shape ``(n, k)``, ordered as
            ``PARAM_NAMES[model]``.
        times: Times since the reference event, shape ``(m,)``.

    Returns:
        Log10 concentrations, shape ``(n, m)``. Under the gamma model,
        non-positive times yield NaN because ``ln(t)`` is undefined there.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    times = np.asarray(times, dtype=float)
    return log10_concentration_rowwise(
        model, params, np.broadcast_to(times, (params.shape[0], times.size))
    )


def log10_concentration_rowwise(
    model: str, params: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """
    Evaluate the model giving each individual its own row of times.

    Args:
        model: ``"exponential"`` or ``"gamma"``.
        params: Natural-scale parameters, shape ``(n, k)``.
        times: Times, shape ``(n, m)`` — row ``i`` belongs to individual ``i``.

    Returns:
        Log10 concentrations, shape ``(n, m)``.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    times = np.asarray(times, dtype=float)
    if model == "exponential":
        a0 = params[:, 0:1]
        c0 = params[:, 1:2]
        return (c0 - a0 * times) / LN10
    a0 = params[:, 0:1]
    b0 = params[:, 1:2]
    c0 = params[:, 2:3]
    return (c0 + b0 * _safe_log(times) - a0 * times) / LN10


def log10_concentration_pointwise(
    model: str, params: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """
    Evaluate the model once per observation.

    Args:
        model: ``"exponential"`` or ``"gamma"``.
        params: Natural-scale parameters, shape ``(n_obs, k)`` — row ``j`` holds
            the parameters of the subject that observation ``j`` belongs to.
        times: Observation times, shape ``(n_obs,)``.

    Returns:
        Log10 concentrations, shape ``(n_obs,)``.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    times = np.asarray(times, dtype=float)
    if model == "exponential":
        return (params[:, 1] - params[:, 0] * times) / LN10
    return (params[:, 2] + params[:, 1] * _safe_log(times) - params[:, 0] * times) / LN10


def peak_day(model: str, params: np.ndarray) -> np.ndarray:
    """
    Time of peak shedding, in days after the reference event.

    The gamma model peaks at ``b0 / a0``. The exponential model is monotonically
    decreasing, so its maximum is at the reference event and this returns 0.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    if model == "exponential":
        return np.zeros(params.shape[0])
    return params[:, 1] / params[:, 0]


def half_life_days(model: str, params: np.ndarray) -> np.ndarray:
    """
    Half-life of the late-phase decline, ``ln(2) / a0``.

    Exact for the exponential model; asymptotic for the gamma model, whose
    decline approaches rate ``a0`` once ``t`` is well past the peak.
    """
    validate_model(model)
    params = np.atleast_2d(np.asarray(params, dtype=float))
    return np.log(2.0) / params[:, 0]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_shedding_models.py -v`
Expected: PASS, 11 tests.

- [ ] **Step 6: Format and commit**

```bash
black .
git add shedding_hub/shedding_models.py tests/test_shedding_models.py pyproject.toml
git commit -m "feat: add exponential and gamma shedding curve models"
```

---

### Task 2: Extract fittable observations from a dataset

**Files:**
- Create: `shedding_hub/shedding_fit.py`
- Test: `tests/test_shedding_fit.py`

**Interfaces:**
- Consumes: `shedding_hub.shedding_models.{PARAM_NAMES, validate_model}`
- Produces:
  - `Observations` dataclass with fields: `subject_index: np.ndarray` (int), `times: np.ndarray`, `values: np.ndarray` (log10; NaN where censored), `censored: np.ndarray` (bool), `censoring_limit: float`, `subject_ids: list`, `n_subjects: int`, `n_excluded_subjects: int`, `n_dropped_measurements: int`
  - `prepare_observations(dataset: dict, analyte: str, model: str, *, min_observations: int | None = None) -> Observations`
  - `CENSORING_MARGIN: float` = `0.01`
  - `SheddingDataError(ValueError)` — raised when an analyte cannot be fitted; carries a `.reason` string from `{"ct_units", "too_few_subjects", "no_positive_measurements", "unknown_analyte"}`

**Why the censoring-limit fallback matters:** the tutorial hand-set `censlim = 1.96` for the multi-subject fit because one observation (1.97) fell below the declared limit of 2, while the single-subject fit kept `censlim = 2`. The rule below reproduces both automatically.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_shedding_fit.py`:

```python
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from shedding_hub.shedding_fit import (
    SheddingDataError,
    prepare_observations,
)


@pytest.fixture
def simple_dataset():
    """Two subjects, one analyte, with a censored point and a qualitative one."""
    return {
        "dataset_id": "test_study",
        "analytes": {
            "stool": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": 1e6},
                    {"analyte": "stool", "time": 2, "value": 1e5},
                    {"analyte": "stool", "time": 3, "value": "negative"},
                ]
            },
            {
                "measurements": [
                    {"analyte": "stool", "time": 1, "value": 1e7},
                    {"analyte": "stool", "time": 2, "value": 1e6},
                    {"analyte": "stool", "time": 3, "value": "positive"},
                    {"analyte": "stool", "time": "unknown", "value": 1e4},
                ]
            },
        ],
    }


def test_extracts_values_on_log10_scale(simple_dataset):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.n_subjects == 2
    np.testing.assert_allclose(obs.values[0], 6.0)
    np.testing.assert_allclose(obs.values[1], 5.0)


def test_negative_becomes_censored(simple_dataset):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censored.sum() == 1
    assert np.isnan(obs.values[obs.censored][0])


def test_censoring_limit_uses_declared_loq_when_below_smallest_positive(
    simple_dataset,
):
    obs = prepare_observations(simple_dataset, "stool", "exponential")
    # LOQ 100 -> log10 = 2, which is below the smallest observed positive (5.0)
    assert obs.censoring_limit == pytest.approx(2.0)


def test_censoring_limit_falls_back_below_smallest_positive(simple_dataset):
    # Declare a limit above every observed value; the fallback must kick in.
    simple_dataset["analytes"]["stool"]["limit_of_quantification"] = 1e8
    with pytest.warns(UserWarning, match="censoring limit"):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.censoring_limit == pytest.approx(5.0 - 0.01)


def test_qualitative_and_unknown_time_are_dropped_with_warning(simple_dataset):
    with pytest.warns(UserWarning):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    # 7 measurements total, minus 1 qualitative "positive", minus 1 unknown time
    assert obs.times.size == 5
    assert obs.n_dropped_measurements == 2


def test_ct_analyte_is_rejected():
    dataset = {
        "dataset_id": "ct_study",
        "analytes": {
            "swab": {
                "specimen": "saliva",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "cycle threshold",
                "limit_of_quantification": "unknown",
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "swab", "time": 1, "value": 20.0},
                    {"analyte": "swab", "time": 2, "value": 25.0},
                ]
            }
        ],
    }
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(dataset, "swab", "exponential")
    assert excinfo.value.reason == "ct_units"


def test_gamma_drops_non_positive_times(simple_dataset):
    simple_dataset["participants"][0]["measurements"].append(
        {"analyte": "stool", "time": 0, "value": 1e6}
    )
    simple_dataset["participants"][0]["measurements"].append(
        {"analyte": "stool", "time": -2, "value": 1e6}
    )
    with pytest.warns(UserWarning):
        gamma_obs = prepare_observations(simple_dataset, "stool", "gamma")
    exponential_obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert (gamma_obs.times > 0).all()
    assert exponential_obs.times.size == gamma_obs.times.size + 2


def test_subjects_below_min_observations_are_excluded(simple_dataset):
    simple_dataset["participants"].append(
        {"measurements": [{"analyte": "stool", "time": 1, "value": 1e5}]}
    )
    with pytest.warns(UserWarning, match="excluded"):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert obs.n_subjects == 2
    assert obs.n_excluded_subjects == 1


def test_no_usable_subject_raises(simple_dataset):
    for participant in simple_dataset["participants"]:
        participant["measurements"] = participant["measurements"][:1]
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(simple_dataset, "stool", "exponential")
    assert excinfo.value.reason == "too_few_subjects"


def test_unknown_analyte_raises(simple_dataset):
    with pytest.raises(SheddingDataError) as excinfo:
        prepare_observations(simple_dataset, "sputum", "exponential")
    assert excinfo.value.reason == "unknown_analyte"


def test_subject_index_is_contiguous_after_exclusions(simple_dataset):
    simple_dataset["participants"].insert(
        0, {"measurements": [{"analyte": "stool", "time": 1, "value": 1e5}]}
    )
    with pytest.warns(UserWarning):
        obs = prepare_observations(simple_dataset, "stool", "exponential")
    assert set(obs.subject_index.tolist()) == {0, 1}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_shedding_fit.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'shedding_hub.shedding_fit'`

- [ ] **Step 3: Implement observation extraction**

Create `shedding_hub/shedding_fit.py`:

```python
"""
Fit shedding-curve models to Shedding Hub datasets.

Fitting is done per analyte by joint maximum likelihood over every subject's
parameters plus one shared measurement-error standard deviation, using a
left-censored normal likelihood so that ``negative`` measurements contribute the
information they carry (that the concentration was below the limit) instead of
being discarded. Roughly 37% of measurements in the repository are ``negative``;
dropping them biases decay rates slow and inflates simulated late-phase shedding.
"""

import math
import warnings
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .shedding_models import PARAM_NAMES, validate_model

CENSORING_MARGIN = 0.01
NEGATIVE_VALUE = "negative"
QUALITATIVE_VALUES = (
    "positive",
    "weak positive",
    "strong positive",
    "inconclusive",
)


class SheddingDataError(ValueError):
    """
    Raised when an analyte cannot be fitted.

    Attributes:
        reason: Machine-readable cause, one of ``ct_units``, ``too_few_subjects``,
            ``no_positive_measurements``, ``unknown_analyte``. The catalog builder
            records this so a missing study reads as unsuitable, not as a bug.
    """

    def __init__(self, message: str, reason: str):
        super().__init__(message)
        self.reason = reason


@dataclass
class Observations:
    """Model-ready observations for a single analyte."""

    subject_index: np.ndarray
    times: np.ndarray
    values: np.ndarray
    censored: np.ndarray
    censoring_limit: float
    subject_ids: list = field(default_factory=list)
    n_subjects: int = 0
    n_excluded_subjects: int = 0
    n_dropped_measurements: int = 0


def _is_ct_unit(unit: Any) -> bool:
    if unit is None:
        return False
    lowered = str(unit).strip().lower()
    return "cycle" in lowered or lowered == "ct"


def _numeric_limit(value: Any) -> float | None:
    """Return a positive numeric limit, or None if unusable (e.g. ``unknown``)."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value) if value > 0 else None


def _resolve_censoring_limit(
    analyte_spec: dict, observed_log10: np.ndarray
) -> float:
    """
    Resolve the log10 censoring limit.

    Prefers the declared limit of quantification, then the limit of detection.
    Falls back to just below the smallest observed positive value when neither is
    usable, or when the declared limit is not strictly below every observation —
    the likelihood is only coherent if censored points really do sit below the
    limit.
    """
    declared = None
    for key in ("limit_of_quantification", "limit_of_detection"):
        limit = _numeric_limit(analyte_spec.get(key))
        if limit is not None:
            declared = math.log10(limit)
            break

    if observed_log10.size == 0:
        if declared is None:
            raise SheddingDataError(
                "No positive measurements and no declared limit, so the censoring "
                "limit cannot be resolved.",
                "no_positive_measurements",
            )
        return declared

    smallest = float(observed_log10.min())
    if declared is None or declared >= smallest:
        fallback = smallest - CENSORING_MARGIN
        warnings.warn(
            "Falling back to a censoring limit of "
            f"{fallback:.4g} (log10) because the declared limit "
            f"({'none' if declared is None else f'{declared:.4g}'}) is missing or "
            f"not below the smallest observed value ({smallest:.4g}).",
            UserWarning,
            stacklevel=2,
        )
        return fallback
    return declared


def prepare_observations(
    dataset: dict,
    analyte: str,
    model: str,
    *,
    min_observations: int | None = None,
) -> Observations:
    """
    Extract model-ready observations for one analyte of one dataset.

    Args:
        dataset: Dataset dictionary from ``load_dataset``.
        analyte: Key into ``dataset["analytes"]``.
        model: ``"exponential"`` or ``"gamma"``.
        min_observations: Minimum usable measurements a subject must have to be
            retained. Defaults to the number of per-subject parameters (3 for
            gamma, 2 for exponential). ``sigma`` is shared across subjects, so a
            subject does not need residual degrees of freedom of its own.

    Returns:
        An ``Observations`` instance with subject indices renumbered contiguously
        from zero over the retained subjects.

    Raises:
        SheddingDataError: The analyte is unknown, uses cycle-threshold units, has
            no positive measurements, or leaves no subject with enough data.
    """
    validate_model(model)
    if not dataset or not isinstance(dataset, dict):
        raise ValueError("Dataset must be a non-empty dictionary")
    for key in ("analytes", "participants"):
        if key not in dataset:
            raise ValueError(f"Dataset missing required key: {key}")

    analytes = dataset["analytes"]
    if analyte not in analytes:
        raise SheddingDataError(
            f"Analyte {analyte!r} not in dataset; available: {sorted(analytes)}.",
            "unknown_analyte",
        )
    analyte_spec = analytes[analyte]

    if _is_ct_unit(analyte_spec.get("unit")):
        raise SheddingDataError(
            f"Analyte {analyte!r} is reported in {analyte_spec.get('unit')!r}. "
            "Cycle-threshold values are inversely related to concentration and "
            "already on a log scale, so neither shedding model applies. Select a "
            "concentration analyte instead.",
            "ct_units",
        )

    if min_observations is None:
        min_observations = len(PARAM_NAMES[model])

    per_subject: list[dict[str, list]] = []
    subject_ids: list = []
    n_dropped = 0

    for position, participant in enumerate(dataset["participants"]):
        times: list[float] = []
        values: list[float] = []
        censored: list[bool] = []
        limits: list[float | None] = []
        for measurement in participant.get("measurements") or []:
            if measurement.get("analyte") != analyte:
                continue
            time = measurement.get("time")
            if not isinstance(time, (int, float)) or isinstance(time, bool):
                n_dropped += 1
                continue
            time = float(time)
            if model == "gamma" and time <= 0:
                n_dropped += 1
                continue
            value = measurement.get("value")
            if isinstance(value, str):
                if value == NEGATIVE_VALUE:
                    times.append(time)
                    values.append(np.nan)
                    censored.append(True)
                    limits.append(
                        _numeric_limit(measurement.get("limit_of_quantification"))
                    )
                else:
                    n_dropped += 1
                continue
            if not isinstance(value, (int, float)) or value <= 0:
                n_dropped += 1
                continue
            times.append(time)
            values.append(math.log10(float(value)))
            censored.append(False)
            limits.append(None)

        if times:
            per_subject.append(
                {
                    "times": times,
                    "values": values,
                    "censored": censored,
                    "limits": limits,
                }
            )
            subject_ids.append(participant.get("patient_id", position + 1))

    retained = [s for s in per_subject if len(s["times"]) >= min_observations]
    n_excluded = len(per_subject) - len(retained)
    retained_ids = [
        subject_id
        for subject_id, subject in zip(subject_ids, per_subject)
        if len(subject["times"]) >= min_observations
    ]

    if n_excluded:
        warnings.warn(
            f"{n_excluded} subject(s) excluded from the {analyte!r} fit for having "
            f"fewer than {min_observations} usable measurements.",
            UserWarning,
            stacklevel=2,
        )
    if n_dropped:
        warnings.warn(
            f"{n_dropped} measurement(s) dropped from the {analyte!r} fit "
            "(qualitative result, unknown time, or a non-positive time under the "
            "gamma model).",
            UserWarning,
            stacklevel=2,
        )
    if not retained:
        raise SheddingDataError(
            f"No subject has at least {min_observations} usable measurements for "
            f"analyte {analyte!r}.",
            "too_few_subjects",
        )

    subject_index = np.concatenate(
        [np.full(len(s["times"]), i, dtype=int) for i, s in enumerate(retained)]
    )
    times_array = np.concatenate([np.asarray(s["times"], float) for s in retained])
    values_array = np.concatenate([np.asarray(s["values"], float) for s in retained])
    censored_array = np.concatenate(
        [np.asarray(s["censored"], bool) for s in retained]
    )

    observed = values_array[~censored_array]
    if observed.size == 0:
        raise SheddingDataError(
            f"Analyte {analyte!r} has no positive measurements to fit.",
            "no_positive_measurements",
        )

    censoring_limit = _resolve_censoring_limit(analyte_spec, observed)

    return Observations(
        subject_index=subject_index,
        times=times_array,
        values=values_array,
        censored=censored_array,
        censoring_limit=censoring_limit,
        subject_ids=retained_ids,
        n_subjects=len(retained),
        n_excluded_subjects=n_excluded,
        n_dropped_measurements=n_dropped,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_shedding_fit.py -v`
Expected: PASS, 11 tests.

- [ ] **Step 5: Format and commit**

```bash
black .
git add shedding_hub/shedding_fit.py tests/test_shedding_fit.py
git commit -m "feat: extract model-ready observations with censoring from datasets"
```

---

### Task 3: Censored maximum-likelihood fitting

**Files:**
- Modify: `shedding_hub/shedding_fit.py` (append)
- Test: `tests/test_shedding_fit.py` (append)

**Interfaces:**
- Consumes: `Observations`, `prepare_observations`, `shedding_models` helpers.
- Produces:
  - `SheddingFit` dataclass with fields `model`, `method`, `population_mean` (`np.ndarray`, shape `(k,)`), `population_cov` (`(k, k)`), `sigma` (float), `subject_params` (`pd.DataFrame | None`), `censoring_limit`, `dataset_id`, `analyte`, `biomarker`, `specimen`, `reference_event`, `unit`, `gene_target`, `dose`, `vaccine_type`, `n_subjects`, `n_measurements`, `n_censored`, `n_excluded_subjects`, `n_dropped_measurements`, `converged`, `log_likelihood`, `aic`
  - properties `param_names`, `median_params` (`np.ndarray`, `exp(population_mean)`), `peak_day`, `peak_log10`, `half_life_days`
  - method `sample_params(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray]` returning natural-scale params `(n, k)` and a `(n,)` array of source dataset ids
  - `fit_shedding_model(dataset: dict, *, analyte: str, model: str = "gamma", min_observations: int | None = None) -> SheddingFit`

**Note on AIC:** `aic` compares models fitted to the *same* observations. The gamma model drops non-positive times while the exponential model keeps them, so for datasets containing such times the two rows are not directly comparable. The `n_measurements` column exposes this; the docstring says so.

- [ ] **Step 1: Create the shared synthetic-dataset fixture**

Later tasks need the same builder. Put it in `tests/conftest.py` as a factory
fixture rather than importing across test modules — `tests/` has no
`__init__.py`, so `from tests.test_shedding_fit import ...` would not resolve.

Create `tests/conftest.py`:

```python
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest


@pytest.fixture
def make_synthetic_dataset():
    """
    Build a dataset by simulating from known population parameters.

    Returns a factory so tests can vary the truth they fit against. Values below
    ``loq`` are written as ``negative``, reproducing real left-censoring.
    """

    def _make(
        model,
        mu,
        cov,
        sigma=0.3,
        n_subjects=40,
        seed=0,
        times=None,
        loq=1e2,
        dataset_id="synthetic",
    ):
        from shedding_hub.shedding_models import log10_concentration

        rng = np.random.default_rng(seed)
        times = np.arange(1.0, 15.0) if times is None else np.asarray(times, float)
        theta = rng.multivariate_normal(np.asarray(mu, float), np.asarray(cov, float), size=n_subjects)
        truth = log10_concentration(model, np.exp(theta), times)
        noisy = truth + rng.normal(0.0, sigma, size=truth.shape)
        limit = np.log10(loq)

        participants = []
        for row in noisy:
            measurements = []
            for time, value in zip(times, row):
                if value < limit:
                    measurements.append(
                        {"analyte": "stool", "time": float(time), "value": "negative"}
                    )
                else:
                    measurements.append(
                        {
                            "analyte": "stool",
                            "time": float(time),
                            "value": float(10.0**value),
                        }
                    )
            participants.append({"measurements": measurements})

        return {
            "dataset_id": dataset_id,
            "analytes": {
                "stool": {
                    "specimen": "stool",
                    "biomarker": "SARS-CoV-2",
                    "reference_event": "symptom onset",
                    "unit": "gc/mL",
                    "limit_of_quantification": loq,
                    "limit_of_detection": "unknown",
                }
            },
            "participants": participants,
        }

    return _make
```

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_shedding_fit.py`:

```python
from shedding_hub.shedding_fit import SheddingFit, fit_shedding_model


def test_recovers_known_exponential_population(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.09, 0.04]), sigma=0.3, n_subjects=60
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    np.testing.assert_allclose(fit.population_mean, mu, atol=0.15)
    assert fit.sigma == pytest.approx(0.3, abs=0.15)
    assert fit.converged


def test_recovers_known_gamma_population(make_synthetic_dataset):
    mu = np.array([np.log(0.5), np.log(1.5), np.log(12.0)])
    dataset = make_synthetic_dataset(
        "gamma", mu, np.diag([0.04, 0.04, 0.04]), sigma=0.3, n_subjects=60
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="gamma")
    np.testing.assert_allclose(fit.population_mean, mu, atol=0.25)
    assert fit.converged


def test_censored_fit_beats_dropping_negatives(make_synthetic_dataset):
    """The property the whole design rests on."""
    mu = np.array([np.log(0.6), np.log(14.0)])
    # A high limit censors much of the decay phase.
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), sigma=0.2, n_subjects=60, loq=1e4
    )

    censored_fit = fit_shedding_model(dataset, analyte="stool", model="exponential")

    dropped = {
        **dataset,
        "participants": [
            {
                "measurements": [
                    m for m in p["measurements"] if m["value"] != "negative"
                ]
            }
            for p in dataset["participants"]
        ],
    }
    naive_fit = fit_shedding_model(dropped, analyte="stool", model="exponential")

    true_decay = mu[0]
    censored_error = abs(censored_fit.population_mean[0] - true_decay)
    naive_error = abs(naive_fit.population_mean[0] - true_decay)
    assert censored_error < naive_error
    # Dropping negatives should understate the decay rate.
    assert naive_fit.population_mean[0] < censored_fit.population_mean[0]


def test_fit_carries_metadata_and_counts(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=10
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert isinstance(fit, SheddingFit)
    assert fit.dataset_id == "synthetic"
    assert fit.analyte == "stool"
    assert fit.biomarker == "SARS-CoV-2"
    assert fit.specimen == "stool"
    assert fit.reference_event == "symptom onset"
    assert fit.unit == "gc/mL"
    assert fit.method == "mle"
    assert fit.n_subjects == 10
    # 10 subjects x 14 time points, none dropped by the exponential model.
    assert fit.n_measurements == 140
    assert 0 < fit.n_censored < fit.n_measurements
    assert fit.param_names == ("a0", "c0")


def test_median_params_are_exp_of_population_mean(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=10
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    np.testing.assert_allclose(fit.median_params, np.exp(fit.population_mean))


def test_gamma_peak_day_is_b_over_a(make_synthetic_dataset):
    mu = np.array([np.log(0.5), np.log(2.0), np.log(12.0)])
    dataset = make_synthetic_dataset(
        "gamma", mu, np.diag([0.04, 0.04, 0.04]), n_subjects=20
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="gamma")
    assert fit.peak_day == pytest.approx(fit.median_params[1] / fit.median_params[0])


def test_sample_params_shape_and_source(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=20
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    rng = np.random.default_rng(1)
    params, sources = fit.sample_params(rng, 25)
    assert params.shape == (25, 2)
    assert (params > 0).all()
    assert sources.shape == (25,)
    assert set(sources.tolist()) == {"synthetic"}


def test_subject_params_has_one_row_per_subject(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=12
    )
    fit = fit_shedding_model(dataset, analyte="stool", model="exponential")
    assert len(fit.subject_params) == 12
    assert list(fit.subject_params.columns) == ["subject_id", "a0", "c0"]
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/test_shedding_fit.py -v -k "recovers or censored_fit or metadata or median_params or peak_day or sample_params or subject_params"`
Expected: FAIL — `ImportError: cannot import name 'fit_shedding_model'`

- [ ] **Step 4: Implement fitting**

Append to `shedding_hub/shedding_fit.py`:

```python
import pandas as pd
from scipy import optimize
from scipy.stats import norm

from .shedding_models import (
    half_life_days,
    log10_concentration,
    log10_concentration_pointwise,
    peak_day,
)

_MIN_PARAM = 1e-6
_THETA_BOUNDS = (-25.0, 25.0)
_LOG_SIGMA_BOUNDS = (-10.0, 5.0)


@dataclass
class SheddingFit:
    """
    A fitted shedding model for one analyte of one dataset.

    ``population_mean`` and ``population_cov`` describe the distribution of
    ``theta = log(params)`` across subjects; drawing from that multivariate normal
    and exponentiating produces a new plausible individual, which is what makes
    simulation possible.

    Two-stage estimation does not shrink individual estimates toward the
    population mean the way a hierarchical Bayesian fit does, so
    ``population_cov`` absorbs within-subject estimation error and overestimates
    true between-subject variance. Simulated cohorts are therefore somewhat more
    dispersed than reality, the more so when subjects have few observations.
    """

    model: str
    method: str
    population_mean: np.ndarray
    population_cov: np.ndarray
    sigma: float
    subject_params: pd.DataFrame | None
    censoring_limit: float
    dataset_id: str
    analyte: str
    biomarker: str | None
    specimen: str | None
    reference_event: str | None
    unit: str | None
    gene_target: str | None
    dose: int | None
    vaccine_type: str | None
    n_subjects: int
    n_measurements: int
    n_censored: int
    n_excluded_subjects: int
    n_dropped_measurements: int
    converged: bool
    log_likelihood: float
    aic: float

    @property
    def param_names(self) -> tuple[str, ...]:
        return PARAM_NAMES[self.model]

    @property
    def median_params(self) -> np.ndarray:
        """
        Parameters of the median individual.

        Because ``theta`` is normal, the parameters are lognormal, so
        ``exp(population_mean)`` is exactly their median — not their mean. Note
        that the median individual's trajectory is not the population's mean
        trajectory; to aggregate load across a cohort, simulate rather than
        scaling this up.
        """
        return np.exp(self.population_mean)

    @property
    def peak_day(self) -> float:
        return float(peak_day(self.model, self.median_params[None, :])[0])

    @property
    def peak_log10(self) -> float:
        """Log10 concentration of the median individual at its peak."""
        if self.model == "exponential":
            return float(
                log10_concentration(
                    self.model, self.median_params[None, :], np.array([0.0])
                )[0, 0]
            )
        return float(
            log10_concentration(
                self.model, self.median_params[None, :], np.array([self.peak_day])
            )[0, 0]
        )

    @property
    def half_life_days(self) -> float:
        return float(half_life_days(self.model, self.median_params[None, :])[0])

    def sample_params(
        self, rng: np.random.Generator, n: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Draw ``n`` individuals' natural-scale parameters.

        Returns:
            ``(params, sources)`` where ``params`` has shape ``(n, k)`` and
            ``sources`` names the dataset each individual came from — constant
            here, but varying for a mixture ensemble, so both share this
            interface.
        """
        if n < 1:
            raise ValueError("n_individuals must be at least 1")
        theta = rng.multivariate_normal(self.population_mean, self.population_cov, n)
        return np.exp(theta), np.full(n, self.dataset_id, dtype=object)


def _initial_theta(model: str, observations: Observations) -> np.ndarray:
    """
    Initialize per-subject log-parameters by ordinary least squares.

    Fits each subject's uncensored points; subjects with too few uncensored
    points fall back to a pooled fit across all subjects, so an all-censored or
    nearly all-censored subject still starts somewhere sensible.
    """
    k = len(PARAM_NAMES[model])
    uncensored = ~observations.censored

    def design(times: np.ndarray) -> np.ndarray:
        if model == "exponential":
            return np.column_stack([np.ones_like(times), -times])
        return np.column_stack([np.ones_like(times), np.log(times), -times])

    def solve(times: np.ndarray, values: np.ndarray) -> np.ndarray:
        coefficients, *_ = np.linalg.lstsq(design(times), values * LN10, rcond=None)
        if model == "exponential":
            c0, a0 = coefficients
            params = np.array([a0, c0])
        else:
            c0, b0, a0 = coefficients
            params = np.array([a0, b0, c0])
        return np.log(np.clip(params, _MIN_PARAM, None))

    pooled = solve(
        observations.times[uncensored], observations.values[uncensored]
    )

    theta = np.tile(pooled, (observations.n_subjects, 1))
    for i in range(observations.n_subjects):
        mask = uncensored & (observations.subject_index == i)
        if mask.sum() >= k:
            try:
                theta[i] = solve(observations.times[mask], observations.values[mask])
            except np.linalg.LinAlgError:
                pass
    return np.clip(theta, *_THETA_BOUNDS)


def _negative_log_likelihood(
    x: np.ndarray, model: str, observations: Observations
) -> float:
    """
    Negative log likelihood with left-censored observations.

    Uncensored points contribute a normal density; censored points contribute
    ``Phi((L - mu) / sigma)``, the probability of falling below the limit. This is
    the direct analogue of Stan's ``normal_lcdf`` term.
    """
    k = len(PARAM_NAMES[model])
    n = observations.n_subjects
    theta = x[: n * k].reshape(n, k)
    log_sigma = x[-1]
    sigma = math.exp(log_sigma)

    params = np.exp(theta)[observations.subject_index]
    predicted = log10_concentration_pointwise(model, params, observations.times)
    if not np.all(np.isfinite(predicted)):
        return np.inf

    total = 0.0
    uncensored = ~observations.censored
    if uncensored.any():
        residual = (observations.values[uncensored] - predicted[uncensored]) / sigma
        total += 0.5 * float(np.sum(residual**2))
        total += float(uncensored.sum()) * (log_sigma + 0.5 * math.log(2 * math.pi))
    if observations.censored.any():
        z = (observations.censoring_limit - predicted[observations.censored]) / sigma
        total -= float(np.sum(norm.logcdf(z)))
    return total if np.isfinite(total) else np.inf


def fit_shedding_model(
    dataset: dict,
    *,
    analyte: str,
    model: str = "gamma",
    min_observations: int | None = None,
) -> SheddingFit:
    """
    Fit a shedding model to one analyte by censored maximum likelihood.

    Args:
        dataset: Dataset dictionary from ``load_dataset``.
        analyte: Key into ``dataset["analytes"]``.
        model: ``"exponential"`` or ``"gamma"``.
        min_observations: Minimum usable measurements per subject; defaults to the
            number of per-subject parameters.

    Returns:
        A ``SheddingFit``.

    Raises:
        SheddingDataError: The analyte cannot be fitted (see ``reason``).

    Note:
        ``aic`` is only comparable between models fitted to the same
        observations. The gamma model drops non-positive times while the
        exponential model keeps them, so compare ``n_measurements`` before
        comparing ``aic`` across models.
    """
    validate_model(model)
    observations = prepare_observations(
        dataset, analyte, model, min_observations=min_observations
    )

    k = len(PARAM_NAMES[model])
    n = observations.n_subjects
    x0 = np.concatenate([_initial_theta(model, observations).ravel(), [math.log(0.5)]])
    bounds = [_THETA_BOUNDS] * (n * k) + [_LOG_SIGMA_BOUNDS]

    result = optimize.minimize(
        _negative_log_likelihood,
        x0,
        args=(model, observations),
        method="L-BFGS-B",
        bounds=bounds,
    )
    if not result.success:
        warnings.warn(
            f"Optimizer did not converge for analyte {analyte!r} "
            f"({result.message}). The fit is returned with converged=False.",
            UserWarning,
            stacklevel=2,
        )

    theta = result.x[: n * k].reshape(n, k)
    sigma = float(np.exp(result.x[-1]))
    population_mean = theta.mean(axis=0)
    population_cov = (
        np.cov(theta, rowvar=False, ddof=1) if n > 1 else np.zeros((k, k))
    )
    population_cov = np.atleast_2d(population_cov)

    subject_params = pd.DataFrame(
        np.exp(theta), columns=list(PARAM_NAMES[model])
    )
    subject_params.insert(0, "subject_id", observations.subject_ids)

    log_likelihood = -float(result.fun)
    n_parameters = n * k + 1
    analyte_spec = dataset["analytes"][analyte]
    specimen = analyte_spec.get("specimen")
    if isinstance(specimen, list):
        specimen = "+".join(specimen)

    return SheddingFit(
        model=model,
        method="mle",
        population_mean=population_mean,
        population_cov=population_cov,
        sigma=sigma,
        subject_params=subject_params,
        censoring_limit=observations.censoring_limit,
        dataset_id=dataset.get("dataset_id", "unknown"),
        analyte=analyte,
        biomarker=analyte_spec.get("biomarker"),
        specimen=specimen,
        reference_event=analyte_spec.get("reference_event"),
        unit=analyte_spec.get("unit"),
        gene_target=analyte_spec.get("gene_target"),
        dose=analyte_spec.get("dose"),
        vaccine_type=analyte_spec.get("vaccine_type"),
        n_subjects=n,
        n_measurements=int(observations.times.size),
        n_censored=int(observations.censored.sum()),
        n_excluded_subjects=observations.n_excluded_subjects,
        n_dropped_measurements=observations.n_dropped_measurements,
        converged=bool(result.success),
        log_likelihood=log_likelihood,
        aic=2.0 * n_parameters - 2.0 * log_likelihood,
    )
```

Also add `LN10` to the imports at the top of the file — change the existing import line to:

```python
from .shedding_models import LN10, PARAM_NAMES, validate_model
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_shedding_fit.py -v`
Expected: PASS, all tests including the eight new ones.

- [ ] **Step 6: Format and commit**

```bash
black .
git add shedding_hub/shedding_fit.py tests/test_shedding_fit.py tests/conftest.py
git commit -m "feat: fit shedding models by censored maximum likelihood"
```

---

### Task 4: Validate the port against the published Rstan tutorial

**Files:**
- Create: `tests/test_shedding_tutorial_agreement.py`

**Interfaces:**
- Consumes: `fit_shedding_model`, `shedding_hub.load_dataset`.
- Produces: nothing (validation only).

**Context:** The tutorial fits subject 3 of `woelfel2020virological` stool with a censored exponential model and reports posterior means `a0 = 0.74`, `c0 = 20.37`, `sig_obs = 0.92`. Priors there are flat (`normal(0, 100)`), so the maximum likelihood estimate should land close. The repository data matches the tutorial exactly: 14 positives, negatives at t = 20, 22, 23, and `limit_of_quantification: 100` (log10 = 2, the tutorial's `censorlimit`). This test is the strongest available evidence the Python port is faithful.

- [ ] **Step 1: Write the failing test**

Create `tests/test_shedding_tutorial_agreement.py`:

```python
"""
Validate the Python port against the published Rstan tutorial.

Reference: https://shedding-hub.github.io/tutorials/Bayesian-workflow-Rstan.html
Subject 3 of woelfel2020virological, stool, censored exponential model, reported
posterior means a0 = 0.74, c0 = 20.37, sig_obs = 0.92 under flat priors.
"""

import matplotlib

matplotlib.use("Agg")

import pathlib

import pytest

import shedding_hub as sh
from shedding_hub.shedding_fit import fit_shedding_model

DATA = pathlib.Path(__file__).parent.parent / "data"


@pytest.fixture
def woelfel_subject_3():
    dataset = sh.load_dataset("woelfel2020virological", local=str(DATA))
    dataset["participants"] = [dataset["participants"][2]]
    return dataset


def test_subject_3_data_matches_the_tutorial(woelfel_subject_3):
    """Guard the fixture: if the dataset changes, the comparison is void."""
    measurements = [
        m
        for m in woelfel_subject_3["participants"][0]["measurements"]
        if m["analyte"] == "stool"
    ]
    positives = [m for m in measurements if m["value"] != "negative"]
    negatives = [m for m in measurements if m["value"] == "negative"]
    assert len(positives) == 14
    assert sorted(m["time"] for m in negatives) == [20, 22, 23]
    assert woelfel_subject_3["analytes"]["stool"]["limit_of_quantification"] == 100


def test_exponential_fit_agrees_with_published_posterior(woelfel_subject_3):
    fit = fit_shedding_model(
        woelfel_subject_3, analyte="stool", model="exponential"
    )
    a0, c0 = fit.median_params
    assert a0 == pytest.approx(0.74, abs=0.15)
    assert c0 == pytest.approx(20.37, abs=2.0)
    assert fit.sigma == pytest.approx(0.92, abs=0.3)
    assert fit.censoring_limit == pytest.approx(2.0)
```

- [ ] **Step 2: Run the test**

Run: `pytest tests/test_shedding_tutorial_agreement.py -v`
Expected: PASS. If it fails, the port is wrong — do not loosen the tolerances to make it pass. Check, in order: that `c0` is treated as natural-log scale (`log10 c = c0 / ln(10)`), that censored points use `norm.logcdf` rather than being dropped, and that `censoring_limit` resolved to 2.0 rather than falling back.

- [ ] **Step 3: Commit**

```bash
black .
git add tests/test_shedding_tutorial_agreement.py
git commit -m "test: validate censored fit against published Rstan tutorial"
```

---

### Task 5: Simulate individuals from a fit

**Files:**
- Create: `shedding_hub/shedding_simulate.py`
- Test: `tests/test_shedding_simulate.py`

**Interfaces:**
- Consumes: `SheddingFit.sample_params`, `shedding_models.log10_concentration_rowwise`.
- Produces: `simulate_shedding(source, *, n_individuals: int, times, incubation_period=None, include_measurement_error: bool = False, seed=None) -> pd.DataFrame` with columns `individual_id`, `time`, `log10_value`, `value`, `detected`, `source_dataset_id`.

`source` is anything exposing `.model`, `.censoring_limit`, `.sigma`, `.reference_event`, `.unit`, and `.sample_params(rng, n)` — satisfied by `SheddingFit` now and `SheddingEnsemble` in Task 7.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_shedding_simulate.py`:

```python
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from shedding_hub.shedding_fit import fit_shedding_model
from shedding_hub.shedding_simulate import simulate_shedding


@pytest.fixture
def exponential_fit(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    dataset = make_synthetic_dataset(
        "exponential", mu, np.diag([0.04, 0.04]), n_subjects=25
    )
    return fit_shedding_model(dataset, analyte="stool", model="exponential")


@pytest.fixture
def gamma_fit(make_synthetic_dataset):
    mu = np.array([np.log(0.5), np.log(2.0), np.log(12.0)])
    dataset = make_synthetic_dataset(
        "gamma", mu, np.diag([0.04, 0.04, 0.04]), n_subjects=25
    )
    return fit_shedding_model(dataset, analyte="stool", model="gamma")


def test_returns_tidy_frame_of_expected_shape(exponential_fit):
    times = np.arange(0.0, 10.0)
    traj = simulate_shedding(
        exponential_fit, n_individuals=20, times=times, seed=0
    )
    assert isinstance(traj, pd.DataFrame)
    assert list(traj.columns) == [
        "individual_id",
        "time",
        "log10_value",
        "value",
        "detected",
        "source_dataset_id",
    ]
    assert len(traj) == 20 * len(times)
    assert traj["individual_id"].nunique() == 20


def test_is_reproducible_with_a_seed(exponential_fit):
    a = simulate_shedding(exponential_fit, n_individuals=10, times=[1, 2], seed=7)
    b = simulate_shedding(exponential_fit, n_individuals=10, times=[1, 2], seed=7)
    c = simulate_shedding(exponential_fit, n_individuals=10, times=[1, 2], seed=8)
    pd.testing.assert_frame_equal(a, b)
    assert not np.allclose(a["log10_value"], c["log10_value"])


def test_measurement_error_is_off_by_default(exponential_fit):
    clean = simulate_shedding(
        exponential_fit, n_individuals=200, times=[3.0], seed=3
    )
    noisy = simulate_shedding(
        exponential_fit,
        n_individuals=200,
        times=[3.0],
        seed=3,
        include_measurement_error=True,
    )
    assert noisy["log10_value"].std() > clean["log10_value"].std()


def test_detected_reflects_the_censoring_limit(exponential_fit):
    traj = simulate_shedding(
        exponential_fit, n_individuals=50, times=np.arange(0.0, 40.0), seed=1
    )
    finite = traj.dropna(subset=["log10_value"])
    expected = finite["log10_value"] >= exponential_fit.censoring_limit
    pd.testing.assert_series_equal(
        finite["detected"], expected, check_names=False
    )


def test_values_below_the_limit_are_not_clipped(exponential_fit):
    traj = simulate_shedding(
        exponential_fit, n_individuals=50, times=np.arange(0.0, 60.0), seed=2
    )
    below = traj[~traj["detected"]].dropna(subset=["log10_value"])
    assert (below["log10_value"] < exponential_fit.censoring_limit).all()
    assert below["log10_value"].min() < exponential_fit.censoring_limit - 1


def test_gamma_non_positive_times_are_nan_and_undetected(gamma_fit):
    traj = simulate_shedding(gamma_fit, n_individuals=5, times=[-1.0, 0.0, 5.0], seed=0)
    non_positive = traj[traj["time"] <= 0]
    assert non_positive["log10_value"].isna().all()
    assert not non_positive["detected"].any()


def test_scalar_incubation_shifts_the_curve(exponential_fit):
    unshifted = simulate_shedding(
        exponential_fit, n_individuals=30, times=[0.0], seed=5
    )
    shifted = simulate_shedding(
        exponential_fit,
        n_individuals=30,
        times=[5.0],
        incubation_period=5.0,
        seed=5,
    )
    np.testing.assert_allclose(
        unshifted["log10_value"].to_numpy(),
        shifted["log10_value"].to_numpy(),
    )


def test_array_incubation_must_match_n_individuals(exponential_fit):
    with pytest.raises(ValueError, match="incubation_period"):
        simulate_shedding(
            exponential_fit,
            n_individuals=5,
            times=[1.0],
            incubation_period=np.array([1.0, 2.0]),
        )


def test_callable_incubation_is_drawn_per_individual(exponential_fit):
    def incubation(rng, n):
        return rng.uniform(2.0, 8.0, size=n)

    traj = simulate_shedding(
        exponential_fit,
        n_individuals=40,
        times=[6.0],
        incubation_period=incubation,
        seed=4,
    )
    # Different offsets mean different effective times, hence spread in values.
    assert traj["log10_value"].std() > 0


def test_result_attrs_record_the_time_origin(exponential_fit):
    native = simulate_shedding(exponential_fit, n_individuals=3, times=[1.0], seed=0)
    assert native.attrs["time_origin"] == "symptom onset"
    assert native.attrs["incubation_applied"] is False

    infection = simulate_shedding(
        exponential_fit,
        n_individuals=3,
        times=[1.0],
        incubation_period=5.0,
        seed=0,
    )
    assert infection.attrs["time_origin"] == "infection"
    assert infection.attrs["incubation_applied"] is True
    assert infection.attrs["unit"] == "gc/mL"


def test_n_individuals_must_be_positive(exponential_fit):
    with pytest.raises(ValueError):
        simulate_shedding(exponential_fit, n_individuals=0, times=[1.0])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_shedding_simulate.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'shedding_hub.shedding_simulate'`

- [ ] **Step 3: Implement simulation**

Create `shedding_hub/shedding_simulate.py`:

```python
"""
Simulate shedding trajectories for synthetic infected individuals.

Intended for agent-based models: draw a cohort of individuals from a fitted
population distribution, then evaluate each one's shedding curve at whatever
times the simulation needs.
"""

from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd

from .shedding_models import log10_concentration_rowwise


def _resolve_incubation(
    incubation_period: Any, rng: np.random.Generator, n: int
) -> tuple[np.ndarray, bool]:
    """Return per-individual offsets in days, and whether a shift was applied."""
    if incubation_period is None:
        return np.zeros(n), False
    if callable(incubation_period):
        offsets = np.asarray(incubation_period(rng, n), dtype=float)
    elif np.isscalar(incubation_period):
        offsets = np.full(n, float(incubation_period))
    else:
        offsets = np.asarray(incubation_period, dtype=float)
    if offsets.shape != (n,):
        raise ValueError(
            "incubation_period must be None, a scalar, a callable, or an array of "
            f"length n_individuals ({n}); got shape {offsets.shape}."
        )
    return offsets, True


def simulate_shedding(
    source,
    *,
    n_individuals: int,
    times: Sequence[float] | np.ndarray,
    incubation_period: float | np.ndarray | Callable | None = None,
    include_measurement_error: bool = False,
    seed: int | None = None,
) -> pd.DataFrame:
    """
    Simulate shedding trajectories for synthetic individuals.

    Args:
        source: A ``SheddingFit`` or ``SheddingEnsemble``. Each individual's
            parameters are drawn from its population distribution.
        n_individuals: Number of individuals to simulate.
        times: Times at which to evaluate each trajectory. Measured from the
            fit's reference event, or from infection when ``incubation_period``
            is given.
        incubation_period: Days from infection to the reference event. ``None``
            leaves output in reference-event time. A scalar applies to everyone;
            an array of length ``n_individuals`` or a callable
            ``(rng, n) -> array`` gives each individual its own, which adds
            realistic timing variability across the cohort.
        include_measurement_error: Add ``N(0, sigma)`` assay noise on the log10
            scale. Off by default: an agent-based model wants the true shed
            concentration, and assay noise is a property of sampling rather than
            of the host.
        seed: Seed for a ``numpy`` generator, making runs reproducible.

    Returns:
        A tidy DataFrame with columns ``individual_id``, ``time``,
        ``log10_value``, ``value``, ``detected``, ``source_dataset_id``.
        ``detected`` is whether the value reaches the censoring limit. Values
        below the limit are reported as-is rather than clipped, so downstream
        mass-balance calculations stay correct. Under the gamma model,
        non-positive times yield NaN.

        ``result.attrs`` records ``time_origin``, ``incubation_applied``,
        ``model``, ``unit``, ``biomarker``, and ``specimen``.

    Note:
        The exponential model is defined for negative times but grows without
        bound going backwards, so simulating before the reference event
        extrapolates into implausible concentrations. Prefer ``times >= 0``.
    """
    if n_individuals < 1:
        raise ValueError("n_individuals must be at least 1")

    rng = np.random.default_rng(seed)
    params, sources = source.sample_params(rng, n_individuals)
    times = np.asarray(times, dtype=float)
    offsets, incubation_applied = _resolve_incubation(
        incubation_period, rng, n_individuals
    )

    shifted = times[None, :] - offsets[:, None]
    log10_values = log10_concentration_rowwise(source.model, params, shifted)

    if include_measurement_error:
        log10_values = log10_values + rng.normal(
            0.0, source.sigma, size=log10_values.shape
        )

    n_times = times.size
    frame = pd.DataFrame(
        {
            "individual_id": np.repeat(np.arange(n_individuals), n_times),
            "time": np.tile(times, n_individuals),
            "log10_value": log10_values.ravel(),
            "source_dataset_id": np.repeat(sources, n_times),
        }
    )
    frame["value"] = np.power(10.0, frame["log10_value"])
    detected = frame["log10_value"] >= source.censoring_limit
    frame["detected"] = detected.fillna(False).astype(bool)
    frame = frame[
        [
            "individual_id",
            "time",
            "log10_value",
            "value",
            "detected",
            "source_dataset_id",
        ]
    ]
    frame.attrs = {
        "time_origin": "infection" if incubation_applied else source.reference_event,
        "incubation_applied": incubation_applied,
        "model": source.model,
        "unit": source.unit,
        "biomarker": getattr(source, "biomarker", None),
        "specimen": getattr(source, "specimen", None),
    }
    return frame
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_shedding_simulate.py -v`
Expected: PASS, 11 tests.

- [ ] **Step 5: Format and commit**

```bash
black .
git add shedding_hub/shedding_simulate.py tests/test_shedding_simulate.py
git commit -m "feat: simulate shedding trajectories from a fitted model"
```

---

### Task 6: Catalog of fits with a browsable table

**Files:**
- Create: `shedding_hub/shedding_catalog.py`
- Test: `tests/test_shedding_catalog.py`

**Interfaces:**
- Consumes: `SheddingFit`, `fit_shedding_model`, `SheddingDataError`, `MODELS`.
- Produces:
  - `fit_to_row(fit: SheddingFit) -> dict` — the shared row builder used by both the catalog table and ensemble components
  - `SheddingCatalog` dataclass: fields `fits: list[SheddingFit]`, `skipped: pd.DataFrame`; property `table -> pd.DataFrame`; methods `select(**keys) -> SheddingFit`, `to_dict()`, `from_dict(payload)` (classmethod)
  - `fit_shedding_models(datasets, *, models=MODELS, min_observations=None) -> SheddingCatalog`
  - `load_shedding_catalog(path: str | None = None) -> SheddingCatalog`

Serialization omits `subject_params` (the catalog only needs `mu`, `Sigma`, `sigma` to simulate); a loaded fit has `subject_params is None`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_shedding_catalog.py`:

```python
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from shedding_hub.shedding_catalog import (
    SheddingCatalog,
    fit_shedding_models,
    load_shedding_catalog,
)


@pytest.fixture
def two_study_catalog(make_synthetic_dataset):
    mu = np.array([np.log(0.6), np.log(18.0)])
    cov = np.diag([0.04, 0.04])
    a = make_synthetic_dataset(
        "exponential", mu, cov, n_subjects=20, seed=1, dataset_id="study_a"
    )
    b = make_synthetic_dataset(
        "exponential", mu, cov, n_subjects=20, seed=2, dataset_id="study_b"
    )
    return fit_shedding_models([a, b], models=("exponential",))


def test_table_has_one_row_per_fit(two_study_catalog):
    table = two_study_catalog.table
    assert len(table) == 2
    assert set(table["dataset_id"]) == {"study_a", "study_b"}


def test_table_reports_medians_and_derived_quantities(two_study_catalog):
    table = two_study_catalog.table
    for column in [
        "dataset_id",
        "analyte",
        "biomarker",
        "specimen",
        "reference_event",
        "unit",
        "model",
        "n_subjects",
        "n_measurements",
        "pct_censored",
        "a_median",
        "sigma",
        "peak_day",
        "peak_log10",
        "half_life_days",
        "aic",
        "converged",
    ]:
        assert column in table.columns
    fit = two_study_catalog.fits[0]
    row = table.iloc[0]
    assert row["a_median"] == pytest.approx(np.exp(fit.population_mean[0]))
    assert row["half_life_days"] == pytest.approx(np.log(2.0) / row["a_median"])


def test_gamma_table_has_b_median_and_peak_day(make_synthetic_dataset):
    mu = np.array([np.log(0.5), np.log(2.0), np.log(12.0)])
    dataset = make_synthetic_dataset("gamma", mu, np.diag([0.04] * 3), n_subjects=20)
    catalog = fit_shedding_models([dataset], models=("gamma",))
    row = catalog.table.iloc[0]
    assert row["peak_day"] == pytest.approx(row["b_median"] / row["a_median"])


def test_select_returns_one_fit(two_study_catalog):
    fit = two_study_catalog.select(dataset_id="study_a")
    assert fit.dataset_id == "study_a"


def test_select_raises_on_ambiguous_match(two_study_catalog):
    with pytest.raises(ValueError, match="matched 2"):
        two_study_catalog.select(analyte="stool")


def test_select_raises_on_no_match(two_study_catalog):
    with pytest.raises(ValueError, match="matched no"):
        two_study_catalog.select(dataset_id="study_z")


def test_ct_analyte_is_recorded_in_skipped():
    dataset = {
        "dataset_id": "ct_study",
        "analytes": {
            "swab": {
                "specimen": "saliva",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "cycle threshold",
                "limit_of_quantification": "unknown",
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {
                "measurements": [
                    {"analyte": "swab", "time": 1, "value": 20.0},
                    {"analyte": "swab", "time": 2, "value": 25.0},
                ]
            }
        ],
    }
    catalog = fit_shedding_models([dataset], models=("exponential",))
    assert catalog.table.empty
    assert (catalog.skipped["reason"] == "ct_units").all()
    assert set(catalog.skipped["dataset_id"]) == {"ct_study"}


def test_cross_sectional_study_is_skipped():
    dataset = {
        "dataset_id": "cross_sectional",
        "analytes": {
            "stool": {
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "reference_event": "symptom onset",
                "unit": "gc/mL",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
            }
        },
        "participants": [
            {"measurements": [{"analyte": "stool", "time": 0, "value": 1e5}]}
            for _ in range(30)
        ],
    }
    catalog = fit_shedding_models([dataset], models=("exponential",))
    assert catalog.table.empty
    assert (catalog.skipped["reason"] == "too_few_subjects").all()


def test_round_trip_serialization(two_study_catalog, tmp_path):
    payload = two_study_catalog.to_dict()
    restored = SheddingCatalog.from_dict(payload)
    assert len(restored.fits) == len(two_study_catalog.fits)
    original = two_study_catalog.fits[0]
    copy = restored.select(
        dataset_id=original.dataset_id, model=original.model
    )
    np.testing.assert_allclose(copy.population_mean, original.population_mean)
    np.testing.assert_allclose(copy.population_cov, original.population_cov)
    assert copy.sigma == pytest.approx(original.sigma)
    assert copy.subject_params is None


def test_restored_fit_can_still_simulate(two_study_catalog):
    from shedding_hub.shedding_simulate import simulate_shedding

    restored = SheddingCatalog.from_dict(two_study_catalog.to_dict())
    traj = simulate_shedding(
        restored.fits[0], n_individuals=5, times=[1.0, 2.0], seed=0
    )
    assert len(traj) == 10


def test_shipped_catalog_covers_every_dataset():
    """CI staleness check: adding a dataset without regenerating must fail."""
    import pathlib

    data_dir = pathlib.Path(__file__).parent.parent / "data"
    on_disk = {
        path.name
        for path in data_dir.iterdir()
        if path.is_dir() and not path.name.startswith(".")
    }
    catalog = load_shedding_catalog()
    accounted = set(catalog.table["dataset_id"]) | set(
        catalog.skipped["dataset_id"]
    )
    missing = on_disk - accounted
    assert not missing, (
        f"Datasets absent from the shipped catalog: {sorted(missing)}. "
        "Run `make catalog` to regenerate."
    )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_shedding_catalog.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'shedding_hub.shedding_catalog'`

The final test (`test_shipped_catalog_covers_every_dataset`) will keep failing until Task 8 generates the catalog file. That is expected and is the point of the check.

- [ ] **Step 3: Implement the catalog**

Create `shedding_hub/shedding_catalog.py`:

```python
"""
A browsable catalog of fitted shedding models, and cross-study ensembles.

The catalog is the surface a modeller browses: one row per (dataset, analyte,
model), summarising each fit by its median individual, so a row reads as "peaks
day 4.2 at 6.8 log10 gc/mL, declines with a 1.5 day half-life" rather than as raw
log-parameters.
"""

import pathlib
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import yaml

from .shedding_fit import SheddingDataError, SheddingFit, fit_shedding_model
from .shedding_models import MODELS, PARAM_NAMES

CATALOG_PATH = pathlib.Path(__file__).parent / "data" / "shedding_catalog.yaml"

_KEY_COLUMNS = (
    "dataset_id",
    "analyte",
    "biomarker",
    "specimen",
    "reference_event",
    "unit",
    "gene_target",
    "dose",
    "vaccine_type",
    "model",
)


def fit_to_row(fit: SheddingFit) -> dict:
    """
    Summarize a fit as one table row describing its median individual.

    Because ``theta = log(params)`` is normal, the parameters are lognormal and
    ``exp(mu)`` is exactly their median. These are therefore labelled ``_median``,
    which is the accurate name rather than a compromise.
    """
    row = {column: getattr(fit, column) for column in _KEY_COLUMNS}
    medians = fit.median_params
    for name, value in zip(PARAM_NAMES[fit.model], medians):
        row[f"{name[0]}_median"] = float(value)
    row.update(
        {
            "n_subjects": fit.n_subjects,
            "n_measurements": fit.n_measurements,
            "pct_censored": (
                100.0 * fit.n_censored / fit.n_measurements
                if fit.n_measurements
                else np.nan
            ),
            "sigma": fit.sigma,
            "peak_day": fit.peak_day,
            "peak_log10": fit.peak_log10,
            "half_life_days": fit.half_life_days,
            "aic": fit.aic,
            "converged": fit.converged,
        }
    )
    return row


def _fits_to_frame(fits: list[SheddingFit]) -> pd.DataFrame:
    if not fits:
        return pd.DataFrame(columns=list(_KEY_COLUMNS))
    return pd.DataFrame([fit_to_row(fit) for fit in fits])


def _fit_to_payload(fit: SheddingFit) -> dict:
    """Serialize a fit, omitting per-subject parameters to keep the file small."""
    payload = {column: getattr(fit, column) for column in _KEY_COLUMNS}
    payload.update(
        {
            "method": fit.method,
            "population_mean": [float(v) for v in fit.population_mean],
            "population_cov": [
                [float(v) for v in row] for row in fit.population_cov
            ],
            "sigma": float(fit.sigma),
            "censoring_limit": float(fit.censoring_limit),
            "n_subjects": int(fit.n_subjects),
            "n_measurements": int(fit.n_measurements),
            "n_censored": int(fit.n_censored),
            "n_excluded_subjects": int(fit.n_excluded_subjects),
            "n_dropped_measurements": int(fit.n_dropped_measurements),
            "converged": bool(fit.converged),
            "log_likelihood": float(fit.log_likelihood),
            "aic": float(fit.aic),
        }
    )
    return payload


def _fit_from_payload(payload: dict) -> SheddingFit:
    return SheddingFit(
        model=payload["model"],
        method=payload.get("method", "mle"),
        population_mean=np.asarray(payload["population_mean"], dtype=float),
        population_cov=np.asarray(payload["population_cov"], dtype=float),
        sigma=float(payload["sigma"]),
        subject_params=None,
        censoring_limit=float(payload["censoring_limit"]),
        dataset_id=payload["dataset_id"],
        analyte=payload["analyte"],
        biomarker=payload.get("biomarker"),
        specimen=payload.get("specimen"),
        reference_event=payload.get("reference_event"),
        unit=payload.get("unit"),
        gene_target=payload.get("gene_target"),
        dose=payload.get("dose"),
        vaccine_type=payload.get("vaccine_type"),
        n_subjects=int(payload["n_subjects"]),
        n_measurements=int(payload["n_measurements"]),
        n_censored=int(payload["n_censored"]),
        n_excluded_subjects=int(payload["n_excluded_subjects"]),
        n_dropped_measurements=int(payload["n_dropped_measurements"]),
        converged=bool(payload["converged"]),
        log_likelihood=float(payload["log_likelihood"]),
        aic=float(payload["aic"]),
    )


@dataclass
class SheddingCatalog:
    """A collection of fitted models with a browsable summary table."""

    fits: list[SheddingFit] = field(default_factory=list)
    skipped: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(
            columns=["dataset_id", "analyte", "model", "reason", "message"]
        )
    )

    @property
    def table(self) -> pd.DataFrame:
        """One row per fit, summarising its median individual."""
        return _fits_to_frame(self.fits)

    def select(self, **keys) -> SheddingFit:
        """
        Return the single fit matching ``keys``.

        Raises:
            ValueError: If no fit matches, or if more than one does. Never picks
                silently — the error lists the candidates and the columns that
                would tell them apart.
        """
        matches = [
            fit
            for fit in self.fits
            if all(getattr(fit, key, None) == value for key, value in keys.items())
        ]
        if not matches:
            raise ValueError(
                f"select({keys}) matched no fits. "
                f"Available combinations are listed in `catalog.table`; "
                f"{len(self.fits)} fit(s) are loaded."
            )
        if len(matches) > 1:
            candidates = _fits_to_frame(matches)
            distinguishing = [
                column
                for column in _KEY_COLUMNS
                if column not in keys and candidates[column].nunique() > 1
            ]
            raise ValueError(
                f"select({keys}) matched {len(matches)} fits. Narrow it with "
                f"{distinguishing or 'more specific keys'}. Candidates:\n"
                f"{candidates[list(_KEY_COLUMNS)].to_string(index=False)}"
            )
        return matches[0]

    def ensemble(self, *, dataset_ids=None, weights="n_subjects", method="mixture", **keys):
        """Build a ``SheddingEnsemble`` from the matching fits. See Task 7."""
        from .shedding_ensemble import build_ensemble

        return build_ensemble(
            self, dataset_ids=dataset_ids, weights=weights, method=method, **keys
        )

    def to_dict(self) -> dict:
        return {
            "fits": [_fit_to_payload(fit) for fit in self.fits],
            "skipped": self.skipped.to_dict(orient="records"),
        }

    @classmethod
    def from_dict(cls, payload: dict) -> "SheddingCatalog":
        skipped = pd.DataFrame(
            payload.get("skipped") or [],
            columns=["dataset_id", "analyte", "model", "reason", "message"],
        )
        return cls(
            fits=[_fit_from_payload(item) for item in payload.get("fits", [])],
            skipped=skipped,
        )


def fit_shedding_models(
    datasets,
    *,
    models=MODELS,
    min_observations: int | None = None,
) -> SheddingCatalog:
    """
    Fit every analyte of every dataset, for every requested model.

    Analytes that cannot be fitted are recorded in ``catalog.skipped`` with a
    reason rather than raising, so one unsuitable analyte does not abort a
    repository-wide build — and so a missing study reads as unsuitable rather
    than as a bug.

    Args:
        datasets: Dataset dictionaries from ``load_dataset``.
        models: Model names to fit. Defaults to both.
        min_observations: Passed through to the fitter.

    Returns:
        A ``SheddingCatalog``.
    """
    fits: list[SheddingFit] = []
    skipped: list[dict] = []

    for dataset in datasets:
        dataset_id = dataset.get("dataset_id", "unknown")
        for analyte in dataset.get("analytes", {}):
            for model in models:
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", UserWarning)
                        fits.append(
                            fit_shedding_model(
                                dataset,
                                analyte=analyte,
                                model=model,
                                min_observations=min_observations,
                            )
                        )
                except SheddingDataError as error:
                    skipped.append(
                        {
                            "dataset_id": dataset_id,
                            "analyte": analyte,
                            "model": model,
                            "reason": error.reason,
                            "message": str(error),
                        }
                    )
                except (ValueError, np.linalg.LinAlgError) as error:
                    skipped.append(
                        {
                            "dataset_id": dataset_id,
                            "analyte": analyte,
                            "model": model,
                            "reason": "did_not_converge",
                            "message": str(error),
                        }
                    )

    return SheddingCatalog(
        fits=fits,
        skipped=pd.DataFrame(
            skipped, columns=["dataset_id", "analyte", "model", "reason", "message"]
        ),
    )


def load_shedding_catalog(path: str | None = None) -> SheddingCatalog:
    """
    Load the catalog of precomputed estimates shipped with the package.

    Args:
        path: Optional path to a catalog YAML. Defaults to the shipped file.

    Returns:
        A ``SheddingCatalog``. Loaded fits carry ``subject_params is None``
        because per-subject values are not serialized; everything needed to
        simulate (``mu``, ``Sigma``, ``sigma``) is present.
    """
    catalog_path = pathlib.Path(path) if path else CATALOG_PATH
    if not catalog_path.is_file():
        raise FileNotFoundError(
            f"No shedding catalog at {catalog_path}. Run `make catalog` to build it."
        )
    with catalog_path.open(encoding="utf-8") as stream:
        payload = yaml.safe_load(stream)
    return SheddingCatalog.from_dict(payload)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_shedding_catalog.py -v -k "not shipped_catalog"`
Expected: PASS, 10 tests. `test_shipped_catalog_covers_every_dataset` still fails until Task 8; that is expected.

- [ ] **Step 5: Format and commit**

```bash
black .
git add shedding_hub/shedding_catalog.py tests/test_shedding_catalog.py
git commit -m "feat: add browsable catalog of fitted shedding models"
```

---

### Task 7: Cross-study ensembles

**Files:**
- Create: `shedding_hub/shedding_ensemble.py`
- Test: `tests/test_shedding_ensemble.py`

**Interfaces:**
- Consumes: `SheddingFit`, `SheddingCatalog`, `fit_to_row`.
- Produces:
  - `SheddingEnsemble` dataclass: fields `fits: list[SheddingFit]`, `weights: np.ndarray`, `method: str`; properties `model`, `sigma`, `censoring_limit`, `reference_event`, `unit`, `biomarker`, `specimen`, `components` (DataFrame), `median_params` (moment only); method `sample_params(rng, n)`
  - `make_ensemble(fits, *, weights="n_subjects", method="mixture") -> SheddingEnsemble`
  - `build_ensemble(catalog, *, dataset_ids=None, weights=..., method=..., **keys) -> SheddingEnsemble`

Note the dataclass field is `fits` while `components` is the DataFrame view, matching the spec's `ens.components`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_shedding_ensemble.py`:

```python
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from shedding_hub.shedding_catalog import fit_shedding_models
from shedding_hub.shedding_ensemble import SheddingEnsemble, make_ensemble


@pytest.fixture
def catalog(make_synthetic_dataset):
    cov = np.diag([0.04, 0.04])

    def study(dataset_id, mu, seed):
        return make_synthetic_dataset(
            "exponential",
            np.asarray(mu),
            cov,
            n_subjects=20,
            seed=seed,
            dataset_id=dataset_id,
        )

    return fit_shedding_models(
        [
            study("study_a", [np.log(0.6), np.log(18.0)], 1),
            study("study_b", [np.log(0.4), np.log(20.0)], 2),
            study("study_c", [np.log(0.8), np.log(16.0)], 3),
        ],
        models=("exponential",),
    )


def test_ensemble_over_all_matching_studies(catalog):
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2", specimen="stool")
    assert len(ensemble.fits) == 3
    assert len(ensemble.components) == 3
    assert set(ensemble.components["dataset_id"]) == {
        "study_a",
        "study_b",
        "study_c",
    }


def test_ensemble_restricted_to_named_studies(catalog):
    ensemble = catalog.ensemble(
        biomarker="SARS-CoV-2", dataset_ids=["study_a", "study_c"]
    )
    assert set(ensemble.components["dataset_id"]) == {"study_a", "study_c"}


def test_unmatched_dataset_id_raises(catalog):
    with pytest.raises(ValueError, match="study_z"):
        catalog.ensemble(biomarker="SARS-CoV-2", dataset_ids=["study_a", "study_z"])


def test_single_component_ensemble_matches_the_underlying_fit(catalog):
    fit = catalog.select(dataset_id="study_a")
    ensemble = make_ensemble([fit])
    a = fit.sample_params(np.random.default_rng(0), 50)[0]
    b = ensemble.sample_params(np.random.default_rng(0), 50)[0]
    np.testing.assert_allclose(a, b)


def test_mixture_draws_from_every_study_and_records_the_source(catalog):
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2", weights="equal")
    _, sources = ensemble.sample_params(np.random.default_rng(0), 600)
    assert set(sources.tolist()) == {"study_a", "study_b", "study_c"}
    counts = np.array(
        [(sources == name).sum() for name in ["study_a", "study_b", "study_c"]]
    )
    # Equal weights: each study should take roughly a third.
    assert (counts > 120).all()


def test_weights_default_to_subject_counts(catalog):
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2")
    expected = np.array([fit.n_subjects for fit in ensemble.fits], dtype=float)
    np.testing.assert_allclose(ensemble.weights, expected / expected.sum())


def test_moment_covariance_is_within_plus_between():
    """Hand-computable two-study example."""
    from shedding_hub.shedding_fit import SheddingFit
    import pandas as pd

    def stub(dataset_id, mean, cov, n_subjects):
        return SheddingFit(
            model="exponential",
            method="mle",
            population_mean=np.asarray(mean, float),
            population_cov=np.asarray(cov, float),
            sigma=0.3,
            subject_params=pd.DataFrame(),
            censoring_limit=2.0,
            dataset_id=dataset_id,
            analyte="stool",
            biomarker="SARS-CoV-2",
            specimen="stool",
            reference_event="symptom onset",
            unit="gc/mL",
            gene_target=None,
            dose=None,
            vaccine_type=None,
            n_subjects=n_subjects,
            n_measurements=100,
            n_censored=10,
            n_excluded_subjects=0,
            n_dropped_measurements=0,
            converged=True,
            log_likelihood=-1.0,
            aic=2.0,
        )

    a = stub("a", [0.0, 0.0], np.eye(2), 10)
    b = stub("b", [2.0, 0.0], np.eye(2), 10)
    ensemble = make_ensemble([a, b], weights="equal", method="moment")

    np.testing.assert_allclose(ensemble.population_mean, [1.0, 0.0])
    # within = I; between = weighted cov of means = [[1, 0], [0, 0]]
    expected = np.eye(2) + np.array([[1.0, 0.0], [0.0, 0.0]])
    np.testing.assert_allclose(ensemble.population_cov, expected)
    np.testing.assert_allclose(ensemble.median_params, np.exp([1.0, 0.0]))


def test_mixed_units_raise(catalog):
    fits = list(catalog.fits)
    fits[1].unit = "gc/dry gram"
    with pytest.raises(ValueError, match="unit"):
        make_ensemble(fits)


def test_mixed_reference_events_raise(catalog):
    fits = list(catalog.fits)
    fits[1].reference_event = "enrollment"
    with pytest.raises(ValueError, match="reference_event"):
        make_ensemble(fits)


def test_two_analytes_from_one_study_raise(catalog):
    fits = list(catalog.fits)
    fits[1].dataset_id = "study_a"
    fits[1].analyte = "stool_orf1a"
    with pytest.raises(ValueError, match="more than one"):
        make_ensemble(fits)


def test_ensemble_can_be_simulated(catalog):
    from shedding_hub.shedding_simulate import simulate_shedding

    ensemble = catalog.ensemble(biomarker="SARS-CoV-2")
    traj = simulate_shedding(ensemble, n_individuals=30, times=[1.0, 5.0], seed=0)
    assert len(traj) == 60
    assert traj["source_dataset_id"].nunique() > 1


def test_mixture_has_no_median_params(catalog):
    ensemble = catalog.ensemble(biomarker="SARS-CoV-2", method="mixture")
    with pytest.raises(ValueError, match="mixture"):
        _ = ensemble.median_params
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_shedding_ensemble.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'shedding_hub.shedding_ensemble'`

- [ ] **Step 3: Implement the ensemble**

Create `shedding_hub/shedding_ensemble.py`:

```python
"""
Combine per-study shedding fits into a cross-study ensemble.

Two methods are available. ``mixture`` draws a study per simulated individual and
then draws that individual's parameters from the study's fit, preserving
between-study heterogeneity and keeping the distribution multimodal when studies
genuinely disagree. ``moment`` collapses the components into a single Gaussian
whose covariance is within-study plus between-study variance.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .shedding_catalog import SheddingCatalog, _fits_to_frame
from .shedding_fit import SheddingFit

_COMPATIBILITY_KEYS = (
    "model",
    "unit",
    "reference_event",
    "biomarker",
    "specimen",
)


@dataclass
class SheddingEnsemble:
    """An ensemble of per-study fits sharing a biomarker, specimen, and unit."""

    fits: list[SheddingFit]
    weights: np.ndarray
    method: str

    @property
    def components(self) -> pd.DataFrame:
        """One row per contributing fit, with the same columns as the catalog."""
        frame = _fits_to_frame(self.fits)
        frame.insert(0, "weight", self.weights)
        return frame

    @property
    def model(self) -> str:
        return self.fits[0].model

    @property
    def sigma(self) -> float:
        """Weighted mean measurement-error SD across components."""
        return float(np.sum(self.weights * np.array([f.sigma for f in self.fits])))

    @property
    def censoring_limit(self) -> float:
        """The most conservative (highest) limit among the components."""
        return float(max(fit.censoring_limit for fit in self.fits))

    @property
    def reference_event(self) -> str | None:
        return self.fits[0].reference_event

    @property
    def unit(self) -> str | None:
        return self.fits[0].unit

    @property
    def biomarker(self) -> str | None:
        return self.fits[0].biomarker

    @property
    def specimen(self) -> str | None:
        return self.fits[0].specimen

    @property
    def population_mean(self) -> np.ndarray:
        """Weighted mean of the component means."""
        means = np.array([fit.population_mean for fit in self.fits])
        return np.sum(self.weights[:, None] * means, axis=0)

    @property
    def population_cov(self) -> np.ndarray:
        """
        Moment-matched covariance: within-study plus between-study.

        Defined for both methods, but only used for sampling under
        ``method="moment"``.
        """
        means = np.array([fit.population_mean for fit in self.fits])
        covs = np.array([fit.population_cov for fit in self.fits])
        within = np.sum(self.weights[:, None, None] * covs, axis=0)
        deviations = means - self.population_mean
        between = np.einsum(
            "s,si,sj->ij", self.weights, deviations, deviations
        )
        return within + between

    @property
    def median_params(self) -> np.ndarray:
        """
        Parameters of the median individual.

        Only well defined for ``method="moment"``. A mixture has no closed-form
        median, so rather than report a misleading number, simulate with
        ``simulate_shedding`` and take empirical quantiles of the result — which
        is the quantity you actually want.
        """
        if self.method != "moment":
            raise ValueError(
                "median_params is only defined for method='moment'. This is a "
                "mixture ensemble, which has no closed-form median: simulate and "
                "take empirical quantiles instead."
            )
        return np.exp(self.population_mean)

    def sample_params(
        self, rng: np.random.Generator, n: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """Draw ``n`` individuals' natural-scale parameters."""
        if n < 1:
            raise ValueError("n_individuals must be at least 1")
        if self.method == "moment":
            theta = rng.multivariate_normal(
                self.population_mean, self.population_cov, n
            )
            return np.exp(theta), np.full(n, "ensemble", dtype=object)

        if len(self.fits) == 1:
            # Skip the categorical draw entirely so a one-study ensemble consumes
            # the generator exactly as the underlying fit would, making the two
            # interchangeable for a given seed.
            return self.fits[0].sample_params(rng, n)

        choices = rng.choice(len(self.fits), size=n, p=self.weights)
        k = self.fits[0].population_mean.size
        theta = np.empty((n, k))
        sources = np.empty(n, dtype=object)
        for index, fit in enumerate(self.fits):
            mask = choices == index
            count = int(mask.sum())
            if not count:
                continue
            theta[mask] = rng.multivariate_normal(
                fit.population_mean, fit.population_cov, count
            )
            sources[mask] = fit.dataset_id
        return np.exp(theta), sources


def _resolve_weights(fits: list[SheddingFit], weights) -> np.ndarray:
    if weights == "equal":
        raw = np.ones(len(fits))
    elif weights == "n_subjects":
        raw = np.array([fit.n_subjects for fit in fits], dtype=float)
    else:
        raw = np.asarray(weights, dtype=float)
        if raw.shape != (len(fits),):
            raise ValueError(
                f"weights must be 'n_subjects', 'equal', or an array of length "
                f"{len(fits)}; got shape {raw.shape}."
            )
    if raw.sum() <= 0:
        raise ValueError("Ensemble weights must sum to a positive value.")
    return raw / raw.sum()


def make_ensemble(
    fits: list[SheddingFit],
    *,
    weights="n_subjects",
    method: str = "mixture",
) -> SheddingEnsemble:
    """
    Assemble an ensemble from explicit fits.

    Accepts fits from anywhere — the shipped catalog, a fresh fit on private
    data, or a mixture of both.

    Args:
        fits: Component fits. A single-component ensemble is legal and behaves
            identically to the underlying fit, so callers can keep one code path
            regardless of how many studies they selected.
        weights: ``"n_subjects"`` (default), ``"equal"``, or an explicit array.
        method: ``"mixture"`` (default) or ``"moment"``.

    Returns:
        A ``SheddingEnsemble``.

    Raises:
        ValueError: If the fits disagree on model, unit, reference event,
            biomarker, or specimen, or if one study contributes more than one
            analyte.
    """
    fits = list(fits)
    if not fits:
        raise ValueError("An ensemble needs at least one fit.")
    if method not in ("mixture", "moment"):
        raise ValueError(
            f"Unknown ensemble method {method!r}; choose 'mixture' or 'moment'."
        )

    for key in _COMPATIBILITY_KEYS:
        values = {getattr(fit, key) for fit in fits}
        if len(values) > 1:
            raise ValueError(
                f"Ensemble components disagree on {key}: {sorted(map(str, values))}. "
                "Estimates are only comparable within one of these."
            )

    dataset_ids = [fit.dataset_id for fit in fits]
    duplicated = {name for name in dataset_ids if dataset_ids.count(name) > 1}
    if duplicated:
        offenders = _fits_to_frame(
            [fit for fit in fits if fit.dataset_id in duplicated]
        )
        raise ValueError(
            f"Study/studies {sorted(duplicated)} contribute more than one analyte "
            "to this ensemble, which would enter their subjects twice. Narrow the "
            "selection (for example by gene_target or analyte). Candidates:\n"
            f"{offenders[['dataset_id', 'analyte', 'gene_target']].to_string(index=False)}"
        )

    return SheddingEnsemble(
        fits=fits, weights=_resolve_weights(fits, weights), method=method
    )


def build_ensemble(
    catalog: SheddingCatalog,
    *,
    dataset_ids=None,
    weights="n_subjects",
    method: str = "mixture",
    **keys,
) -> SheddingEnsemble:
    """
    Build an ensemble from the catalog fits matching ``keys``.

    Args:
        catalog: The catalog to draw components from.
        dataset_ids: Optional list restricting components to named studies. A
            name with no matching fit raises rather than being dropped, which
            would silently shrink the ensemble.
        weights: Passed to ``make_ensemble``.
        method: Passed to ``make_ensemble``.
        **keys: Attribute filters, e.g. ``biomarker="SARS-CoV-2"``.

    Returns:
        A ``SheddingEnsemble``.
    """
    matches = [
        fit
        for fit in catalog.fits
        if all(getattr(fit, key, None) == value for key, value in keys.items())
    ]
    if dataset_ids is not None:
        requested = list(dataset_ids)
        available = {fit.dataset_id for fit in matches}
        missing = [name for name in requested if name not in available]
        if missing:
            raise ValueError(
                f"No fit matching {keys} for dataset_id(s) {missing}. "
                f"Matching studies are {sorted(available)}."
            )
        matches = [fit for fit in matches if fit.dataset_id in set(requested)]
    if not matches:
        raise ValueError(
            f"No fits match {keys}. Browse `catalog.table` for available "
            "combinations."
        )
    return make_ensemble(matches, weights=weights, method=method)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_shedding_ensemble.py -v`
Expected: PASS, 12 tests.

- [ ] **Step 5: Format and commit**

```bash
black .
git add shedding_hub/shedding_ensemble.py tests/test_shedding_ensemble.py
git commit -m "feat: add cross-study shedding ensembles"
```

---

### Task 8: Build and ship the catalog

**Files:**
- Create: `scripts/build_shedding_catalog.py`
- Create: `shedding_hub/data/shedding_catalog.yaml` (generated)
- Modify: `Makefile`
- Modify: `pyproject.toml`

**Interfaces:**
- Consumes: `fit_shedding_models`, `load_dataset`.
- Produces: the shipped catalog file, satisfying `test_shipped_catalog_covers_every_dataset` from Task 6.

- [ ] **Step 1: Write the build script**

Create `scripts/build_shedding_catalog.py`:

```python
"""
Regenerate the precomputed shedding-model catalog shipped with the package.

Run via `make catalog`. Fitting every analyte of every dataset takes a while,
which is exactly why the result is precomputed rather than fitted on demand.
"""

import argparse
import pathlib
import sys
import warnings

import yaml

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from shedding_hub import load_dataset  # noqa: E402
from shedding_hub.shedding_catalog import (  # noqa: E402
    CATALOG_PATH,
    fit_shedding_models,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data", default=str(REPO_ROOT / "data"), help="Directory of datasets."
    )
    parser.add_argument(
        "--output", default=str(CATALOG_PATH), help="Catalog file to write."
    )
    args = parser.parse_args()

    data_dir = pathlib.Path(args.data)
    dataset_ids = sorted(
        path.name
        for path in data_dir.iterdir()
        if path.is_dir() and not path.name.startswith(".")
    )

    datasets = []
    for dataset_id in dataset_ids:
        print(f"loading {dataset_id}", flush=True)
        datasets.append(load_dataset(dataset_id, local=str(data_dir)))

    print(f"fitting {len(datasets)} dataset(s)", flush=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        catalog = fit_shedding_models(datasets)

    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(catalog.to_dict(), stream, sort_keys=False)

    print(f"wrote {len(catalog.fits)} fit(s) to {output}")
    print(f"skipped {len(catalog.skipped)} analyte/model combination(s)")
    if not catalog.skipped.empty:
        print(catalog.skipped["reason"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Add the make target**

In `Makefile`, change the first line to include `catalog`:

```make
.PHONY : backup_data assert_data_unchanged extraction catalog
```

and append at the end of the file:

```make
# Refit every analyte in data/ and rewrite the shipped catalog. Slow by design;
# run it whenever datasets are added or changed.
catalog :
	python scripts/build_shedding_catalog.py
```

- [ ] **Step 3: Declare the catalog as package data**

Append to `pyproject.toml`:

```toml
[tool.setuptools.package-data]
shedding_hub = ["data/*.yaml"]
```

- [ ] **Step 4: Build the catalog**

Run: `python scripts/build_shedding_catalog.py`
Expected: progress lines per dataset, then a summary such as `wrote N fit(s) to .../shedding_hub/data/shedding_catalog.yaml` and a table of skip reasons. This takes several minutes.

- [ ] **Step 5: Verify the coverage check now passes**

Run: `pytest tests/test_shedding_catalog.py -v`
Expected: PASS, all 11 tests including `test_shipped_catalog_covers_every_dataset`.

- [ ] **Step 6: Sanity-check the catalog by hand**

Run:

```bash
python -c "
import shedding_hub as sh
cat = sh.shedding_catalog.load_shedding_catalog()
t = cat.table
print(t.shape)
print(t[t.dataset_id=='woelfel2020virological'][['analyte','model','n_subjects','peak_day','peak_log10','half_life_days']].to_string(index=False))
print(cat.skipped['reason'].value_counts().to_string())
"
```

Expected: woelfel stool/sputum rows with plausible values — peak within roughly the first two weeks and `peak_log10` in the 4–9 range for `gc/mL`. If `peak_day` is negative or `peak_log10` exceeds ~15, stop and investigate the fit rather than shipping the catalog.

- [ ] **Step 7: Commit**

```bash
black .
git add scripts/build_shedding_catalog.py shedding_hub/data/shedding_catalog.yaml Makefile pyproject.toml
git commit -m "feat: build and ship precomputed shedding catalog"
```

---

### Task 9: Plotting, package exports, and README

**Files:**
- Modify: `shedding_hub/shedding_simulate.py` (append)
- Modify: `shedding_hub/__init__.py`
- Modify: `README.md`
- Test: `tests/test_shedding_simulate.py` (append)

**Interfaces:**
- Produces: `plot_simulated_shedding(traj, *, source=None, observed=None, quantiles=(0.05, 0.5, 0.95), figsize=(8, 6)) -> matplotlib.figure.Figure`, and the public exports.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_shedding_simulate.py`:

```python
from matplotlib.figure import Figure

from shedding_hub.shedding_simulate import plot_simulated_shedding


def test_plot_returns_a_figure(exponential_fit):
    traj = simulate_shedding(
        exponential_fit, n_individuals=40, times=np.arange(0.0, 20.0), seed=0
    )
    fig = plot_simulated_shedding(traj, source=exponential_fit)
    assert isinstance(fig, Figure)


def test_plot_rejects_an_empty_frame():
    with pytest.raises(ValueError, match="empty"):
        plot_simulated_shedding(pd.DataFrame())


def test_public_exports_are_available():
    import shedding_hub as sh

    for name in [
        "fit_shedding_model",
        "fit_shedding_models",
        "load_shedding_catalog",
        "make_ensemble",
        "simulate_shedding",
        "plot_simulated_shedding",
        "SheddingFit",
        "SheddingEnsemble",
        "SheddingCatalog",
    ]:
        assert hasattr(sh, name), name
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_shedding_simulate.py -v -k "plot or exports"`
Expected: FAIL — `ImportError: cannot import name 'plot_simulated_shedding'`

- [ ] **Step 3: Implement the plot**

Append to `shedding_hub/shedding_simulate.py`:

```python
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def plot_simulated_shedding(
    traj: pd.DataFrame,
    *,
    source=None,
    observed: dict | None = None,
    quantiles: tuple[float, float, float] = (0.05, 0.5, 0.95),
    figsize: tuple[float, float] = (8, 6),
) -> Figure:
    """
    Plot the median and a credible band of simulated trajectories.

    Args:
        traj: Output of ``simulate_shedding``.
        source: Optional fit or ensemble, used to draw the censoring limit.
        observed: Optional dataset dictionary; its measurements are overlaid as
            points so simulated and real trajectories can be compared.
        quantiles: Lower, middle, and upper quantiles for the band.
        figsize: Figure size in inches.

    Returns:
        The figure. It is closed in the pyplot state so notebooks do not display
        it twice, matching the convention in ``shedding_peak.py``.
    """
    if traj.empty:
        raise ValueError("Simulation result is empty, cannot create plot")

    lower, middle, upper = quantiles
    summary = (
        traj.dropna(subset=["log10_value"])
        .groupby("time")["log10_value"]
        .quantile([lower, middle, upper])
        .unstack()
    )

    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(
        summary.index,
        summary[lower],
        summary[upper],
        alpha=0.25,
        color="tab:blue",
        label=f"{int((upper - lower) * 100)}% of simulated individuals",
    )
    ax.plot(
        summary.index, summary[middle], color="tab:blue", lw=2, label="Median"
    )

    if source is not None:
        ax.axhline(
            source.censoring_limit,
            ls=":",
            color="gray",
            label="Limit of quantification",
        )

    if observed is not None:
        analyte = getattr(source, "analyte", None)
        times, values = [], []
        for participant in observed.get("participants", []):
            for measurement in participant.get("measurements") or []:
                if analyte is not None and measurement.get("analyte") != analyte:
                    continue
                time = measurement.get("time")
                value = measurement.get("value")
                if isinstance(time, (int, float)) and isinstance(
                    value, (int, float)
                ):
                    times.append(float(time))
                    values.append(np.log10(float(value)))
        if times:
            ax.scatter(
                times, values, s=18, color="black", alpha=0.5, label="Observed"
            )

    origin = traj.attrs.get("time_origin", "reference event")
    unit = traj.attrs.get("unit", "")
    ax.set_xlabel(f"Days after {origin}")
    ax.set_ylabel(f"log10 concentration ({unit})" if unit else "log10 concentration")
    ax.set_title("Simulated shedding trajectories")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.close(fig)
    return fig
```

- [ ] **Step 4: Add public exports**

Replace `shedding_hub/__init__.py` imports and `__all__` — append these imports after the existing ones:

```python
from .shedding_models import MODELS, PARAM_NAMES

from .shedding_fit import SheddingDataError, SheddingFit, fit_shedding_model

from .shedding_catalog import (
    SheddingCatalog,
    fit_shedding_models,
    load_shedding_catalog,
)

from .shedding_ensemble import SheddingEnsemble, make_ensemble

from .shedding_simulate import plot_simulated_shedding, simulate_shedding
```

and add these entries to the end of the `__all__` list:

```python
    "MODELS",
    "PARAM_NAMES",
    "SheddingDataError",
    "SheddingFit",
    "SheddingCatalog",
    "SheddingEnsemble",
    "fit_shedding_model",
    "fit_shedding_models",
    "load_shedding_catalog",
    "make_ensemble",
    "simulate_shedding",
    "plot_simulated_shedding",
```

- [ ] **Step 5: Run the tests**

Run: `pytest tests/test_shedding_simulate.py -v`
Expected: PASS, 14 tests.

- [ ] **Step 6: Add the README section**

In `README.md`, after the "Statistical Summaries" section and before "Visualization", insert:

````markdown
### Simulating Shedding

Simulate shedding trajectories for synthetic infected individuals — intended for
agent-based models of wastewater surveillance. Browse the catalog of fitted
estimates, pick one study or an ensemble across studies, then simulate.

```python
>>> import numpy as np
>>> import shedding_hub as sh
>>> catalog = sh.load_shedding_catalog()
>>> catalog.table[['dataset_id', 'specimen', 'model', 'peak_day']].head()  # doctest: +SKIP
>>> fit = catalog.select(
...     dataset_id='woelfel2020virological', analyte='stool', model='gamma'
... )  # doctest: +SKIP
>>> traj = sh.simulate_shedding(
...     fit, n_individuals=100, times=np.arange(0, 30), seed=42
... )  # doctest: +SKIP
>>> list(traj.columns)  # doctest: +SKIP
['individual_id', 'time', 'log10_value', 'value', 'detected', 'source_dataset_id']

```

Pass `incubation_period` to express times as days since infection rather than
days since the study's reference event:

```python
>>> traj = sh.simulate_shedding(
...     fit, n_individuals=100, times=np.arange(0, 30),
...     incubation_period=5.0, seed=42
... )  # doctest: +SKIP
>>> traj.attrs['time_origin']  # doctest: +SKIP
'infection'

```

To pool evidence across studies, build an ensemble. Each simulated individual is
drawn from one contributing study, so between-study variation is preserved:

```python
>>> ensemble = catalog.ensemble(
...     biomarker='SARS-CoV-2', specimen='stool',
...     reference_event='symptom onset', unit='gc/mL', model='gamma',
... )  # doctest: +SKIP
>>> traj = sh.simulate_shedding(
...     ensemble, n_individuals=1000, times=np.arange(0, 30), seed=42
... )  # doctest: +SKIP

```

Estimates come from a censored maximum-likelihood fit, so `negative`
measurements inform the fit rather than being discarded. Because the two-stage
fit does not shrink individual estimates toward the population mean, simulated
cohorts are somewhat more dispersed than reality.
````

- [ ] **Step 7: Verify README doctests pass, then check the example by hand**

The catalog-dependent lines carry `# doctest: +SKIP` because which rows exist
depends on which fits converged, and a README example that silently breaks CI
whenever a dataset changes is worse than one that is verified deliberately.

Run: `python -m doctest -o ELLIPSIS -o NORMALIZE_WHITESPACE README.md`
Expected: no output (no failures).

Then confirm the skipped example actually works against the built catalog:

```bash
python -c "
import numpy as np, shedding_hub as sh
catalog = sh.load_shedding_catalog()
fit = catalog.select(dataset_id='woelfel2020virological', analyte='stool', model='gamma')
traj = sh.simulate_shedding(fit, n_individuals=100, times=np.arange(0, 30), seed=42)
print(list(traj.columns))
print(traj.attrs['time_origin'])
"
```

Expected: the column list from the README, then `symptom onset`. If `select`
raises because that fit did not converge, pick a row that exists from
`catalog.table` and update the README example to match.

- [ ] **Step 8: Run the whole suite**

Run: `pytest -v && black --check .`
Expected: all tests pass and formatting is clean.

- [ ] **Step 9: Commit**

```bash
git add shedding_hub/shedding_simulate.py shedding_hub/__init__.py README.md tests/test_shedding_simulate.py
git commit -m "feat: add simulation plot, public exports, and README section"
```

---

## Self-Review Notes

Checked against `docs/superpowers/specs/2026-07-18-shedding-simulation-design.md`:

- **Spec coverage.** Models → Task 1. Data preparation, censoring, all exclusion rules → Task 2. Censored MLE, two-stage population, `SheddingFit` → Task 3. Tutorial validation → Task 4. Simulation, incubation shift, measurement error, `detected` → Task 5. Catalog, table, `select`, `skipped`, serialization → Task 6. Ensemble mixture/moment, three selection levels, `make_ensemble`, compatibility rules → Task 7. Shipped catalog, `make catalog`, package data, CI coverage check → Task 8. Plot, exports, README → Task 9.
- **Deviation from spec.** The spec named a single `shedding_hub/simulate.py` and `tests/test_simulate.py` while permitting a split; this plan splits into four modules and five test files, since the feature spans four distinct responsibilities. The spec's `Files` section anticipated this.
- **Known open item.** `aic` is not comparable between models when the gamma fit dropped non-positive times that the exponential fit kept. Surfaced via `n_measurements` and documented in `fit_shedding_model`; not silently hidden.
- **Type consistency.** `sample_params(rng, n) -> (params, sources)` is implemented identically on `SheddingFit` (Task 3) and `SheddingEnsemble` (Task 7), which is what lets `simulate_shedding` (Task 5) accept either. `fit_to_row` is defined once in Task 6 and reused by `SheddingEnsemble.components` in Task 7. The `SheddingEnsemble` dataclass field is `fits`; `components` is the DataFrame view, matching the spec's `ens.components`.
- **Circular import.** `shedding_catalog.ensemble()` imports `build_ensemble` from `shedding_ensemble` inside the method body, because `shedding_ensemble` imports from `shedding_catalog` at module level. Keep that import local.

Four defects found during self-review and fixed inline:

1. **Cross-module test imports would not resolve.** Tasks 5–7 originally did `from tests.test_shedding_fit import _synthetic_dataset`, but `tests/` has no `__init__.py`, so `tests` is not an importable package. Replaced with a `make_synthetic_dataset` factory fixture in `tests/conftest.py` (Task 3, Step 1), which pytest makes available to every test module without any import.
2. **A tautological assertion.** `assert fit.n_measurements == len(fit.subject_params) * 0 + fit.n_measurements` asserted nothing. Replaced with the actual expected count (`140` = 10 subjects × 14 times) plus a real bound on `n_censored`.
3. **Single-component ensembles would not match their underlying fit.** `sample_params` called `rng.choice` before `rng.multivariate_normal`, consuming generator state the bare fit does not, so identical seeds gave different draws and the "one code path" property the spec promises would have been false. Added a short-circuit that delegates straight to the sole fit.
4. **Fragile README doctests.** README doctests run in CI, and the example depends on which fits converge — a dataset change could break the build for unrelated reasons. The catalog-dependent lines now carry `# doctest: +SKIP`, with Task 9 Step 7 verifying the example by hand against the built catalog instead.
