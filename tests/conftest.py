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
        theta = rng.multivariate_normal(
            np.asarray(mu, float), np.asarray(cov, float), size=n_subjects
        )
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


@pytest.fixture(scope="session")
def shipped_catalog():
    """
    The shipped catalog, loaded once for the whole session.

    Parsing it costs ~2.2s, and the selection tests need it dozens of times.
    Session scope is safe because nothing in these tests mutates a catalog.
    """
    from shedding_hub import load_shedding_catalog

    return load_shedding_catalog()


@pytest.fixture(scope="session")
def woelfel_dataset():
    """A real dataset for plot smoke tests, read from the repo rather than the network."""
    from shedding_hub import load_dataset

    return load_dataset("woelfel2020virological", local="./data")


@pytest.fixture
def ct_dataset():
    """One analyte in cycle threshold, three subjects that rise then fall."""
    curve = [
        (1, 30.0),
        (3, 26.0),
        (5, 25.0),
        (8, 27.0),
        (12, 30.0),
        (18, 33.0),
        (25, 36.0),
        (30, "negative"),
    ]
    shifts = (0.0, 1.5, -1.0)

    participants = []
    for shift in shifts:
        measurements = []
        for time, value in curve:
            shifted = value if isinstance(value, str) else value + shift
            measurements.append({"analyte": "swab", "time": time, "value": shifted})
        participants.append({"measurements": measurements})

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
        "participants": participants,
    }
