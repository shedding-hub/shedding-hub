from .util import check_dataset, folded_str, literal_str, load_dataset, normalize_str
from .shedding_duration import (
    calc_shedding_duration,
    calc_shedding_durations,
    plot_shedding_duration,
    plot_shedding_durations,
)

from .shedding_peak import (
    calc_shedding_peak,
    calc_shedding_peaks,
    plot_shedding_peak,
    plot_shedding_peaks,
)

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

from .stats import (
    calc_shedding_summary,
    calc_detection_summary,
    calc_clearance_summary,
    calc_value_summary,
    calc_dataset_summary,
    compare_datasets,
)

from .shedding_models import MODELS, PARAM_NAMES

from .shedding_fit import SheddingDataError, SheddingFit, fit_shedding_model

from .shedding_catalog import (
    SheddingCatalog,
    fit_shedding_models,
    load_shedding_catalog,
)

from .shedding_ensemble import SheddingEnsemble, make_ensemble

from .shedding_select import (
    REFERENCE_EVENT_CLASSES,
    Selection,
    classify_reference_event,
    shedding_for,
    shedding_options,
)

from .shedding_simulate import plot_simulated_shedding, simulate_shedding

__all__ = [
    "check_dataset",
    "folded_str",
    "literal_str",
    "load_dataset",
    "normalize_str",
    "calc_shedding_duration",
    "calc_shedding_durations",
    "plot_shedding_duration",
    "plot_shedding_durations",
    "calc_shedding_peak",
    "calc_shedding_peaks",
    "plot_shedding_peak",
    "plot_shedding_peaks",
    "plot_time_course",
    "plot_time_courses",
    "plot_shedding_heatmap",
    "plot_mean_trajectory",
    "plot_catalog_fits",
    "plot_fit_diagnostic",
    "plot_clearance_curve",
    "plot_detection_probability",
    "plot_value_distribution_by_time",
    "calc_shedding_summary",
    "calc_detection_summary",
    "calc_clearance_summary",
    "calc_value_summary",
    "calc_dataset_summary",
    "compare_datasets",
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
    "REFERENCE_EVENT_CLASSES",
    "Selection",
    "classify_reference_event",
    "shedding_for",
    "shedding_options",
    "simulate_shedding",
    "plot_simulated_shedding",
]
