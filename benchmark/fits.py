"""Downstream equivalence of drafts and truth (Analysis F, benchmark/PLAN.md).

For each study and each analyte of truth, the package's three curve families
(`exponential`, `gamma`, `gamma_shifted`) are fitted with
`shedding_hub.fit_shedding_model` to truth and to the draft's matching
analyte, and the population summaries are compared: peak day, peak height
(log10, or cycles for a cycle-threshold analyte) and half-life. The report
says how far apart the two fits are, not what the curves mean.

Truth is built as a dataset: the reference, with each difference that
adjudication found in the draft's favour (verdict `S0`) applied to it:

    wrong analyte field             the draft's value
    wrong value, time or status     the draft's reading
    omission                        the reading removed
    extra reading                   the draft's reading added
    missing participant or analyte  removed
    extra participant or analyte    the draft's added, with analyte keys
                                    mapped to the reference's

Every other verdict leaves the reference as it is. Attributes are left alone,
since the fits do not use them. A correct value written by the adjudicator
for "both wrong" is free text and is not applied.

Fits are not reproducible across platforms, so this runs only in the pinned
environment: Linux, Python 3.11, and the numpy and scipy versions pinned in
requirements.txt, as the Build workflow uses. `--any-environment` overrides
the check for testing; the environment is recorded with the results either
way. One way to run it on another machine (Docker):

    docker run --rm -v "$PWD":/w -v RUNS:/runs -w /w python:3.11 bash -c \\
      "pip install -q -r requirements.txt && pip install -q -e . && \\
       python benchmark/fits.py /runs/RUN_ID --comparison /runs/RUN_ID.cmp \\
       [--adjudication /runs/RUN_ID.adj] --out /runs/RUN_ID.fits"

Usage (from the repository root):
    python benchmark/fits.py RUN_DIR --comparison CMP_DIR [--adjudication ADJ_DIR]
        --out OUT_DIR [--era manual]
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import platform
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(REPO_ROOT))
import adjudication as adj  # noqa: E402
import align  # noqa: E402
import discrepancies as disc  # noqa: E402
import metrics  # noqa: E402
from shedding_hub import SheddingDataError, fit_shedding_model  # noqa: E402
from shedding_hub.shedding_models import MODELS  # noqa: E402

PINNED_PYTHON = "3.11"
PINNED_PACKAGES = ("numpy", "scipy")
SUMMARIES = ("peak_day", "peak_log10", "half_life_days")

FIT_COLUMNS = [
    "study_id",
    "analyte",
    "draft_analyte",
    "model",
    "side",
    "status",
    "converged",
    "value_type",
    "n_subjects",
    "n_measurements",
    *SUMMARIES,
]
PAIR_COLUMNS = [
    "study_id",
    "analyte",
    "model",
    "value_type",
    "peak_day_difference",
    "peak_log10_difference",
    "half_life_ratio",
]
SUMMARY_COLUMNS = [
    "model",
    "value_type",
    "quantity",
    "estimate",
    "ci_low",
    "ci_high",
    "analytes",
    "studies",
]
CONVERGENCE_COLUMNS = ["model", "both", "truth_only", "draft_only", "neither"]


# --------------------------------------------------------------- environment


def pinned_versions(requirements: Path) -> dict:
    pins = {}
    for line in requirements.read_text(encoding="utf-8").splitlines():
        found = re.match(r"^([A-Za-z0-9_.-]+)==([^\s;#]+)", line)
        if found and found.group(1).lower() in PINNED_PACKAGES:
            pins[found.group(1).lower()] = found.group(2)
    return pins


def environment() -> dict:
    import scipy

    return {
        "system": platform.system(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }


def environment_problems(actual: dict, pins: dict) -> list[str]:
    """How the running environment differs from the pinned one."""
    problems = []
    if actual["system"] != "Linux":
        problems.append(f"system is {actual['system']}, not Linux")
    if not actual["python"].startswith(PINNED_PYTHON + "."):
        problems.append(f"Python is {actual['python']}, not {PINNED_PYTHON}")
    for name, version in pins.items():
        if actual.get(name) != version:
            problems.append(f"{name} is {actual.get(name)}, not {version}")
    return problems


# --------------------------------------------------------------------- truth


def _parts(path: str) -> list[str]:
    return [p.replace("~1", "/").replace("~0", "~") for p in path.split("/")[1:]]


def truth_dataset(s0: dict, s2: dict, rows: list[dict], verdict: dict, tol) -> dict:
    """The reference with every difference judged in the draft's favour applied."""
    study = adj.Study(s0, s2, tol)
    truth = copy.deepcopy(s2)
    analytes, people = truth["analytes"], truth["participants"]
    drop_readings, drop_people, drop_analytes = set(), set(), set()
    added_people = []

    def new_key(k0: str) -> str:
        key, n = k0, 2
        while key in analytes:
            key, n = f"{k0}_{n}", n + 1
        return key

    for row in rows:
        if verdict.get(row["id"]) != "S0":
            continue
        kind, level = row["type"], row["level"]
        s0_parts, s2_parts = _parts(row["s0_path"]), _parts(row["s2_path"])
        if level == "analyte":
            k0, name = _parts(row["s0_anchor"])[1:3]
            k2 = study.pairs[k0]
            if name in study.a0[k0]:
                analytes[k2][name] = study.a0[k0][name]
            else:
                analytes[k2].pop(name, None)
        elif kind in ("wrong value", "wrong non-detect status", "wrong time"):
            i0, j0, leaf = int(s0_parts[1]), int(s0_parts[3]), s0_parts[4]
            j2_owner, j2 = int(s2_parts[1]), int(s2_parts[3])
            reading = study.p0[i0]["measurements"][j0]
            people[j2_owner]["measurements"][j2][leaf] = reading[leaf]
        elif kind == "omission":
            drop_readings.add((int(s2_parts[1]), int(s2_parts[3])))
        elif kind == "extra" and level == "measurement":
            i0, j0 = int(s0_parts[1]), int(s0_parts[3])
            reading = copy.deepcopy(study.p0[i0]["measurements"][j0])
            k2 = study.pairs.get(reading.get("analyte"))
            if k2 is not None and study.people.get(i0) is not None:
                reading["analyte"] = k2
                people[study.people[i0]]["measurements"].append(reading)
        elif kind == "missing" and level == "participant":
            drop_people.add(int(s2_parts[1]))
        elif kind == "extra" and level == "participant":
            participant = copy.deepcopy(study.p0[int(s0_parts[1])])
            participant["measurements"] = [
                {**m, "analyte": study.pairs[m.get("analyte")]}
                for m in participant.get("measurements") or []
                if m.get("analyte") in study.pairs
            ]
            added_people.append(participant)
        elif kind == "missing analyte":
            drop_analytes.add(s2_parts[1])
        elif kind == "extra analyte":
            k0 = s0_parts[1]
            key = new_key(k0)
            analytes[key] = copy.deepcopy(study.a0[k0])
            for i0, participant in enumerate(study.p0):
                j2 = study.people.get(i0)
                if j2 is None:
                    continue
                for m in participant.get("measurements") or []:
                    if m.get("analyte") == k0:
                        people[j2]["measurements"].append({**m, "analyte": key})

    for j, participant in enumerate(people):
        participant["measurements"] = [
            m
            for k, m in enumerate(participant["measurements"])
            if (j, k) not in drop_readings and m.get("analyte") not in drop_analytes
        ]
    truth["participants"] = [
        p for j, p in enumerate(people) if j not in drop_people
    ] + added_people
    for key in drop_analytes:
        analytes.pop(key, None)
    return truth


# ---------------------------------------------------------------------- fits


def fit_one(dataset: dict | None, analyte: str | None, model: str) -> dict:
    """The summaries of one fit, or why there is none."""
    if dataset is None or analyte is None:
        return {"status": "no analyte"}
    try:
        # The package warns about every excluded subject; the counts it needs
        # are in the fit itself.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fit = fit_shedding_model(dataset, analyte=analyte, model=model)
    except SheddingDataError as exc:
        return {"status": f"not fitted: {exc.reason}"}
    except Exception as exc:  # a failed fit is a result, not a crash
        return {"status": f"error: {type(exc).__name__}"}
    out = {
        "status": "fitted",
        "converged": bool(fit.converged),
        "value_type": fit.value_type,
        "n_subjects": fit.n_subjects,
        "n_measurements": fit.n_measurements,
    }
    for name in SUMMARIES:
        try:
            out[name] = float(getattr(fit, name))
        except Exception:
            out[name] = math.nan
    return out


def usable(fit: dict) -> bool:
    return fit.get("status") == "fitted" and bool(fit.get("converged"))


def compare_fits(truth: dict, draft: dict) -> dict:
    """Differences of the summaries of two usable fits."""
    same_type = truth["value_type"] == draft["value_type"]

    def difference(name):
        return draft.get(name, math.nan) - truth.get(name, math.nan)

    half = math.nan
    if same_type and truth["half_life_days"] > 0 and draft["half_life_days"] > 0:
        half = draft["half_life_days"] / truth["half_life_days"]
    return {
        "value_type": truth["value_type"] if same_type else "mixed",
        "peak_day_difference": difference("peak_day"),
        # Peak height and half-life compare only within one value type.
        "peak_log10_difference": difference("peak_log10") if same_type else math.nan,
        "half_life_ratio": half,
    }


def median_by_study(values_by_study: dict, seed: int, replicates: int):
    """Median over analytes, with a 95% interval resampling studies."""
    studies = [v for v in values_by_study.values() if v]
    pooled = [x for v in studies for x in v]
    if not pooled:
        return math.nan, math.nan, math.nan
    rng = np.random.default_rng(seed)
    medians = []
    for _ in range(replicates):
        draw = rng.integers(0, len(studies), len(studies))
        medians.append(np.median([x for i in draw for x in studies[i]]))
    low, high = np.percentile(medians, [2.5, 97.5])
    return float(np.median(pooled)), float(low), float(high)


def summarize(pairs: list[dict], seed: int, replicates: int) -> list[dict]:
    """Per model and value type, so that log10 units and cycles never mix."""
    out = []
    quantities = (
        ("peak day difference, days (median)", "peak_day_difference", None),
        ("absolute peak day difference, days (median)", "peak_day_difference", abs),
        ("peak height difference (median)", "peak_log10_difference", None),
        ("half-life ratio, draft / truth (median)", "half_life_ratio", None),
    )
    groups = sorted({(p["model"], p["value_type"]) for p in pairs})
    for model, value_type in groups:
        mine = [
            p for p in pairs if (p["model"], p["value_type"]) == (model, value_type)
        ]
        for label, column, transform in quantities:
            by_study = defaultdict(list)
            for p in mine:
                value = p[column]
                if not math.isnan(value):
                    by_study[p["study_id"]].append(
                        transform(value) if transform else value
                    )
            estimate, low, high = median_by_study(by_study, seed, replicates)
            if not by_study:
                continue
            out.append(
                {
                    "model": model,
                    "value_type": value_type,
                    "quantity": label,
                    "estimate": estimate,
                    "ci_low": low,
                    "ci_high": high,
                    "analytes": sum(len(v) for v in by_study.values()),
                    "studies": len(by_study),
                }
            )
    return out


def convergence(fits: list[dict]) -> list[dict]:
    """How many analytes give a usable fit in both versions, one, or neither."""
    sides = defaultdict(dict)
    for f in fits:
        sides[(f["study_id"], f["analyte"], f["model"])][f["side"]] = usable(f)
    out = []
    for model in MODELS:
        counts = {"both": 0, "truth_only": 0, "draft_only": 0, "neither": 0}
        for (_, _, m), s in sides.items():
            if m != model:
                continue
            t, d = s.get("truth", False), s.get("draft", False)
            key = (
                "both"
                if t and d
                else "truth_only" if t else "draft_only" if d else "neither"
            )
            counts[key] += 1
        out.append({"model": model, **counts})
    return out


def study_fits(study_id: str, s0: dict | None, truth: dict, tol) -> list[dict]:
    """Every fit of one study: each truth analyte, each model, both sides."""
    pairing = adj.Study(s0, truth, tol) if s0 is not None else None
    out = []
    for k2 in truth.get("analytes", {}):
        k0 = pairing.back.get(k2) if pairing else None
        for model in MODELS:
            for side, dataset, key in (("truth", truth, k2), ("draft", s0, k0)):
                out.append(
                    {
                        "study_id": study_id,
                        "analyte": k2,
                        "draft_analyte": k0 or "",
                        "model": model,
                        "side": side,
                        **fit_one(dataset, key, model),
                    }
                )
    return out


def pair_fits(fits: list[dict]) -> list[dict]:
    by_key = defaultdict(dict)
    for f in fits:
        by_key[(f["study_id"], f["analyte"], f["model"])][f["side"]] = f
    out = []
    for (study_id, analyte, model), s in sorted(by_key.items()):
        if usable(s.get("truth", {})) and usable(s.get("draft", {})):
            out.append(
                {
                    "study_id": study_id,
                    "analyte": analyte,
                    "model": model,
                    **compare_fits(s["truth"], s["draft"]),
                }
            )
    return out


# ------------------------------------------------------------------------ CLI


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--comparison", type=Path, required=True)
    ap.add_argument("--adjudication", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--era", default="manual", help="manual, ai_assisted or all")
    ap.add_argument("--ref", default=disc.RELEASE)
    ap.add_argument("--seed", type=int, default=metrics.SEED)
    ap.add_argument("--replicates", type=int, default=metrics.REPLICATES)
    ap.add_argument("--any-environment", action="store_true")
    args = ap.parse_args(argv)

    actual = environment()
    problems = environment_problems(
        actual, pinned_versions(REPO_ROOT / "requirements.txt")
    )
    if problems and not args.any_environment:
        for problem in problems:
            print(f"NOT THE PINNED ENVIRONMENT: {problem}", file=sys.stderr)
        return 2

    rows = metrics.read_csv(args.comparison / "discrepancies.csv")
    key, answers = [], {}
    if args.adjudication:
        key = metrics.read_csv(args.adjudication / "key.csv")
        answers, wrong = metrics.read_answers(
            args.adjudication / "discrepancy_packet.csv", metrics.DISCREPANCY_OUTCOMES
        )
        if wrong:
            for problem in wrong:
                print(f"UNKNOWN OUTCOME {problem}", file=sys.stderr)
            return 1
    verdict = metrics.verdicts(rows, key, answers)
    by_study = defaultdict(list)
    for row in rows:
        by_study[row["study_id"]].append(row)

    manifest = disc.read_manifest()
    fits = []
    for study_dir in sorted(d for d in args.run_dir.iterdir() if d.is_dir()):
        info = manifest.get(study_dir.name)
        if info is None or (args.era != "all" and info["era"] != args.era):
            continue
        s2 = disc.load_reference(info["reference_yaml"], args.ref)
        tol = align.tolerances_for(info["input_csv_origin"])
        draft = study_dir / "s0.yaml"
        s0 = None
        if draft.is_file():
            try:
                s0 = yaml.load(
                    draft.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader
                )
            except yaml.YAMLError:
                s0 = None
        s0 = s0 if isinstance(s0, dict) else None
        truth = (
            truth_dataset(s0, s2, by_study[study_dir.name], verdict, tol)
            if s0 is not None
            else s2
        )
        fits += study_fits(study_dir.name, s0, truth, tol)

    pairs = pair_fits(fits)
    args.out.mkdir(parents=True, exist_ok=True)
    disc.write_csv(args.out / "fits.csv", FIT_COLUMNS, fits)
    disc.write_csv(args.out / "fit_pairs.csv", PAIR_COLUMNS, pairs)
    disc.write_csv(
        args.out / "fit_summary.csv",
        SUMMARY_COLUMNS,
        summarize(pairs, args.seed, args.replicates),
    )
    disc.write_csv(
        args.out / "fit_convergence.csv", CONVERGENCE_COLUMNS, convergence(fits)
    )
    (args.out / "environment.json").write_text(
        json.dumps(
            {**actual, "pinned": not problems, "differences": problems}, indent=2
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        f"{len({f['study_id'] for f in fits})} studies, {len(fits)} fits, "
        f"{len(pairs)} compared -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
