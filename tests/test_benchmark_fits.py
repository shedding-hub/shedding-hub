"""Analysis F: truth is built from the verdicts, and equal data fit equally."""

import copy
import csv
import math
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import align  # noqa: E402
import discrepancies as disc  # noqa: E402
import fits  # noqa: E402
import metrics  # noqa: E402


def scored_data_rows(rows):
    """Scored differences other than attributes, which truth leaves alone."""
    return [
        r
        for r in rows
        if r["scored"] == "yes" and r["type"] not in metrics.ATTRIBUTE_TYPES
    ]


def load(study):
    path = REPO_ROOT / "data" / study / f"{study}.yaml"
    return yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader)


def numeric_reading(participant):
    for j, m in enumerate(participant["measurements"]):
        if align.is_number(m["value"]) and not any(
            k != j and (o["time"] == m["time"] or o["value"] == m["value"])
            for k, o in enumerate(participant["measurements"])
        ):
            return j
    return None


def tan_draft():
    """Every kind of measurement and participant difference, one of each."""
    reference = load("tan2021early")
    draft = copy.deepcopy(reference)
    people = [
        (i, j)
        for i, p in enumerate(reference["participants"])
        if len(p["measurements"]) >= 3 and (j := numeric_reading(p)) is not None
    ]
    (a, ja), (b, jb), (c, jc), (d, jd), (e, je) = people[:5]
    ms = draft["participants"]
    ms[a]["measurements"][ja]["value"] += 9
    ms[b]["measurements"][jb]["time"] += 40
    del ms[c]["measurements"][jc]
    ms[d]["measurements"].append(dict(ms[d]["measurements"][0], time=900, value=21.0))
    ms[e]["measurements"][je]["value"] = "negative"
    key = next(iter(draft["analytes"]))
    draft["analytes"][key]["limit_of_detection"] = 40
    newcomer = copy.deepcopy(ms[0])
    for m in newcomer["measurements"]:
        m["time"] += 500
    del draft["participants"][-1]
    draft["participants"].append(newcomer)
    return draft, reference


def antinori_draft():
    """One analyte dropped and one added."""
    reference = load("antinori2022epidemiological")
    draft = copy.deepcopy(reference)
    gone = "plasma"
    del draft["analytes"][gone]
    for p in draft["participants"]:
        p["measurements"] = [m for m in p["measurements"] if m["analyte"] != gone]
    draft["analytes"]["urine"] = dict(reference["analytes"]["stool"], specimen="urine")
    first = draft["participants"][0]["measurements"][0]
    draft["participants"][0]["measurements"].append(dict(first, analyte="urine"))
    return draft, reference


def build_truth(study, draft, reference, verdict_for):
    rows = disc.compare(study, draft, reference, align.REPORTED).rows
    verdict = {r["id"]: verdict_for(r) for r in rows if r["scored"] == "yes"}
    return fits.truth_dataset(draft, reference, rows, verdict, align.REPORTED), rows


def test_the_tan_draft_covers_every_measurement_and_participant_kind():
    draft, reference = tan_draft()
    rows = disc.compare("tan2021early", draft, reference, align.REPORTED).rows
    found = {(r["level"], r["type"]) for r in rows}
    assert {
        ("measurement", "wrong value"),
        ("measurement", "wrong time"),
        ("measurement", "omission"),
        ("measurement", "extra"),
        ("measurement", "wrong non-detect status"),
        ("analyte", "wrong limit_of_detection"),
        ("participant", "missing"),
        ("participant", "extra"),
    } <= found


def test_with_every_difference_in_the_drafts_favour_truth_is_the_draft():
    for study, (draft, reference) in (
        ("tan2021early", tan_draft()),
        ("antinori2022epidemiological", antinori_draft()),
    ):
        truth, rows = build_truth(study, draft, reference, lambda r: "S0")
        assert len({r["type"] for r in scored_data_rows(rows)}) >= 2
        left = disc.compare(study, draft, truth, align.REPORTED).rows
        assert scored_data_rows(left) == []


def test_any_other_verdict_leaves_the_reference():
    draft, reference = tan_draft()
    for verdict in ("S2", "both", "neither", "curator", "reference"):
        truth, _ = build_truth("tan2021early", draft, reference, lambda r: verdict)
        assert truth == reference


def test_identical_data_fit_identically_across_renamed_analytes():
    reference = load("tan2021early")
    draft = copy.deepcopy(reference)
    draft["analytes"] = {"renamed": draft["analytes"].pop("swab_SARSCoV2_N")}
    for p in draft["participants"]:
        for m in p["measurements"]:
            m["analyte"] = "renamed"
    rows = fits.study_fits("tan2021early", draft, reference, align.REPORTED)
    assert {r["draft_analyte"] for r in rows} == {"renamed"}
    pairs = fits.pair_fits(rows)
    assert len(pairs) == len(fits.MODELS)
    for p in pairs:
        assert p["peak_day_difference"] == 0
        assert p["peak_log10_difference"] == 0
        assert p["half_life_ratio"] == 1
    assert {c["model"]: c["both"] for c in fits.convergence(rows)} == {
        m: 1 for m in fits.MODELS
    }


def test_a_missing_draft_fits_one_side_only():
    reference = load("tan2021early")
    rows = fits.study_fits("tan2021early", None, reference, align.REPORTED)
    assert {r["status"] for r in rows if r["side"] == "draft"} == {"no analyte"}
    assert fits.pair_fits(rows) == []
    counts = {c["model"]: c for c in fits.convergence(rows)}
    assert counts["gamma"]["truth_only"] == 1


def test_units_never_mix_in_the_summary():
    def pair(study, value_type, height):
        return {
            "study_id": study,
            "analyte": "a",
            "model": "gamma",
            "value_type": value_type,
            "peak_day_difference": 1.0,
            "peak_log10_difference": height,
            "half_life_ratio": 1.5,
        }

    rows = fits.summarize(
        [pair("s1", "ct", 2.0), pair("s2", "concentration", 0.1)], 1, 50
    )
    height = {
        r["value_type"]: r["estimate"]
        for r in rows
        if r["quantity"].startswith("peak height")
    }
    assert height == {"ct": 2.0, "concentration": 0.1}


def test_median_by_study_resamples_studies():
    assert fits.median_by_study({"s": [1.0, 2.0, 3.0]}, 1, 50) == (2.0, 2.0, 2.0)
    estimate, low, high = fits.median_by_study(
        {"a": [0.0], "b": [10.0], "c": [5.0]}, 1, 200
    )
    assert estimate == 5.0 and low <= estimate <= high
    assert all(math.isnan(x) for x in fits.median_by_study({}, 1, 10))


def test_the_environment_check():
    pins = fits.pinned_versions(REPO_ROOT / "requirements.txt")
    assert set(pins) == {"numpy", "scipy"}
    good = {"system": "Linux", "python": "3.11.9", **pins}
    assert fits.environment_problems(good, pins) == []
    bad = {**good, "system": "Windows", "python": "3.12.1", "numpy": "1.26.4"}
    assert len(fits.environment_problems(bad, pins)) == 3


def test_the_command_line_refuses_an_unpinned_environment(tmp_path, monkeypatch):
    unpinned = {"system": "Windows", "python": "3.12.1", "numpy": "1", "scipy": "1"}
    monkeypatch.setattr(fits, "environment", lambda: unpinned)
    argv = [str(tmp_path), "--comparison", str(tmp_path), "--out", str(tmp_path / "o")]
    assert fits.main(argv) == 2
    assert not (tmp_path / "o").exists()


def test_the_command_line_fits_a_run(tmp_path):
    draft, _ = tan_draft()
    run = tmp_path / "run"
    (run / "tan2021early").mkdir(parents=True)
    (run / "tan2021early" / "s0.yaml").write_text(
        yaml.safe_dump(draft), encoding="utf-8"
    )
    cmp_dir, out = tmp_path / "cmp", tmp_path / "out"
    assert disc.main([str(run), "--out", str(cmp_dir), "--ref", "HEAD"]) == 0
    argv = [str(run), "--comparison", str(cmp_dir), "--out", str(out), "--ref", "HEAD"]
    assert fits.main(argv + ["--any-environment", "--replicates", "20"]) == 0
    rows = list(csv.DictReader((out / "fits.csv").open(encoding="utf-8")))
    assert len(rows) == 2 * len(fits.MODELS)
    assert {r["side"] for r in rows} == {"truth", "draft"}
    for name in ("fit_pairs.csv", "fit_summary.csv", "fit_convergence.csv"):
        assert (out / name).is_file()
    assert '"pinned"' in (out / "environment.json").read_text(encoding="utf-8")
