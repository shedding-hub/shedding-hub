import collections
import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import manifest as mf  # noqa: E402


def _rows():
    with (REPO_ROOT / "benchmark" / "manifest.csv").open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _dataset(method, sources, models=None, version=None):
    curation = {"method": method, "data_source": sources}
    if models:
        curation["models"] = models
    if version:
        curation["pipeline_version"] = version
    return {"doi": "10.1/x", "curation": curation}


def test_manifest_covers_the_analysis_set_once():
    rows = _rows()
    ids = [r["study_id"] for r in rows]
    assert len(ids) == len(set(ids)) == 144
    assert not set(mf.EXCLUDED) & set(ids)
    assert set(mf.EXCLUDED) == {"jones2021estimating", "cdc2024nhphrn"}
    assert list(rows[0]) == mf.COLUMNS


def test_manifest_era_and_pipeline_counts():
    rows = _rows()
    assert collections.Counter(r["era"] for r in rows) == {
        "manual": 37,
        "ai_assisted": 107,
    }
    versions = collections.Counter(r["pipeline_version"] for r in rows)
    assert versions == {"": 37, "v1": 16, "v2": 91}
    # Analysis A is the manual era, and only it.
    assert all((r["in_analysis_A"] == "yes") == (r["era"] == "manual") for r in rows)
    # A pipeline version and an extraction model go together.
    assert all(bool(r["pipeline_version"]) == bool(r["extraction_model"]) for r in rows)


def test_manifest_origin_is_consistent_with_the_csv_flag():
    for r in _rows():
        if r["input_csv_origin"] == "table_supplement":
            assert r["has_figure_csv"] == "yes"
        if r["input_csv_origin"] == "none":
            assert r["has_figure_csv"] == "no"


def test_origin_precedence_and_mixed_note():
    row = mf.manifest_row(
        "s",
        _dataset("ai_assisted", ["main_text", "figure", "author_shared"]),
        "1",
        True,
    )
    assert row["input_csv_origin"] == "author_shared"
    assert "author_shared + figure_digitized" in row["notes"]


def test_origin_for_a_csv_without_figure_or_repository_source():
    with_csv = mf.manifest_row("s", _dataset("ai_assisted", ["table"]), "1", True)
    without = mf.manifest_row("s", _dataset("ai_assisted", ["table"]), "1", False)
    assert with_csv["input_csv_origin"] == "table_supplement"
    assert without["input_csv_origin"] == "none"
    assert without["notes"] == ""


def test_missing_csv_and_missing_pmid_are_noted():
    row = mf.manifest_row("s", _dataset("manual", ["data_repository"]), "", False)
    assert row["input_csv_origin"] == "data_repository"
    assert row["notes"] == "no input CSV yet; no PMID"
    assert row["in_analysis_A"] == "yes"


def test_extraction_model_is_read_from_the_models_list():
    dataset = _dataset(
        "ai_assisted",
        ["table"],
        models=["gpt-5.2 (extraction)", "claude-opus-4-8 (review)"],
        version="v2",
    )
    row = mf.manifest_row("s", dataset, "1", True)
    assert (row["pipeline_version"], row["extraction_model"]) == ("v2", "gpt-5.2")


def test_ascii_id_strips_accents():
    assert mf.ascii_id("coppée2023temporal") == "coppee2023temporal"
    assert mf.ascii_id("suñer2023viral") == "suner2023viral"


def test_reference_is_the_released_file_unless_overridden():
    rows = {r["study_id"]: r for r in _rows()}
    for study, row in rows.items():
        if study not in mf.REFERENCE_OVERRIDES:
            assert row["reference_yaml"] == f"data/{study}/{study}.yaml"
    path, reason = mf.REFERENCE_OVERRIDES["obara2008single"]
    assert rows["obara2008single"]["reference_yaml"] == path
    assert reason in rows["obara2008single"]["notes"]
    assert rows["obara2008single"]["in_analysis_A"] == "yes"


def test_every_study_has_full_text_evidence_and_manual_studies_have_a_csv():
    rows = _rows()
    assert {r["evidence_tier"] for r in rows} <= {"pmc", "publisher", "manual-pdf"}
    assert all(r["has_figure_csv"] == "yes" for r in rows if r["era"] == "manual")
