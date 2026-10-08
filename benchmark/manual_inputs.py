"""
Build the curator-style input CSVs for the manual-era benchmark studies.

In production, a curator prepared one CSV per study for the extraction agent:
one row per measurement, starting with `PatientID,time,value`, plus whatever
other columns the source carried. The manual-era studies never had such a
file, so the benchmark builds one from the raw files stored with each dataset
at the frozen release.

A CSV is made in two steps. A builder lays the raw file out, and
`manual_curation.py` then applies the two things a production curator did
before handing a file to the agent: writing non-detects as `negative` and
undoing log10, each by the rule in the study's extraction script.

What a builder may do is limited to what laying out a CSV requires:

- pick the rows that belong to the study and the raw columns that hold the
  measurements, the time and the participant;
- reshape to one row per measurement;
- rename those three columns to `PatientID`, `time` and `value`;
- carry other raw columns along under their raw names;
- number participants 1, 2, 3, ... where neither the source file nor the paper
  gives an identifier, and put anything that describes a participant (mother
  or neonate, for example) in a column of its own;
- drop cells that hold no measurement.

A builder may not change a value. Beyond the two curation steps, times are
not re-aligned to a reference event, units are not converted and no label is
renamed. Nothing is taken from the released YAML.

Each builder returns the table and a note. The note records which raw columns
were used and anything in the raw file that already reflects a curator's
judgement; the curation step adds what it changed. Notes are written to
`benchmark/manual_inputs/NOTES.csv`.
"""

import argparse
import csv
import io
import pathlib
import subprocess
import sys

import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO_ROOT / "benchmark" / "manual_inputs"


def raw(release: str, study: str, name: str) -> io.BytesIO:
    """A raw file stored with a dataset, read from the release tag."""
    blob = subprocess.run(
        ["git", "show", f"{release}:data/{study}/{name}"],
        capture_output=True,
        cwd=REPO_ROOT,
        check=True,
    ).stdout
    return io.BytesIO(blob)


def lay_out(table: pd.DataFrame, patient: str, time: str, value: str) -> pd.DataFrame:
    """Put the three standard columns first, under their standard names."""
    renamed = table.rename(columns={patient: "PatientID", time: "time", value: "value"})
    others = [c for c in renamed.columns if c not in ("PatientID", "time", "value")]
    return renamed[["PatientID", "time", "value", *others]]


# The combined dataset's LOD is in the unit of its own `value` column, viral
# copies per mL, while the measurements taken from it are cycle thresholds.
LOD_COLUMN = "LOD_viral_copies_per_mL"


def number_participants(labels) -> list[str]:
    """Generic ids 1, 2, 3, ... in order of first appearance."""
    order = {}
    return [str(order.setdefault(label, len(order) + 1)) for label in labels]


def zuo2020alterations(release):
    table = pd.read_csv(raw(release, "zuo2020alterations", "data.csv"), dtype=str)
    note = (
        "data.csv, digitized from Fig. 3: ID, day and value as PatientID, time and "
        "value; Sex and Age kept."
    )
    return lay_out(table, "ID", "day", "value"), note


def iwakiri2009quantitative(release):
    table = pd.read_csv(raw(release, "iwakiri2009quantitative", "data.csv"), dtype=str)
    table["Specimens"] = table["Specimens"].str.strip()
    note = (
        "data.csv, transcribed from Table 1: Specimens, days_after_onset and value "
        "as PatientID, time and value; Outbreak, Sex and Age kept. The raw time "
        "column is named days_after_onset, a curator-supplied reference event."
    )
    return lay_out(table, "Specimens", "days_after_onset", "value"), note


def fajnzylber2020sars(release):
    sheet = pd.read_excel(
        raw(release, "fajnzylber2020sars", "41467_2020_19057_MOESM3_ESM.xlsx"),
        sheet_name="Inpatient",
        dtype=str,
    )
    loads = [c for c in sheet.columns if c.endswith("_VL")]
    keep = ["PID", "Sx_Onset_to_SC", "Age", "Sex", "Race_Ethnicity"]
    table = sheet.melt(
        id_vars=keep, value_vars=loads, var_name="specimen", value_name="value"
    )
    table = table[table["value"].notna()]
    # Stable order: participant, then time, then the sheet's column order.
    table["_time"] = pd.to_numeric(table["Sx_Onset_to_SC"])
    table["_col"] = table["specimen"].map({c: i for i, c in enumerate(loads)})
    table = table.sort_values(["PID", "_time", "_col"], kind="stable")
    table = table.drop(columns=["_time", "_col"])
    note = (
        "Supplement xlsx, sheet Inpatient: one row per non-empty viral-load cell of "
        f"the columns {', '.join(loads)}; the column name is kept in specimen. PID "
        "and Sx_Onset_to_SC as PatientID and time; Age, Sex and Race_Ethnicity "
        "kept. Clinical and laboratory columns are not carried. The raw time column "
        "is days from symptom onset to sample collection, as the authors supplied it."
    )
    return lay_out(table, "PID", "Sx_Onset_to_SC", "value"), note


def kissler2021viral(release):
    table = pd.read_csv(raw(release, "kissler2021viral", "ct_dat_clean.csv"), dtype=str)
    quantities = ["CT.Mean", "log10_GEperML"]
    keep = ["Person.ID", "Date.Index", "Symptomatic", "Novel.Persistent.Infection"]
    order = table[keep].drop_duplicates().reset_index(drop=True).reset_index()
    long = table.melt(
        id_vars=keep, value_vars=quantities, var_name="analyte", value_name="value"
    )
    long = long[long["value"].notna()]
    long = long.merge(order, on=keep, how="left")
    long["_col"] = long["analyte"].map({c: i for i, c in enumerate(quantities)})
    long = long.sort_values(["index", "_col"], kind="stable")
    long = long.drop(columns=["index", "_col"])
    note = (
        "ct_dat_clean.csv from the authors' repository: one row per sample and "
        "quantity, for the columns CT.Mean and log10_GEperML; the column name is "
        "kept in analyte. Person.ID and Date.Index as PatientID and time; "
        "Symptomatic and Novel.Persistent.Infection kept. N1_CT_Value, CT.T1 and "
        "Adjusted are not carried."
    )
    return lay_out(long, "Person.ID", "Date.Index", "value"), note


def salvatore2020epidemiological(release):
    sheet = pd.read_excel(
        raw(release, "salvatore2020epidemiological", "CombinedDataset.xlsx"),
        sheet_name="Viral_Load",
        dtype=str,
    )
    table = sheet[sheet["StudyNum"] == "17"]
    table = table[["PatientID", "Day", "Ctvalue", "Age", "Sex", "LOD"]]
    table = table[table["Ctvalue"].notna()]
    table = table.rename(columns={"LOD": LOD_COLUMN})
    note = (
        "CombinedDataset.xlsx (Challenger et al., a third-party combined dataset), "
        "sheet Viral_Load, rows with StudyNum 17: Day and Ctvalue as time and "
        "value; Age and Sex kept, and LOD kept as LOD_viral_copies_per_mL, since "
        "the combined dataset states its limit of detection in viral copies per "
        "mL and not in cycles. The combined dataset was already harmonized across "
        "studies by its authors."
    )
    return lay_out(table, "PatientID", "Day", "Ctvalue"), note


BUILDERS = {
    f.__name__: f
    for f in (
        fajnzylber2020sars,
        iwakiri2009quantitative,
        kissler2021viral,
        salvatore2020epidemiological,
        zuo2020alterations,
    )
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("studies", nargs="*", help="Studies to build (default: all)")
    parser.add_argument("--release", default="v1.1.0", help="Frozen release tag.")
    args = parser.parse_args()

    from manual_curation import curate

    # Builders are grouped by kind of source, one module per group.
    for module in (
        "manual_builders_challenger",
        "manual_builders_repository",
        "manual_builders_figure",
        "manual_builders_legacy",
    ):
        try:
            BUILDERS.update(__import__(module).BUILDERS)
        except ModuleNotFoundError:
            pass

    selected = args.studies or sorted(BUILDERS)
    OUT.mkdir(parents=True, exist_ok=True)
    notes_path = OUT / "NOTES.csv"
    notes = {}
    if notes_path.exists():
        with notes_path.open(encoding="utf-8", newline="") as fh:
            notes = {row["study_id"]: row for row in csv.DictReader(fh)}

    for study in selected:
        table, note = BUILDERS[study](args.release)
        table, curated = curate(study, table)
        note = f"{note} {curated}"
        table.to_csv(OUT / f"{study}.csv", index=False, lineterminator="\n")
        notes[study] = {
            "study_id": study,
            "rows": len(table),
            "participants": table["PatientID"].nunique(),
            "columns": ";".join(table.columns),
            "note": note,
        }
        print(
            f"{study:32s} {len(table):5d} rows  {table['PatientID'].nunique():4d} ids"
        )

    with notes_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["study_id", "rows", "participants", "columns", "note"],
            lineterminator="\n",
        )
        writer.writeheader()
        for study in sorted(notes):
            writer.writerow(notes[study])
    return 0


if __name__ == "__main__":
    sys.exit(main())
