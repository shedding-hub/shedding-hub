"""
Input CSV builders for studies curated from the Challenger combined dataset,
and for two studies whose raw files a curator had already tabulated.

`CombinedDataset.xlsx` is a third-party dataset (Challenger et al.) that
pools individual-level viral load data from many studies and was harmonized
by its authors. Its sheet `Viral_Load` has one row per sample. The curators
read `Ctvalue` as the measurement and `Day` as the time for every study taken
from it. The sheet's own `value` column (viral copies per mL, set to 1 below
the limit of detection) is carried along as `CombinedDataset_value`, because
the name `value` is taken by the measurement.
"""

import io
import subprocess

import pandas as pd

from manual_inputs import REPO_ROOT, lay_out, raw

# study -> StudyNum in the combined dataset, as each extraction script filters.
STUDY_NUMBERS = {
    "alsharrah2020clinical": 15,
    "gautret2020hydroxychloroquine": 5,
    "salvatore2020epidemiological": 17,
    "shrestha2020distribution": 12,
    "tan2021early": 16,
    "team2020clinical": 7,
    "yilmaz2020upper": 14,
    "young2020epidemiologic": 6,
}

COMBINED_NOTE = (
    "CombinedDataset.xlsx (Challenger et al., a third-party combined dataset "
    "already harmonized across studies), sheet Viral_Load, rows with StudyNum "
    "{number}: Day and Ctvalue as time and value; Age, Sex and LOD kept; the "
    "sheet's own value column (viral copies per mL, 1 below the limit of "
    "detection) kept as CombinedDataset_value. Day is days since symptom onset "
    "as the combined dataset defines it."
)


def combined_rows(release: str, study: str, number: int) -> pd.DataFrame:
    sheet = pd.read_excel(
        raw(release, study, "CombinedDataset.xlsx"), sheet_name="Viral_Load", dtype=str
    )
    table = sheet[sheet["StudyNum"] == str(number)]
    table = table[["PatientID", "Day", "Ctvalue", "Age", "Sex", "LOD", "value"]]
    table = table.rename(columns={"value": "CombinedDataset_value"})
    return lay_out(table, "PatientID", "Day", "Ctvalue")


def combined_builder(study: str, number: int, extra: str = ""):
    def build(release):
        note = COMBINED_NOTE.format(number=number) + (" " + extra if extra else "")
        return combined_rows(release, study, number), note

    build.__name__ = study
    return build


def kimse2020viral(release):
    combined = combined_rows(release, "kimse2020viral", 4)
    combined["source_file"] = "CombinedDataset.xlsx"
    extra = pd.read_excel(
        raw(release, "kimse2020viral", "asymptomatic.xlsx"), dtype=str
    )
    extra = lay_out(extra, "PatientID", "Day", "value")
    extra["source_file"] = "asymptomatic.xlsx"
    table = pd.concat([combined, extra], ignore_index=True)
    note = (
        COMBINED_NOTE.format(number=4)
        + " Followed by the rows of asymptomatic.xlsx, a curator-digitized file: "
        "Day and value as time and value; viral load, Age and Sex kept. The file "
        "each row came from is in source_file. asymptomatic_demo.xlsx is not used."
    )
    return table, note


def kim2020viral(release):
    sheet = pd.read_excel(raw(release, "kim2020viral", "kim2020viral.xlsx"), dtype=str)
    sheet.columns = [c.strip() for c in sheet.columns]
    sheet = sheet.reset_index()
    quantities = ["Value", "Ctvalue"]
    keep = ["index", "PatientID", "Day", "Age", "Sex", "Type"]
    table = sheet.melt(
        id_vars=keep, value_vars=quantities, var_name="analyte", value_name="value"
    )
    table = table[table["value"].notna()]
    table["_col"] = table["analyte"].map({c: i for i, c in enumerate(quantities)})
    table["index"] = table["index"].astype(int)
    table = table.sort_values(["index", "_col"], kind="stable")
    table = table.drop(columns=["index", "_col"])
    note = (
        "kim2020viral.xlsx, a curator-prepared table: one row per non-empty cell "
        "of the columns Value and Ctvalue; the column name is kept in analyte. "
        "Day as time; Age, Sex and Type kept. The three supplementary xls files "
        "stored with the dataset are not used. The file was tabulated by a "
        "curator from the supplement and the Challenger combined dataset."
    )
    return lay_out(table, "PatientID", "Day", "value"), note


def obara2008single(release):
    blob = subprocess.run(
        [
            "git",
            "show",
            f"{release}:archived data/obara2008single/log10_copy_number_data.xlsx",
        ],
        capture_output=True,
        cwd=REPO_ROOT,
        check=True,
    ).stdout
    table = pd.read_excel(io.BytesIO(blob), dtype=str)
    note = (
        "archived data/obara2008single/log10_copy_number_data.xlsx, the raw file "
        "of the hand-curated version: Employee, Date and Log10 copy number as "
        "PatientID, time and value; Sample ID kept. Dates are as the curator "
        "typed them and carry a placeholder year. employee_pcr_data.xlsx, a "
        "transcription of Table 1, is not used."
    )
    return lay_out(table, "Employee", "Date", "Log10 copy number"), note


EXTRA = {
    "gautret2020hydroxychloroquine": (
        "The treatment table the extraction script typed in from the paper is "
        "not carried."
    ),
    "team2020clinical": (
        "The demographic table the extraction script typed in from the paper is "
        "not carried."
    ),
    "young2020epidemiologic": (
        "The released values are concentrations the curator computed from Ct "
        "with a regression fitted to Viral_Loads.csv; that conversion is not "
        "applied here and Viral_Loads.csv is not used."
    ),
}

BUILDERS = {
    study: combined_builder(study, number, EXTRA.get(study, ""))
    for study, number in STUDY_NUMBERS.items()
}
BUILDERS.update(
    {f.__name__: f for f in (kimse2020viral, kim2020viral, obara2008single)}
)
