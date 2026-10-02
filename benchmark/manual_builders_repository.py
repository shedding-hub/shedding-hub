"""
Input CSV builders for manual-era studies whose values came from a data
repository or a supplementary file.

The rules are the ones stated in `manual_inputs.py`: layout only. Each builder
selects the study's rows and the raw columns that hold the participant, the
time and the measurement, reshapes to one row per measurement, and carries
other raw columns along. No value is changed.

Three layout devices recur here and are worth stating once.

A wide file is made long by taking one row per measurement column. The raw
column name is kept in `analyte`, so nothing about the quantity is renamed.

Where a measurement column has a companion column in the raw file (a
detection flag, a result class), the companion is carried beside the value in
one column named after the raw suffix. The pairing is positional and comes
from the raw column names, not from any rule about what counts as detected.

Where the raw file records a test whose result sits in another column (a
negative PCR with an empty viral-load cell), the row is kept with an empty
`value` and the raw result column carried, rather than writing a non-detect
code into `value`.
"""

import io
import pathlib
import urllib.request

import pandas as pd

from manual_inputs import lay_out, raw

DOWNLOADS = pathlib.Path(__file__).resolve().parent / "manual_inputs" / "_downloads"


def read_csv(source, **kwargs) -> pd.DataFrame:
    """A CSV as text, cell for cell: only an empty cell counts as empty."""
    return pd.read_csv(
        source, dtype=str, keep_default_na=False, na_values=[""], **kwargs
    )


def interleave(parts: list[pd.DataFrame]) -> pd.DataFrame:
    """Stack per-column tables so each raw row's measurements stay together."""
    stacked = pd.concat(
        [part.assign(_row=range(len(part)), _col=i) for i, part in enumerate(parts)]
    )
    stacked = stacked.sort_values(["_row", "_col"], kind="stable")
    return stacked.drop(columns=["_row", "_col"]).reset_index(drop=True)


def download(study: str, name: str, url: str) -> io.BytesIO:
    """A source file that is not stored with the dataset, cached after one fetch."""
    cached = DOWNLOADS / study / name
    if not cached.exists():
        cached.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=60) as response:
            cached.write_bytes(response.read())
    return io.BytesIO(cached.read_bytes())


def arts2023longitudinal(release):
    workbook = raw(release, "arts2023longitudinal", "msphere.00132-23-s0002.xlsx")
    shedding = pd.read_excel(workbook, sheet_name="Shedding Data", dtype=str)
    people = pd.read_excel(workbook, sheet_name="Demographics", dtype=str)
    person_columns = [
        "Sex",
        "Age",
        "Symptomatic",
        "Index (I) or Household member (H)",
        "Vaccination Status",
    ]
    sample_columns = ["ID", "Day", "solid_frac", "extraction_mass_mg"]
    # measurement column -> its companion flag columns in the sheet
    quantities = {
        "N_conc (gc/mg-dw)": {"det": "N_det"},
        "ORF1a_conc (gc/gm-dw)": {"det": "ORF1a_det"},
        "PMMoV_conc (gc/mg-dw)": {"det": "PMMoV_det"},
        "crAss_conc": {"det": "crAss_det", "quant": "crAss_quant"},
    }
    parts = []
    for column, flags in quantities.items():
        part = shedding[sample_columns].copy()
        part["value"] = shedding[column]
        part["analyte"] = column
        part["det"] = shedding[flags["det"]]
        part["quant"] = shedding[flags["quant"]] if "quant" in flags else pd.NA
        parts.append(part)
    table = interleave(parts)
    # A sample contributes a row for a target when the sheet says anything
    # about it: a concentration, or a detection or quantification flag.
    table = table[table[["value", "det", "quant"]].notna().any(axis=1)]
    table = table.merge(people[["ID", *person_columns]], on="ID", how="left")
    note = (
        "Supplement xlsx, sheets Shedding Data and Demographics, joined on ID. One "
        "row per sample and target for the columns N_conc (gc/mg-dw), ORF1a_conc "
        "(gc/gm-dw), PMMoV_conc (gc/mg-dw) and crAss_conc; the column name is kept "
        "in analyte. ID and Day as PatientID and time. The sheet's own flags are "
        "carried beside each value: det holds N_det, ORF1a_det, PMMoV_det or "
        "crAss_det, and quant holds crAss_quant. A row is kept when the sheet "
        "gives a concentration or a flag for that target, so value can be empty "
        "where a target was reported as not detected without a concentration. "
        "solid_frac, extraction_mass_mg, Sex, Age, Symptomatic, Index (I) or "
        "Household member (H) and Vaccination Status kept. The sheet's legend "
        "describes Day as day after acute illness onset, as the authors supplied it."
    )
    return lay_out(table, "ID", "Day", "value"), note


def ke2022daily(release):
    samples = pd.read_excel(
        raw(release, "ke2022daily", "41564_2022_1105_MOESM4_ESM.xlsx"),
        sheet_name="data_samples",
        dtype=str,
    )
    onset = read_csv(raw(release, "ke2022daily", "symptom_date_gap.csv"))
    keep = ["Ind", "Time", "Lineage", "Antigen", "Virus_pos_days", "Age", "Index"]
    parts = []
    for column in ("Nasal_CN", "Saliva_Ct"):
        part = samples[keep].copy()
        part["value"] = samples[column]
        part["analyte"] = column
        parts.append(part)
    table = interleave(parts)
    table = table[table["value"].notna()]
    table = table.merge(onset, on="Ind", how="left")
    note = (
        "Supplement xlsx, sheet data_samples: one row per non-empty cell of the "
        "columns Nasal_CN and Saliva_Ct; the column name is kept in analyte. Ind "
        "and Time as PatientID and time; Lineage, Antigen, Virus_pos_days, Age "
        "and Index kept. first_symptom is joined on Ind from symptom_date_gap.csv, "
        "a file the curator made by reading the first symptom day from Extended "
        "Data Fig. 2 of the article; it is a curator-derived field. Time is the "
        "authors' day index and is not shifted by first_symptom."
    )
    return lay_out(table, "Ind", "Time", "value"), note


def kissler2021densely(release):
    table = read_csv(raw(release, "kissler2021densely", "ct_dat_refined.csv"))
    note = (
        "ct_dat_refined.csv from the authors' repository: one row per row of the "
        "file. PersonID, TestDateIndex and CtT1 as PatientID, time and value; "
        "NovelPersistent, B117Status and PersonIDClean kept. TestDateIndex is the "
        "authors' fractional day index, kept unrounded."
    )
    return lay_out(table, "PersonID", "TestDateIndex", "CtT1"), note


def lavezzo2020suppression(release):
    workbook = raw(
        release, "lavezzo2020suppression", "anonymised_data_public_final.xlsx"
    )
    pcr = pd.read_excel(workbook, sheet_name="RT_PCR_DATA", dtype=str)
    people = pd.read_excel(workbook, sheet_name="anonymised_dataset", dtype=str)
    # The daily swab results sit under date headers; write them as dates.
    dates = {
        c: c.strftime("%Y-%m-%d") for c in people.columns if hasattr(c, "strftime")
    }
    people = people.rename(columns=dates)
    person_columns = [
        "age_group",
        "gender",
        "first_symptoms_date",
        "first_sampling",
        "second_sampling",
        *dates.values(),
    ]
    recorded = {"First_survey": "first_sampling", "Second_survey": "second_sampling"}
    joined = pcr.merge(people[["id", *person_columns]], on="id", how="left")
    parts = []
    for survey, sampling in recorded.items():
        genes = [f"RT_PCR_Genome_Equivalents_{survey}_{gene}" for gene in ("RdRp", "E")]
        # A survey contributes rows for a person when the workbook records a
        # result for that survey, or a quantity for either gene.
        surveyed = joined[sampling].notna() | joined[genes].notna().any(axis=1)
        for column in genes:
            part = joined[
                ["id", "symptomatic_at_first_sampling", "symptomatic_at_follow_up"]
            ].copy()
            part["time"] = survey
            part["value"] = joined[column]
            part["analyte"] = column
            part[person_columns] = joined[person_columns]
            parts.append(part[surveyed])
    table = pd.concat(parts)
    order = {name: i for i, name in enumerate(pcr["id"])}
    table["_id"] = table["id"].map(order)
    table["_survey"] = table["time"].map({s: i for i, s in enumerate(recorded)})
    table = table.sort_values(["_id", "_survey"], kind="stable")
    table = table.drop(columns=["_id", "_survey"])
    note = (
        "Authors' workbook, sheets RT_PCR_DATA and anonymised_dataset, joined on "
        "id. One row per person, survey and gene for the four columns "
        "RT_PCR_Genome_Equivalents_First_survey_RdRp, ..._First_survey_E, "
        "..._Second_survey_RdRp and ..._Second_survey_E; the column name is kept "
        "in analyte. The workbook gives no date for a quantity, only the survey "
        "it belongs to, so time holds the survey name from the column header "
        "(First_survey or Second_survey), not a number. The daily swab results "
        "(one column per date, 2020-02-21 to 2020-03-10, Pos or Neg), "
        "first_sampling, second_sampling and first_symptoms_date are carried so "
        "that a date can be worked out; no date is chosen here. A row is kept "
        "when first_sampling or second_sampling records a result for that survey "
        "or either gene has a quantity, so value is empty where a swab was taken "
        "and nothing was quantified. The _bis columns and the Ct columns are not "
        "carried. age_group, gender, symptomatic_at_first_sampling and "
        "symptomatic_at_follow_up kept."
    )
    return lay_out(table, "id", "time", "value"), note


def liu2024longitudinal(release):
    table = read_csv(raw(release, "liu2024longitudinal", "Liu2024.csv"))
    keep = [
        "subject",
        "day_actual",
        "Age",
        "Gender",
        "Race",
        "Hispanic.or.Latin.Origin",
        "Cohort",
        "inpatient",
        "COVID.19.Confirmed.date",
    ]
    parts = []
    for target in ("mtDNA", "PMMoV", "N1"):
        part = table[keep].copy()
        part["value"] = table[f"gc_dryg_{target}"]
        part["analyte"] = f"gc_dryg_{target}"
        part["dpcr_result_class"] = table[f"dpcr_result_class_{target}"]
        parts.append(part)
    long = interleave(parts)
    long = long[long["value"].notna()]
    note = (
        "Liu2024.csv, contributed by the authors: one row per sample and target "
        "for the columns gc_dryg_mtDNA, gc_dryg_PMMoV and gc_dryg_N1; the column "
        "name is kept in analyte. subject and day_actual as PatientID and time. "
        "dpcr_result_class holds the file's dpcr_result_class_mtDNA, _PMMoV or "
        "_N1 entry for that target, unchanged. Age, Gender, Race, "
        "Hispanic.or.Latin.Origin, Cohort, inpatient and COVID.19.Confirmed.date "
        "kept. day_actual is a day index the authors supplied."
    )
    return lay_out(long, "subject", "day_actual", "value"), note


def natarajan2022gastrointestinal(release):
    study = "natarajan2022gastrointestinal"
    base = (
        "https://github.com/alex-dahlen/lambda_fecal_shedding/raw/7affa71/"
        "Source%20Data/Source%20data_modified.xlsx%20-%20"
    )

    def source(name):
        return download(study, f"{name}.csv", base + name.replace(" ", "%20") + ".csv")

    people = read_csv(source("Baseline data"))
    index = read_csv(source("Index of clinical stool samples"), skiprows=3)
    index = index[
        [
            "RNA sample ID",
            "Subject study ID",
            "Date of sample collection",
            "Date of subject enrollment",
            "Timepoint category*",
            "Stool preservative",
        ]
    ]
    concentration = "Viral RNA concentration (copies/μL)"
    parts = []
    for sheet in (
        "Raw RT_qPCR gRNA data",
        "Raw ddPCR gRNA data",
        "Raw RT-qPCR sgRNA data",
    ):
        assay = read_csv(source(sheet), skiprows=3, encoding="utf-8")
        part = assay[["RNA sample ID"]].copy()
        part["value"] = assay[concentration]
        part["analyte"] = sheet
        part["Target gene"] = assay["Target gene"] if "Target gene" in assay else pd.NA
        part["cQ"] = assay["cQ"] if "cQ" in assay else pd.NA
        parts.append(part)
    table = pd.concat(parts, ignore_index=True)
    # Only wells that hold a participant's stool sample: controls and standards
    # have no entry in the sample index.
    table = table.merge(index, on="RNA sample ID", how="inner")
    table = table.merge(
        people[["Participant ID", "Sex", "Age", "Arm of study"]],
        left_on="Subject study ID",
        right_on="Participant ID",
        how="left",
    ).drop(columns="Participant ID")
    table = table[table["value"].notna()]
    note = (
        "No raw file is stored with this dataset. The builder downloads the five "
        "CSV files the extraction script reads from the authors' repository "
        "(github.com/alex-dahlen/lambda_fecal_shedding, commit 7affa71, folder "
        "Source Data: Baseline data, Index of clinical stool samples, Raw RT_qPCR "
        "gRNA data, Raw ddPCR gRNA data, Raw RT-qPCR sgRNA data) and caches them "
        "under benchmark/manual_inputs/_downloads. One row per assay well whose "
        "RNA sample ID is in the sample index; wells for controls and standards "
        "are left out. Subject study ID, Date of sample collection and Viral RNA "
        "concentration (copies/μL) as PatientID, time and value; time is a "
        "calendar date, and Date of subject enrollment is carried beside it. The "
        "assay table's name is kept in analyte, with Target gene, cQ, Stool "
        "preservative, RNA sample ID and Timepoint category* under their raw "
        "names; the sgRNA table has no Target gene column. Sex, Age and Arm of "
        "study are joined from Baseline data on the participant ID."
    )
    table = table[
        [
            "Subject study ID",
            "Date of sample collection",
            "value",
            "analyte",
            "Target gene",
            "cQ",
            "Stool preservative",
            "RNA sample ID",
            "Date of subject enrollment",
            "Timepoint category*",
            "Sex",
            "Age",
            "Arm of study",
        ]
    ]
    return (
        lay_out(table, "Subject study ID", "Date of sample collection", "value"),
        note,
    )


def tsang2016individual(release):
    table = pd.read_csv(
        raw(release, "tsang2016individual", "data.csv"),
        header=None,
        dtype=str,
        keep_default_na=False,
    )
    # The file has no header. Names follow README.txt, which numbers columns
    # from 1; positions here are from 0.
    named = {
        0: "Household ID",
        1: "Member ID",
        3: "PCR-confirmed infection",
        4: "Date of symptom onset",
        7: "Age",
        8: "Sex",
        9: "Vaccination",
        10: "Date of symptom onset of the index case",
        31: "Influenza virus subtype",
    }
    days = {15 + day: str(day) for day in range(13)}
    table = table.apply(lambda column: column.str.strip())
    table = table[table[3] == "1"]
    index_households = set(table.loc[table[1] == "0", 0])
    table = table[table[0].isin(index_households)]
    people = table[list(named)].rename(columns=named)
    people.insert(0, "PatientID", people["Household ID"] + "-" + people["Member ID"])
    parts = []
    for position, day in days.items():
        part = people.copy()
        part["time"] = day
        part["value"] = table[position]
        parts.append(part)
    long = interleave(parts)
    # README.txt: "in all data file, -1 represents missing".
    long = long[long["value"] != "-1"]
    note = (
        "data.csv from the authors' Dryad deposit, which has no header row; "
        "column names are taken from README.txt. Rows are the people with "
        "PCR-confirmed infection equal to 1 in households whose index case is "
        "also confirmed. One row per person and day for README columns 16 to 28, "
        "'Observed viral shedding from Day 0 to Day 12'; time is that day number, "
        "which README defines as days from symptom onset in the index case, so "
        "it is not each contact's own onset. Cells equal to -1 are left out "
        "because README states that -1 represents missing. PatientID joins "
        "Household ID and Member ID with a hyphen, since neither identifies a "
        "person alone; both are also kept. PCR-confirmed infection, Date of "
        "symptom onset, Date of symptom onset of the index case, Age, Sex, "
        "Vaccination and Influenza virus subtype kept with the file's codes."
    )
    return lay_out(long, "PatientID", "time", "value"), note


def yuan2021sars(release):
    table = read_csv(raw(release, "yuan2021sars", "pone.0247367.s001.csv"))
    # The file writes NA where a row belongs to no case.
    table = table[table["case"].notna() & (table["case"] != "NA")]
    keep = [
        "id",
        "date_collection",
        "specimen",
        "type",
        "age",
        "gender1",
        "date_onset",
        "date_detection",
    ]
    parts = []
    for column in ("ORF1ab", "N"):
        part = table[keep].copy()
        part["value"] = table[column]
        part["analyte"] = column
        parts.append(part)
    long = interleave(parts)
    long = long[long["value"].notna()]
    note = (
        "Supplement csv (S1 file): rows with a case label in the case column (the "
        "file writes NA elsewhere). One row "
        "per sample and target for the columns ORF1ab and N; the column name is "
        "kept in analyte. id and date_collection as PatientID and time; time is a "
        "calendar date, and date_onset and date_detection are carried beside it. "
        "specimen, type, age and gender1 kept. date_onset is NA in the file for "
        "the asymptomatic cases and is carried as written."
    )
    return lay_out(long, "id", "date_collection", "value"), note


BUILDERS = {
    f.__name__: f
    for f in (
        arts2023longitudinal,
        ke2022daily,
        kissler2021densely,
        lavezzo2020suppression,
        liu2024longitudinal,
        natarajan2022gastrointestinal,
        tsang2016individual,
        yuan2021sars,
    )
}
