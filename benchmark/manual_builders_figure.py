from manual_inputs import lay_out, raw

import json

import pandas as pd

# Builders for the manual-era studies whose values were digitized from figures.
#
# The digitized values are taken, in order of preference, from a digitized file
# stored with the dataset, from a WebPlotDigitizer project file stored with the
# dataset, or from numbers the curator typed into the extraction script as
# literal data. The same limits apply as in `manual_inputs.py`: layout only, no
# value is changed. A study whose digitized values exist nowhere but in the
# output of the extraction script has no builder.


def read_csv(release, study, name):
    """A raw CSV with every cell kept as text and any byte-order mark removed."""
    return pd.read_csv(raw(release, study, name), dtype=str, encoding="utf-8-sig")


def melt_in_row_order(table, id_vars, value_vars, var_name):
    """One row per non-empty cell, in the raw row order, then the column order."""
    table = table.reset_index(drop=True).reset_index(names="_row")
    long = table.melt(
        id_vars=["_row", *id_vars],
        value_vars=value_vars,
        var_name=var_name,
        value_name="value",
    )
    long = long[long["value"].notna()]
    long["_col"] = long[var_name].map({c: i for i, c in enumerate(value_vars)})
    long = long.sort_values(["_row", "_col"], kind="stable")
    return long.drop(columns=["_row", "_col"])


def covid2020clinical(release):
    table = read_csv(release, "covid2020clinical", "data.csv")
    note = (
        "data.csv, digitized from the figures: ID, Day and Value as PatientID, time "
        "and value; Specimen kept. Day and Value are raw axis readings. Readings "
        "above 40 mark the figure's inconclusive and negative bands, not Ct values, "
        "and some non-detects are already written as negative, a curator's coding."
    )
    return lay_out(table, "ID", "Day", "Value"), note


def hakki2022onset(release):
    table = read_csv(release, "hakki2022onset", "trajectories.csv")
    long = melt_in_row_order(
        table,
        ["participant", "day", "days_since_peak", "vaccinated", "WGS"],
        ["copy", "pfu"],
        "analyte",
    )
    note = (
        "trajectories.csv from the authors' repository: one row per non-empty cell "
        "of the columns copy and pfu; the column name is kept in analyte. "
        "participant and day as PatientID and time; days_since_peak, vaccinated and "
        "WGS kept. LFD and the *_exist flags are not carried. days_since_peak is the "
        "authors' own alignment to peak viral load."
    )
    return lay_out(long, "participant", "day", "value"), note


def han2020sequential(release):
    rows = []
    for name in ("mother_project.json", "neonate_project.json"):
        project = json.load(raw(release, "han2020sequential", name))
        for dataset in project["datasetColl"]:
            for point in dataset["data"]:
                time, value = point["value"]
                rows.append(
                    {
                        "PatientID": name.removesuffix(".json"),
                        "time": repr(time),
                        "value": repr(value),
                        "analyte": dataset["name"],
                    }
                )
    note = (
        "WebPlotDigitizer project files mother_project.json and "
        "neonate_project.json, digitized from Fig. 1: one row per digitized point, "
        "with the raw x and y axis values as time and value. PatientID is the "
        "project file name. analyte is the dataset name the curator gave each "
        "series in the project file, which already names a specimen and a gene "
        "target."
    )
    return pd.DataFrame(rows), note


def lescure2020clinical(release):
    table = read_csv(release, "lescure2020clinical", "Viral_Loads.csv")
    table = table[table["cov_study"] == "4"][["ID", "dao", "VL", "cens"]]
    # Typed into the extraction script: the curator's reading of Fig. 2 for the
    # one patient the combined file does not include.
    typed = pd.DataFrame(
        {
            "ID": ["E2", "E2", "E2", "E2", "E2"],
            "dao": ["0", "2", "5", "11", "13"],
            "VL": ["positive", "negative", "negative", "negative", "negative"],
        }
    )
    table = pd.concat([table, typed], ignore_index=True)
    note = (
        "Viral_Loads.csv (Goyal et al., a third-party combined dataset), rows with "
        "cov_study 4: ID, dao and VL as PatientID, time and value; cens kept. Ct "
        "and the regression columns are not carried. dao is days after the first "
        "positive sample, as that dataset defines it. The five rows for patient E2 "
        "are not in the file: they are the curator's reading of Fig. 2, typed into "
        "the extraction script, and already written as positive or negative."
    )
    return lay_out(table, "ID", "dao", "VL"), note


def pace2024prevalence(release):
    table = read_csv(release, "pace2024prevalence", "data for pace2024prevalence.csv")
    note = (
        "data for pace2024prevalence.csv, digitized from Figs. 3 and 4: Patient "
        "Number, Day and Ct Value as PatientID, time and value; sex, infant and "
        "Specimen kept."
    )
    return lay_out(table, "Patient Number", "Day", "Ct Value"), note


def peiris2003clinical(release):
    table = read_csv(release, "peiris2003clinical", "data.csv")
    note = (
        "data.csv, digitized from Fig. 4: PatientID, Day and value as PatientID, "
        "time and value. Day and value are raw axis readings."
    )
    return lay_out(table, "PatientID", "Day", "value"), note


def seah2020assessing(release):
    table = read_csv(release, "seah2020assessing", "data for seah2020assessing.csv")
    note = (
        "data for seah2020assessing.csv, digitized from Fig. 1: Patient Number, "
        "Days.Since.Initial.COVID-19.Symptoms and ct value as PatientID, time and "
        "value. The raw time column is named for days since initial symptoms, a "
        "reference event the curator supplied, and its cells read 'Day N'."
    )
    return (
        lay_out(
            table,
            "Patient Number",
            "Days.Since.Initial.COVID-19.Symptoms",
            "ct value",
        ),
        note,
    )


def vetter2020daily(release):
    table = read_csv(release, "vetter2020daily", "msphere.00827-20-st002.csv")
    note = (
        "msphere.00827-20-st002.csv, the article's supplementary table: patient "
        "number, dpo and Log copies/ml as PatientID, time and value; Date, Sample "
        "and Virus isolation kept. dpo is days post onset, as the authors supplied "
        "it."
    )
    return lay_out(table, "patient number", "dpo", "Log copies/ml"), note


def wang2020fecal(release):
    parts = []
    for name in ("A data.xlsx", "B data.xlsx"):
        sheet = pd.read_excel(raw(release, "wang2020fecal", name), dtype=str)
        sheet["source_file"] = name
        parts.append(sheet)
    table = pd.concat(parts, ignore_index=True)
    note = (
        "A data.xlsx and B data.xlsx, digitized from the two panels of Fig. 2: "
        "patients, days and ctvalue as PatientID, time and value; the file name is "
        "kept in source_file. days and ctvalue are raw axis readings."
    )
    return lay_out(table, "patients", "days", "ctvalue"), note


def xing2020prolonged(release):
    # Typed into the extraction script as the digitized series of Fig. 4. The
    # script pairs ct_throat with every time point and ct_fecal with the time
    # points from the third one on.
    time = [1, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36]
    ct_fecal = [
        12.65,
        12.31,
        15.81,
        13.97,
        18,
        24.37,
        28.29,
        34.67,
        37.05,
        38.76,
        40,
        40,
        40,
        40,
        40,
        40,
        40,
    ]
    ct_throat = [
        18.17,
        23.97,
        20.53,
        25.91,
        28.36,
        34.83,
        32.17,
        38.19,
        40,
        40,
        40,
        40,
        40,
        40,
        40,
        40,
        40,
        40,
        40,
    ]
    rows = [
        {"PatientID": "Case1", "time": t, "value": v, "specimen": "ct_throat"}
        for t, v in zip(time, ct_throat)
    ]
    rows += [
        {"PatientID": "Case1", "time": t, "value": v, "specimen": "ct_fecal"}
        for t, v in zip(time[2:], ct_fecal)
    ]
    note = (
        "No digitized file is stored: the values are the lists time, ct_throat and "
        "ct_fecal typed into the extraction script, digitized from Fig. 4. The list "
        "name is kept in specimen. PatientID is the key the script uses for the one "
        "participant. The pairing of ct_fecal with the time points from the third "
        "one on follows the script."
    )
    return pd.DataFrame(rows).astype(str), note


def xu2020characteristics(release):
    table = read_csv(
        release, "xu2020characteristics", "data for xu2020characteristics.csv"
    )
    note = (
        "data for xu2020characteristics.csv, digitized from Fig. 1b: Patient "
        "Number, Day and Ct Value as PatientID, time and value; Sex, Age and "
        "Specimen kept. Day and Ct Value are raw axis readings. The file does not "
        "state the unit of Age."
    )
    return lay_out(table, "Patient Number", "Day", "Ct Value"), note


def yang2020laboratory(release):
    sheet = pd.read_excel(raw(release, "yang2020laboratory", "data.xlsx"), dtype=str)
    sheet.columns = [str(c) for c in sheet.columns]
    per_case = ["Sex", "Age", "Adm_day", "treat_day"]
    days = [c for c in sheet.columns if c not in ("Case", "Resp_tract", *per_case)]
    # The per-case columns are filled on one of a case's two rows only.
    sheet[per_case] = sheet.groupby("Case")[per_case].transform("first")
    long = melt_in_row_order(sheet, ["Case", "Resp_tract", *per_case], days, "day")
    note = (
        "data.xlsx, digitized from Fig. 2 and Supplementary Fig. 1: one row per "
        "non-empty cell of the day columns 1 to 50; Case and the day column name "
        "as PatientID and time; Resp_tract, Sex, Age, Adm_day and treat_day kept. "
        "Sex, Age, Adm_day and treat_day are filled on one row per case in the "
        "file and are repeated here on every row of that case. Non-detects are "
        "coded u in the file."
    )
    return lay_out(long, "Case", "day", "value"), note


BUILDERS = {
    f.__name__: f
    for f in (
        covid2020clinical,
        hakki2022onset,
        han2020sequential,
        lescure2020clinical,
        pace2024prevalence,
        peiris2003clinical,
        seah2020assessing,
        vetter2020daily,
        wang2020fecal,
        xing2020prolonged,
        xu2020characteristics,
        yang2020laboratory,
    )
}
