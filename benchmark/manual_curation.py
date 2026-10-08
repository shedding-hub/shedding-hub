"""
The two steps a production curator took before handing a CSV to the agent.

The production input CSVs of the AI era were not raw. Of the 100 that have a
`value` column, 77 already write non-detects as `negative`, and none of those
that hold concentrations keeps them on a log10 scale. Deciding what counts as
a non-detect, and undoing a logarithm, were the curator's work, not the
agent's. The manual-era CSVs built here follow the same practice, so that the
benchmark asks of the agent what production asked of it.

For each study the rule is the one in the curator's extraction script, applied
to the laid-out table by code:

- non-detects: the cells the script turns into `negative` (or `positive`,
  `inconclusive`) are written that way;
- logarithms: where the script computes `10 ** value`, so does this module.

Nothing else is taken from the scripts. Times are not rounded or re-aligned,
units are not converted, standard curves from cycle threshold to concentration
are not applied, and no label is renamed. Where a script goes on to do one of
those, the note for the study says so.

Each rule returns the table and a sentence for the study's note.
"""

import pandas as pd


def numbers(table: pd.DataFrame) -> pd.Series:
    """The value column as numbers; cells that are not numbers become NaN."""
    return pd.to_numeric(table["value"], errors="coerce")


def blank(table: pd.DataFrame) -> pd.Series:
    return table["value"].isna() | (table["value"].astype(str).str.strip() == "")


def write(table: pd.DataFrame, mask, word: str) -> int:
    table.loc[mask, "value"] = word
    return int(mask.sum())


def to_linear(table: pd.DataFrame, mask) -> int:
    """10 ** value for the chosen cells, written as Python writes a float."""
    table.loc[mask, "value"] = [repr(10 ** float(v)) for v in table.loc[mask, "value"]]
    return int(mask.sum())


def coded(count: int, rule: str) -> str:
    return f"Non-detects: {count} cells written as negative ({rule}), as the extraction script does."


def at_ceiling(threshold: float, exact: bool = False):
    """Cycle thresholds at (or at and above) a ceiling are non-detects."""

    def rule(table):
        value = numbers(table)
        mask = (value == threshold) if exact else (value >= threshold)
        sign = "equal to" if exact else "at or above"
        return table, coded(
            write(table, mask, "negative"), f"value {sign} {threshold:g}"
        )

    return rule


def combined_flag(table):
    mask = pd.to_numeric(table["CombinedDataset_value"], errors="coerce") == 1
    return table, coded(
        write(table, mask, "negative"),
        "the combined dataset's own value column equal to 1, which its legend "
        "defines as below the limit of detection",
    )


def all_log10(table):
    count = to_linear(table, numbers(table).notna())
    return (
        table,
        f"Logarithms: all {count} values converted from log10 to the linear scale, as the extraction script does.",
    )


def not_a_number(table):
    mask = numbers(table).isna()
    return table, coded(write(table, mask, "negative"), "cells that are not numbers")


def alsharrah2020clinical(table):
    value = numbers(table)
    mask = value.isna() | (value == 41)
    return table, coded(
        write(table, mask, "negative"), "value equal to 41, or not a number"
    )


def arts2023longitudinal(table):
    det, quant = table["det"] == "True", table["quant"] == "True"
    crass = table["analyte"] == "crAss_conc"
    positive = write(table, crass & det & ~quant, "positive")
    negative = write(table, (crass & ~det & ~quant) | (~crass & ~det), "negative")
    return table, (
        f"Non-detects: {negative} cells written as negative (det is False) and "
        f"{positive} crAss_conc cells written as positive (detected but not "
        "quantified: det True, quant False), as the extraction script does. The "
        "script keeps the concentration of a non-detect as a limit of blank; "
        "that is not carried."
    )


def covid2020clinical(table):
    value = numbers(table)
    rounded = value.map(lambda v: round(v) if pd.notna(v) else v)
    inconclusive = write(table, (rounded > 40) & (rounded < 50), "inconclusive")
    negative = write(table, value.isna() | (rounded >= 50), "negative")
    return table, (
        f"Non-detects: {negative} cells written as negative (not a number, or a "
        f"reading that rounds to 50 or more) and {inconclusive} as inconclusive "
        "(a reading that rounds to between 41 and 49), as the extraction script "
        "does: readings above 40 mark the figure's bands, not Ct values."
    )


def fajnzylber2020sars(table):
    value = numbers(table)
    negative = write(table, value == 1, "negative")
    linear = to_linear(table, value.notna() & (value != 1))
    return table, (
        f"Non-detects: {negative} cells written as negative (log10 value equal "
        f"to 1). Logarithms: the other {linear} values converted from log10 to "
        "the linear scale. Both as the extraction script does."
    )


def hakki2022onset(table):
    return table, coded(
        write(table, numbers(table) == 1, "negative"), "copy or pfu equal to 1"
    )


def han2020sequential(table):
    limit = 5.7e3  # limit_of_detection in the extraction script
    linear = 10 ** numbers(table)
    count = to_linear(table, linear >= limit)
    negative = write(table, linear < limit, "negative")
    return table, (
        f"Logarithms: {count} values converted from log10 to the linear scale. "
        f"Non-detects: {negative} readings below 5.7e3 on that scale written as "
        "negative; 5.7e3 is the limit of detection the extraction script sets. "
        "The script also rounds time to whole days; that is not applied."
    )


def iwakiri2009quantitative(table):
    value = numbers(table)
    mask = value.isna() | (value == 0)
    return table, coded(
        write(table, mask, "negative"), "cells that are not numbers, or zero"
    )


def ke2022daily(table):
    value = numbers(table)
    mask = ((table["analyte"] == "Nasal_CN") & (value == 48)) | (
        (table["analyte"] == "Saliva_Ct") & (value == 47)
    )
    return table, (
        coded(
            write(table, mask, "negative"),
            "Nasal_CN equal to 48 or Saliva_Ct equal to 47",
        )
        + " The script then converts CN and Ct to concentrations with the "
        "article's calibration formulas; that conversion is not applied."
    )


def kimse2020viral(table):
    value = numbers(table)
    digitized = table["source_file"] == "asymptomatic.xlsx"
    rounded = value.map(lambda v: round(v) if pd.notna(v) else v)
    mask = (~digitized & (value >= 40)) | (digitized & (rounded >= 40))
    return table, (
        coded(
            write(table, mask, "negative"),
            "value at or above 40; for the rows of asymptomatic.xlsx, a value "
            "that rounds to 40 or more",
        )
        + " The script also rounds the values of asymptomatic.xlsx to whole "
        "cycles; that is not applied."
    )


def kissler2021densely(table):
    value = numbers(table)
    mask = value.isna() | (value == 40)
    return table, coded(write(table, mask, "negative"), "value equal to 40, or empty")


def kissler2021viral(table):
    value = numbers(table)
    ct = table["analyte"] == "CT.Mean"
    floor = 2.65760968215156  # the file's lowest log10_GEperML, as the script tests
    negative = write(
        table,
        (ct & (value.isna() | (value >= 40)))
        | (~ct & (value.isna() | (value == floor))),
        "negative",
    )
    linear = to_linear(table, ~ct & value.notna() & (value != floor))
    return table, (
        f"Non-detects: {negative} cells written as negative (CT.Mean at or above "
        f"40; log10_GEperML equal to {floor}). Logarithms: the other {linear} "
        "log10_GEperML values converted to the linear scale. Both as the "
        "extraction script does."
    )


def lavezzo2020suppression(table):
    return table, (
        coded(write(table, blank(table), "negative"), "empty cells")
        + " The script also works out a date for each quantity from the daily "
        "swab columns; no date is chosen here."
    )


def lescure2020clinical(table):
    value, time = numbers(table), pd.to_numeric(table["time"], errors="coerce")
    linear = to_linear(table, (table["cens"] == "0") & value.notna())
    negative = write(table, table["cens"] == "1", "negative")
    positive = write(
        table,
        ((table["PatientID"] == "E3") & time.isin([13, 15]))
        | ((table["PatientID"] == "E5") & (time == 8)),
        "positive",
    )
    return table, (
        f"Logarithms: {linear} values with cens 0 converted from log10 to the "
        f"linear scale. Non-detects: {negative} cells with cens 1 written as "
        f"negative, then {positive} of them (E3 on days 13 and 15, E5 on day 8) "
        "written as positive, the curator's reading of Fig. 2 for samples that "
        "were positive but not quantifiable. All as the extraction script does. "
        "The script also moves the two E3 samples from days 13 and 15 to days 14 "
        "and 16, a correction from Fig. 2; that is not applied."
    )


def liu2024longitudinal(table):
    mask = table["dpcr_result_class"] == "Negative"
    return table, coded(write(table, mask, "negative"), "dpcr_result_class is Negative")


def lui2020viral(table):
    empty = blank(table)
    linear = to_linear(table, ~empty)
    negative = write(table, empty, "negative")
    return table, (
        f"Non-detects: {negative} readings stored as null written as negative; "
        "the file's notes say these are negative samples read from the figures. "
        f"Logarithms: the other {linear} values converted from log10 to the "
        "linear scale."
    )


def natarajan2022gastrointestinal(table):
    return table, (
        coded(write(table, numbers(table) == 0, "negative"), "concentration equal to 0")
        + " The script then multiplies concentrations by 1,000 to go from "
        "copies per microlitre to copies per millilitre; that is not applied."
    )


def pace2024prevalence(table):
    return table, coded(
        write(table, numbers(table) == 0, "negative"), "Ct Value equal to 0"
    )


def vetter2020daily(table):
    value = numbers(table)
    mask = value.isna() | (value == 0)
    return table, coded(write(table, mask, "negative"), "cells reading neg, or zero")


def woelfel2020virological(table):
    value = numbers(table)
    negative = write(table, value < 0.05, "negative")
    linear = to_linear(table, value >= 0.05)
    return table, (
        f"Non-detects: {negative} readings below 0.05 on the log10 axis written "
        f"as negative. Logarithms: the other {linear} values converted from "
        "log10 to the linear scale. Both as the extraction script does. The "
        "script also rounds days; that is not applied."
    )


def xu2020characteristics(table):
    value = numbers(table)
    count = write(table, (value >= 39) & (value <= 41), "negative")
    return table, (
        coded(count, "a digitized Ct between 39 and 41, which the script reads as 40")
        + " The script then converts Ct to copies per mL with a standard curve it "
        "fits; that conversion is not applied."
    )


def young2020epidemiologic(table):
    table, note = at_ceiling(38, exact=True)(table)
    return table, note + (
        " The script then converts Ct to copies per swab with a standard curve "
        "it fits; that conversion is not applied."
    )


def zuo2020alterations(table):
    value = numbers(table)
    negative = write(table, value == 0, "negative")
    linear = to_linear(table, value.notna() & (value != 0))
    return table, (
        f"Logarithms: {linear} values converted from log10 to the linear scale, "
        f"as the extraction script does. {negative} further cells (log10 equal "
        "to 0) written as negative; the file already writes the other "
        "non-detects as negative."
    )


RULES = {
    "alsharrah2020clinical": alsharrah2020clinical,
    "arts2023longitudinal": arts2023longitudinal,
    "covid2020clinical": covid2020clinical,
    "fajnzylber2020sars": fajnzylber2020sars,
    "gautret2020hydroxychloroquine": combined_flag,
    "hakki2022onset": hakki2022onset,
    "han2020sequential": han2020sequential,
    "iwakiri2009quantitative": iwakiri2009quantitative,
    "ke2022daily": ke2022daily,
    "kim2020viral": not_a_number,  # the builder has already spread `ud`
    "kimse2020viral": kimse2020viral,
    "kissler2021densely": kissler2021densely,
    "kissler2021viral": kissler2021viral,
    "lavezzo2020suppression": lavezzo2020suppression,
    "lescure2020clinical": lescure2020clinical,
    "liu2024longitudinal": liu2024longitudinal,
    "lui2020viral": lui2020viral,
    "natarajan2022gastrointestinal": natarajan2022gastrointestinal,
    "obara2008single": all_log10,
    "pace2024prevalence": pace2024prevalence,
    "peiris2003clinical": all_log10,
    "salvatore2020epidemiological": at_ceiling(40),
    "seah2020assessing": not_a_number,
    "shrestha2020distribution": combined_flag,
    "tan2021early": at_ceiling(38, exact=True),
    "team2020clinical": combined_flag,
    "tsang2016individual": all_log10,
    "vetter2020daily": vetter2020daily,
    "wang2020fecal": at_ceiling(40),
    "woelfel2020virological": woelfel2020virological,
    "xing2020prolonged": at_ceiling(40),
    "xu2020characteristics": xu2020characteristics,
    "yang2020laboratory": not_a_number,
    "yilmaz2020upper": at_ceiling(40, exact=True),
    "young2020epidemiologic": young2020epidemiologic,
    "yuan2021sars": at_ceiling(40),
    "zuo2020alterations": zuo2020alterations,
}


def curate(study: str, table: pd.DataFrame):
    """Apply a study's rule. Returns the table and the sentence for its note."""
    table = table.copy()
    table["value"] = table["value"].astype(object)
    return RULES[study](table)
