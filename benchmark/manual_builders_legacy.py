"""
Input CSV builders for two studies with no digitized file in this repository.

woelfel2020virological was digitized by code, from the vector graphics of the
article PDF; its builder repeats that reading step and stops there.

lui2020viral's readings survive only in the earlier SheddingHub repository (github.com/CIDMATH/SheddingHub, commit 31db6f6,
"Uploaded initial files from Hoffmann et al. 2023"), not in this one.

Those files are JSON datasets with one entry per sample under `loads`. A copy
of each file used is kept in `benchmark/manual_inputs/_sources/`. Only the
`loads` entries are laid out. The file's own descriptive fields (assay, unit,
limit of quantification, reference event, severity) are decisions about the
study and are not carried.
"""

import json
from bisect import bisect

import numpy as np
import pandas as pd

from manual_inputs import OUT, lay_out, raw


def lui2020viral(release):
    source = OUT / "_sources" / "Lui2020Viral.json"
    loads = json.loads(source.read_text(encoding="utf-8"))["loads"]
    table = pd.DataFrame(
        {
            "patient": [str(entry["patient"]) for entry in loads],
            "day": [str(entry["day"]) for entry in loads],
            # A reading the file stores as null stays an empty cell.
            "value": [
                "" if entry["value"] is None else repr(entry["value"])
                for entry in loads
            ],
        }
    )
    note = (
        "Lui2020Viral.json from the earlier SheddingHub repository (CIDMATH, "
        "commit 31db6f6), the hand-made file this dataset was first curated in: "
        "one row per entry of loads, with patient, day and value as PatientID, "
        "time and value. Values are as stored there (the file states log10 "
        "gc/mL) and day is as stored there (the file states days since "
        "symptoms). Readings stored as null are left empty; the file's notes say "
        "these are negative samples read from the figures. The file's assay, "
        "unit, limit of quantification, reference event and severity fields are "
        "not carried."
    )
    return lay_out(table, "patient", "day", "value"), note


# woelfel2020virological: digitization parameters, as the extraction script
# sets them. The figure is a 3 x 3 grid of patient panels on page 3 of the PDF.
WOELFEL_XDIVS = [225, 395]
WOELFEL_YDIVS = [-340, -240, -145]
WOELFEL_SERIES = {
    (0.9289997816085815, 0.4899977147579193, 0.1919890195131302): "sputum",
    (1.0, 0.7529869675636292, 0.0): "oropharyngeal_swab",
    (0.49799343943595886, 0.49799343943595886, 0.49799343943595886): "stool",
}
WOELFEL_XLIMS = [
    (3.5, 22.5),
    (2.5, 20.5),
    (2.5, 23.5),
    (3.5, 20.5),
    (3.5, 27.5),
    (5.5, 22.5),
    (3.5, 28.5),
    (1.5, 15.5),
    (7.5, 12.5),
]
WOELFEL_YLIM = (0, 10)


def woelfel2020virological(release):
    """Digitize Fig. 2 from the PDF's vector graphics, in axis units.

    This is the reading step of the extraction script and nothing after it:
    marker positions are taken from the page, assigned to a patient panel and
    a series by position and colour, and mapped to the axes. The script's
    rounding of days, its non-detect threshold and its conversion from log10
    are not applied.
    """
    import pymupdf

    page = pymupdf.Document(
        stream=raw(
            release, "woelfel2020virological", "woelfel2020virological.pdf"
        ).read()
    )[2]

    def panel(x, y):
        col = bisect(WOELFEL_XDIVS, x)
        row = 3 - bisect(WOELFEL_YDIVS, y)
        return None if row > 2 else col + 3 * row

    axes, markers = {}, {}
    for drawing in page.get_drawings():
        rect, color = drawing["rect"], drawing["color"]
        if not color:
            continue
        if (not rect.width or not rect.height) and color == (0, 0, 0):
            continue
        if all(item[0] == "c" for item in drawing["items"]):
            if color == (0, 0, 0):
                continue
            x, y = (rect.x0 + rect.x1) / 2, -(rect.y0 + rect.y1) / 2
            markers.setdefault(panel(x, y), []).append((color, x, y))
            continue
        if color != (0, 0, 0):
            continue
        for _, *points in drawing["items"]:
            xs, ys = [pt.x for pt in points], [-pt.y for pt in points]
            index = panel(xs[0], ys[0])
            if index is not None:
                axes.setdefault(index, []).append((xs, ys))

    rows = []
    for patient in sorted(axes):
        xaxis = max(axes[patient], key=lambda xy: abs(xy[0][0] - xy[0][1]))
        yaxis = max(axes[patient], key=lambda xy: abs(xy[1][0] - xy[1][1]))
        xmap = np.polynomial.Polynomial.fit(xaxis[0], WOELFEL_XLIMS[patient], 1)
        ymap = np.polynomial.Polynomial.fit(yaxis[1], WOELFEL_YLIM, 1)
        for color, x, y in markers[patient]:
            rows.append(
                {
                    "panel": str(patient),
                    "day": repr(float(xmap(x))),
                    "log10": repr(float(ymap(y))),
                    "specimen": WOELFEL_SERIES[color],
                }
            )
    table = pd.DataFrame(rows)
    note = (
        "woelfel2020virological.pdf: marker positions read from the vector "
        "graphics of Fig. 2 with the digitization settings of the extraction "
        "script (panel boundaries, axis limits, and the series each marker "
        "colour stands for), in axis units. PatientID is the panel index, time "
        "is the x-axis reading and value is the y-axis reading; neither is "
        "rounded or converted. The series name is in specimen. The script's "
        "rounding of days, its non-detect threshold and its conversion from "
        "log10 are not applied, and the ages and sexes it types in are not "
        "carried."
    )
    return lay_out(table, "panel", "day", "log10"), note


BUILDERS = {
    "lui2020viral": lui2020viral,
    "woelfel2020virological": woelfel2020virological,
}
