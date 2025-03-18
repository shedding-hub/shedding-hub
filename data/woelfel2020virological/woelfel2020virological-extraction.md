# Extracting Data From Vector Graphics

When raw data are not available, they can often be extracted from figures in publications. Vector graphics, as are often present in modern PDFs, are particularly suitable for data extraction because each visual element can be accessed programmatically. For rasterized figures, other methods, such as [WebPlotDigitizer](https://automeris.io), may be suitable. Here, we extract viral loads for throat swabs and stool and sputum samples from [one of the early virological analyses conducted to study SARS-CoV-2](https://doi.org/10.1038/s41586-020-2196-x). Unfortunately, there is no "standard" format for these vector graphics, and each extraction is likely to take some time to get right. Throat samples were recorded as oropharyngeal samples in the standardized dataset.

The data for nine patients are shown in figure 2 on the third page as a vector graphic. We first download the PDF document, load it using the [`pymupdf`](https://pymupdf.readthedocs.io/en/latest/index.html) library, and extract the third page.

```python
from bisect import bisect
import io
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import pymupdf
import yaml


# Load the document from Nature and get the page with the figure.
with open("woelfel2020virological.pdf", "rb") as fp:
    document = pymupdf.Document(stream=fp.read())
page = document[2]
```

Having obtained the page, we can get all drawings on the page using the [`get_drawings`](https://pymupdf.readthedocs.io/en/latest/page.html#Page.get_drawings) method of the page (see [here](https://pymupdf.readthedocs.io/en/latest/recipes-drawing-and-graphics.html#how-to-extract-drawings) for further details). Each drawing is a vector graphics "path", such as a line or polygon (see [here](https://developer.mozilla.org/en-US/docs/Web/SVG/Tutorial/Paths) for a discussion of paths in the context of SVG files). Before extracting data from the figure, let us first draw all lines and markers in the `drawings` object to verify we have obtained the right information. We also return all lines and markers grouped by their color.

```python
def draw_lines_and_markers(drawings, alpha=0.7, ax=None):
    """
    Draw all lines in the figure that have a color specified.

    Args:
        drawings: Drawings to plot.
        ax: Axes to plot into or the default axes.

    Returns:
        Tuple comprising lines and markers each keyed by color.
    """
    ax = ax or plt.gca()
    ax.set_aspect("equal")

    lines = {}
    markers = {}

    for drawing in drawings:
        rect = drawing["rect"]

        # Skip if there is no color.
        color = drawing["color"]
        if not color:
            continue

        # Skip if the drawing doesnt have a size and the color is black.
        if (not rect.width or not rect.height) and color == (0, 0, 0):
            continue

        # If all of the kinds are "c", we likely just have marker. If the marker is
        # black, we don't care about it and drop it.
        is_marker = all(item[0] == "c" for item in drawing["items"])
        if is_marker:
            if color == (0, 0, 0):
                continue
            x = (rect.x0 + rect.x1) / 2
            y = - (rect.y0 + rect.y1) / 2
            ax.scatter(x, y, color=color, marker=".", alpha=alpha)
            markers.setdefault(color, []).append((x, y))
            continue

        for kind, *points in drawing["items"]:
            # Check we only have lines if the drawing wasn't a marker.
            assert kind == "l"
            # Get the line coordinates. We invert the y coordinate because PDF
            # documents use the convention to measure distance from the top instead
            # of the bottom of the page.
            x = [pt.x for pt in points]
            y = [-pt.y for pt in points]

            # Add the line to the plot and our collection of lines.
            marker = None if color == (0, 0, 0) else "x"
            ax.plot(x, y, color=color, alpha=alpha, marker=marker)
            lines.setdefault(color, []).append((x, y))

    return lines, markers


# Get the drawings and plot them.
fig, ax = plt.subplots()
drawings = page.get_drawings()
lines_by_color, markers_by_color = draw_lines_and_markers(drawings)
fig.tight_layout()
```

We are able to successfully reproduce the figure in Python and are almost ready to extract the data. However, we need to separate the data by patient. Fortunately, the panels nicely fall into a grid, and we can use divisions along the x and y axes to separate lines capturing data for different patients.

```python
fig, ax = plt.subplots()
ax.set_aspect("equal")

# X and Y divisions between different patients. We show the divisions as dotted lines.
xdivs = [225, 395]
for xdiv in xdivs:
    ax.axvline(xdiv, color="magenta")
ydivs = [-340, -240, -145]
for ydiv in ydivs:
    ax.axhline(ydiv, color="magenta")


def get_patient_idx(x, y):
    # Look up the x and y indices for any vertex of the line using the bisect method
    # (see https://docs.python.org/3/library/bisect.html for details).
    col = bisect(xdivs, x)
    row = 3 - bisect(ydivs, y)
    if row > 2:
        return None
    return col + 3 * row


# Run over all lines and group them by patient.
lines_by_patient_and_color = {}
for color, lines in lines_by_color.items():
    for line in lines:
        x, y = line
        idx = get_patient_idx(x[0], y[0])
        if idx is None:
            continue
        # Save the line by patient and color.
        lines_by_patient_and_color.setdefault(idx, {}).setdefault(color, []).append(
            line
        )

        # Demonstrate that we got the line right.
        ax.plot(x, y, color=f"C{idx}")

# Do the same for markers.
markers_by_patient_and_color = {}
for color, markers in markers_by_color.items():
    for marker in markers:
        idx = get_patient_idx(*marker)
        markers_by_patient_and_color.setdefault(idx, {}).setdefault(color, []).append(
            marker
        )

        # Demonstrate that we got the marker right.
        ax.scatter(*marker, color=f"C{idx}", marker=".")
```

We now have all the data grouped by patient and color of the lines which allows us to distinguish stool from sputum samples and throat swaps from the lines that represent the axes. The colors are however cumbersome, and we replace them by semantically meaningful keys below.

```python
kind_by_color = {
    (0, 0, 0): "axes",
    (0.9289997816085815, 0.4899977147579193, 0.1919890195131302): "sputum",
    (1.0, 0.7529869675636292, 0.0): "oropharyngeal_swab",
    (0.49799343943595886, 0.49799343943595886, 0.49799343943595886): "stool",
}

lines_by_patient_and_kind = {
    patient: {kind_by_color[color]: lines for color, lines in by_color.items()}
    for patient, by_color in sorted(lines_by_patient_and_color.items())
}
markers_by_patient_and_kind = {
    patient: {kind_by_color[color]: markers for color, markers in by_color.items()}
    for patient, by_color in sorted(markers_by_patient_and_color.items())
}
```

We are finally ready to extract the data for each patient. We obtain the x and y axes by picking the longest horizontal and vertical lines of the `axes` kind, respectively. This allows us to transform from the graphics coordinates to the data coordinates. The y axis has limits from 0 (representing a negative sample) to 10 on the log10 scale. But the x axes of each subplot are not shared, and we manually extract the limits from the figure in the publication. The attributes are manually defined.

```python
fig, ax = plt.subplots()

attributes = [
    {"age": 34, "sex": "male"},
    {"age": 41, "sex": "male"},
    {"age": 28, "sex": "male"},
    {"age": 33, "sex": "male"},
    {"age": 52, "sex": "male"},
    {"age": 33, "sex": "male"},
    {"age": 58, "sex": "male"},
    {"age": 49, "sex": "male"},
    {"age": 49, "sex": "male"},
]

xlims = [
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

patients = []
for patient, lines_by_kind in lines_by_patient_and_kind.items():
    # Extract the axes.
    lines = lines_by_kind["axes"]
    xaxis = max(lines, key=lambda xy: abs(xy[0][0] - xy[0][1]))
    yaxis = max(lines, key=lambda xy: abs(xy[1][0] - xy[1][1]))
    ax.plot(*xaxis, color="C0")
    ax.plot(*yaxis, color="C1")

    # Linear fit from the graphics coordinates to the data coordinates.
    xinterp = np.polynomial.Polynomial.fit(xaxis[0], xlims[patient], 1)
    yinterp = np.polynomial.Polynomial.fit(yaxis[1], (0, 10), 1)

    measurements = []
    for kind, markers in markers_by_patient_and_kind[patient].items():
        x, y = np.transpose(markers)
        ax.plot(x, y, marker=".")

        # Transform from graphics to data coordinates.
        days = xinterp(x)
        log10s = yinterp(y)

        # Sanity check that the days are close to integers.
        np.testing.assert_allclose(days, np.round(days), atol=0.05)

        # Store data in the format expected by the schema.
        measurements.extend({
            "analyte": kind,
            "time": round(day),
            "value": "negative" if log10 < 0.05 else float(10 ** log10)
        } for day, log10 in zip(days, log10s))
    patients.append({"attributes": attributes[patient], "measurements": measurements})
```

Having extracted the data, we can now plot it on regular axes to verify that the extraction worked as expected.

```python
fig, axes = plt.subplots(3, 3, sharex=True, sharey=True)
for i, (ax, patient) in enumerate(zip(axes.ravel(), patients)):
    # For plotting, we need to group the samples by analyte again.
    xy_by_analyte = {}
    for measurement in patient["measurements"]:
        xy_by_analyte.setdefault(measurement["analyte"], []).append(
            (
                measurement["time"],
                0.1 if measurement["value"] == "negative" else measurement["value"],
            )
        )
    for analyte, xy in sorted(xy_by_analyte.items()):
        ax.plot(*np.transpose(xy), label=analyte, marker=".")
    ax.set_yscale("log")
    ax.axhline(1, color="k")

    print(f"patient {i}", {key: len(value) for key, value in xy_by_analyte.items()})

ax.legend(fontsize="small")
fig.tight_layout()
```

Finally, we write the data as a YAML file which can be copied into the `data` folder of the repository.

```python
# Dump the data to YAML in the correct format.
with open("woelfel2020virological.yaml", "w") as fp:
    fp.write("# yaml-language-server: $schema=../.schema.yaml\n")
    yaml.dump({
        "title": "Virological assessment of hospitalized patients with COVID-2019",
        "doi": "10.1038/s41586-020-2196-x",
        "description": "The authors conducted a virological analysis of nine linked "
        "cases of COVID-19 in Munich in early 2020. They quantified SARS-CoV-2 RNA "
        "gene copies in throat swabs and RNA concentrations in stool and sputum "
        "samples. Abundances were quantified using RT-qPCR assays targeting the E and "
        "RdRP genes as described in 10.2807/1560-7917.ES.2020.25.3.2000045. Values "
        "were programmatically extracted from the figure, resulting in a number of "
        "significant figures far exceeding the performance of the assays. The "
        "demographic data, including age and sex, are derived from the Challenger et "
        "al. BMC Medicine (2022) 20:25 https://doi.org/10.1186/s12916-021-02220-0 "
        "supplement combined dataset.",
        "analytes": {
            "stool": {
                "description": "RNA gene copy concentration in stool samples. The "
                "authors report that \"stool samples were taken and shipped in native "
                "conditions,\" suggesting that results reported as gene copies per "
                "gram refer to wet weight.",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
                "specimen": "stool",
                "biomarker": "SARS-CoV-2",
                "gene_target": "E and RdRP (not further specified by authors)",
                "unit": "gc/mL",
                "reference_event": "symptom onset"
            },
            "sputum": {
                "description": "RNA gene copy concentration in sputum samples. Results "
                "are reported as gene copies per mL.",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
                "specimen": "sputum",
                "biomarker": "SARS-CoV-2",
                "gene_target": "E and RdRP (not further specified by authors)",
                "unit": "gc/mL",
                "reference_event": "symptom onset"
            },
            "throat_swab": {
                "description": "Number of gene copies per throat swab.",
                "limit_of_quantification": 100,
                "limit_of_detection": "unknown",
                "specimen": "throat_swab",
                "biomarker": "SARS-CoV-2",
                "gene_target": "E and RdRP (not further specified by authors)",
                "unit": "gc/swab",
                "reference_event": "symptom onset"
                },
        },
        "participants": patients,
    }, fp, sort_keys=False)
```
