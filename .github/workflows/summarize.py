import argparse
import pandas as pd
import pathlib
import subprocess
import yaml


def with_name(x: pd.Series, name: str) -> pd.Series:
    if isinstance(x, dict):
        x = pd.Series(x)
    x.name = name
    return x


def __main__(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "-o", help="output filename")
    parser.add_argument(
        "filenames",
        nargs="*",
        help="input filenames (defaults to files changed compared with the `main` "
        "branch)",
    )
    args = parser.parse_args(argv)

    # Get all changed data files compared with the main branch or the explicit
    # filenames.
    if args.filenames:
        filenames = args.filenames
    else:
        output = subprocess.check_output(
            ["git", "diff", "--name-only", "origin/main"], text=True
        )
        # Select only changed data files.
        filenames = [
            name
            for name in map(pathlib.Path, output.splitlines())
            if name.parent.name == "data"
            and name.suffix == ".yaml"
            and not name.name.startswith(".")
        ]

    # Iterate over files and summarize their content.
    lines = []
    for filename in filenames:
        with open(filename) as fp:
            data = yaml.safe_load(fp)

        measurements = pd.DataFrame(
            [
                measurement | {"participant": i}
                for i, participant in enumerate(data["participants"])
                for measurement in participant["measurements"]
            ]
        )

        # Summarize changes.
        grouped = measurements.groupby("analyte")
        lines.extend(
            [
                f"🔄 Summary for changed file `{filename}`:",
                "",
                "```",
                pd.DataFrame(
                    [
                        with_name(grouped.count().value, "n_samples"),
                        with_name(
                            grouped.nunique().participant, "n_unique_participants"
                        ),
                        with_name(
                            {
                                key: (subset.value == "negative").sum()
                                for key, subset in grouped
                            },
                            "n_negative",
                        ),
                        with_name(
                            {
                                key: (subset.value == "positive").sum()
                                for key, subset in grouped
                            },
                            "n_positive",
                        ),
                        with_name(
                            {
                                key: (
                                    ~subset.value.isin({"positive", "negative"})
                                ).sum()
                                for key, subset in grouped
                            },
                            "n_quantified",
                        ),
                    ]
                ).T.to_string(),
                "```",
                "",
            ]
        )

    # Either print to stdout or save to file.
    if lines:
        text = "\n".join(map(str, lines))
        if args.output:
            with open(args.output, "w") as fp:
                fp.write(text)
        else:
            print(text)


if __name__ == "__main__":
    __main__()
