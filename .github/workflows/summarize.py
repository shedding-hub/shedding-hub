import argparse
import pandas as pd
import subprocess
import yaml


def with_name(x: pd.Series, name: str) -> pd.Series:
    x.name = name
    return x


def __main__(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", nargs="?")
    args = parser.parse_args(argv)

    # Get all changed data files compared with the main branch.
    output = subprocess.check_output(
        ["git", "diff", "--name-only", "origin/main"], text=True
    )
    filenames = [
        name
        for name in output.splitlines()
        if name.startswith("data/") and name.endswith(".yaml")
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
                    ]
                ).T,
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
