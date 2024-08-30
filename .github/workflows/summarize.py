import pandas as pd
import subprocess
import yaml


def with_name(x: pd.Series, name: str) -> pd.Series:
    x.name = name
    return x


def __main__() -> None:
    # Get all changed data files compared with the main branch.
    output = subprocess.check_output(["git", "diff", "--name-only", "main"], text=True)
    filenames = [
        name
        for name in output.splitlines()
        if name.startswith("data/") and name.endswith(".yaml")
    ]

    # Iterate over files and summarize their content.
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
        lines = [
            f"🔄 Summary for changed file `{filename}`:",
            "",
            "```",
            pd.DataFrame(
                [
                    with_name(grouped.count().value, "n_samples"),
                    with_name(grouped.nunique().participant, "n_unique_participants"),
                ]
            ).T,
            "```",
        ]

        print("\n".join(map(str, lines)))


if __name__ == "__main__":
    __main__()
