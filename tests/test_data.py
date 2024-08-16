import jsonschema
from pathlib import Path
import pytest
import yaml


DATA_PATHS = [
    path for path in Path("data").glob("*.yaml") if path.name != ".schema.yaml"
]
VALID_EXAMPLE_PATHS = list(Path("tests/examples").glob("valid_*.yaml"))
INVALID_EXAMPLE_PATHS = list(Path("tests/examples").glob("invalid_*.yaml"))


def load_and_validate(path: Path):
    """
    Load and validate a dataset.
    """
    with open("data/.schema.yaml") as fp:
        schema = yaml.safe_load(fp)
    with path.open() as fp:
        data = yaml.safe_load(fp)
    jsonschema.validate(data, schema)

    # Ensure there is exactly one of `analyte` and `analytes`.
    has_analyte = "analyte" in data
    has_analytes = "analytes" in data
    if has_analyte == has_analytes:
        raise ValueError("Data must have exactly one of `analyte` or `analytes` field.")

    for i, participant in enumerate(data["participants"]):
        for j, measurement in enumerate(participant["measurements"]):
            if has_analyte and "analyte" in measurement:
                raise ValueError(
                    "Data declared only a single analyte using the top-level `analyte` "
                    "field, and individual measurements must not have an `analyte` "
                    f"field. Measurement {j} for patient {i} has an `analyte` field."
                )
            elif has_analytes and "analyte" not in measurement:
                raise ValueError(
                    "Data declared multiple analytes using the top-level `analytes` "
                    "field, and each individual measurement must have an `analyte` "
                    f"field. Measurement {j} for patient {i} does not has an `analyte` "
                    "field."
                )
            elif has_analytes and measurement["analyte"] not in data["analytes"]:
                raise ValueError(
                    f"Data declared valid analytes {set(data['analytes'])}. "
                    f"Measurement {j} for patient {i} declares the invalid analyte "
                    f"`{measurement['analyte']}`."
                )

    if has_analytes:
        used_analytes = {
            measurement["analyte"]
            for participant in data["participants"]
            for measurement in participant["measurements"]
        }
        unused_analytes = set(data["analytes"]) - used_analytes
        if unused_analytes:
            raise ValueError(f"Data declared unused analytes {unused_analytes}.")


@pytest.mark.parametrize("path", DATA_PATHS, ids=[path.stem for path in DATA_PATHS])
def test_data_validity(path: Path) -> None:
    load_and_validate(path)


@pytest.mark.parametrize(
    "path", VALID_EXAMPLE_PATHS, ids=[path.stem for path in VALID_EXAMPLE_PATHS]
)
def test_valid_examples(path: Path) -> None:
    load_and_validate(path)


@pytest.mark.parametrize(
    "path", INVALID_EXAMPLE_PATHS, ids=[path.stem for path in INVALID_EXAMPLE_PATHS]
)
def test_invalid_examples(path: Path) -> None:
    with pytest.raises((ValueError, jsonschema.ValidationError)):
        load_and_validate(path)
