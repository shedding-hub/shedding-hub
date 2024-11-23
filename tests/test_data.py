import jsonschema
import numpy as np
from pathlib import Path
import pytest
import re
import requests
import yaml


DATA_PATHS = list(Path("data").glob("*/*.yaml"))
VALID_EXAMPLE_PATHS = list(Path("tests/examples").glob("valid_*.yaml"))
INVALID_EXAMPLE_PATHS = list(Path("tests/examples").glob("invalid_*.yaml"))


def load_and_validate(path: Path, skip_filename_check: bool = False):
    """
    Load and validate a dataset.
    """
    if not skip_filename_check:
        assert (
            path.stem == path.parent.stem
        ), "The data filename must match the parent folder."
        assert (
            str(path) == str(path).lower()
        ), "Data paths and filenames should be lowercase."
        assert " " not in str(path), "Data paths should not contain spaces."
        assert re.match(
            r"[a-z]+\d{4}[a-z]+\.yaml", path.name
        ), "File name must match the pattern `[author][year][first word of title]`."

    with open("data/.schema.yaml") as fp:
        schema = yaml.safe_load(fp)
    with path.open() as fp:
        data = yaml.safe_load(fp)
    jsonschema.validate(data, schema)

    # If there is a doi, validate it.
    doi_or_url = False
    doi = data.get("doi")
    if doi:
        response = requests.get(f"https://doi.org/{doi}", allow_redirects=False)
        assert response.status_code == 302, f"doi `{doi}` could not be resolved."
        doi_or_url = True
    url = data.get("url")
    if url:
        response = requests.get(url)
        response.raise_for_status()
        doi_or_url = True
    assert doi_or_url, "At least one of `doi` or `url` must be given."

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

            value = measurement["value"]
            if not isinstance(value, str) and np.isnan(value):
                raise ValueError(f"Measurement {j} for patient {i} has nan `value`.")
            time = measurement.get("time", 0)
            if not isinstance(time, str) and np.isnan(time):
                raise ValueError(f"Measurement {j} for patient {i} has nan `time`.")

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
    load_and_validate(path, skip_filename_check=True)


@pytest.mark.parametrize(
    "path", INVALID_EXAMPLE_PATHS, ids=[path.stem for path in INVALID_EXAMPLE_PATHS]
)
def test_invalid_examples(path: Path) -> None:
    with pytest.raises((AssertionError, ValueError, jsonschema.ValidationError)):
        load_and_validate(path, skip_filename_check=True)
