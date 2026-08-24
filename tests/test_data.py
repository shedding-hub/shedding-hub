import jsonschema
import numpy as np
import os
from pathlib import Path
import pytest
import re
import requests
import time
import yaml

DATA_PATHS = list(Path("data").glob("*/*.yaml"))
VALID_EXAMPLE_PATHS = list(Path("tests/examples").glob("valid_*.yaml"))
INVALID_EXAMPLE_PATHS = list(Path("tests/examples").glob("invalid_*.yaml"))

# Every dataset's doi is resolved against the publisher, so one run makes as
# many requests as there are datasets -- 93 and climbing. Publishers drop
# connections under that, and a single dropped connection failed a whole
# 9-minute job with "RemoteDisconnected('Remote end closed connection without
# response')". Retried rather than re-run: a transient refusal says nothing
# about whether the doi resolves.
DOI_RETRIES = 3
DOI_BACKOFF_SECONDS = 2
REQUEST_TIMEOUT_SECONDS = 30


def _get_with_retry(url: str, **kwargs) -> requests.Response:
    """GET a URL, retrying only transport failures -- never a real HTTP status."""
    for attempt in range(DOI_RETRIES):
        try:
            return requests.get(url, timeout=REQUEST_TIMEOUT_SECONDS, **kwargs)
        except requests.exceptions.RequestException:
            # The last attempt raises: a URL that never answers is a genuine
            # failure of this check, not something to swallow.
            if attempt == DOI_RETRIES - 1:
                raise
            time.sleep(DOI_BACKOFF_SECONDS * (attempt + 1))


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

    # If there is a doi, validate it. Resolving it is a network round trip per
    # dataset, and three CI jobs load this file, so a pull request used to make
    # that trip three times over. Setting SHEDDING_HUB_SKIP_LINK_CHECKS drops
    # only the request -- every offline check below still runs, including the
    # requirement that a doi or url be present at all. The one job that leaves
    # it unset (data-validation) keeps the guarantee for the whole pull request.
    skip_link_checks = os.environ.get("SHEDDING_HUB_SKIP_LINK_CHECKS") == "1"
    doi_or_url = False
    doi = data.get("doi")
    if doi:
        if not skip_link_checks:
            response = _get_with_retry(f"https://doi.org/{doi}", allow_redirects=False)
            assert response.status_code == 302, f"doi `{doi}` could not be resolved."
        doi_or_url = True
    url = data.get("url")
    if url:
        if not skip_link_checks:
            response = _get_with_retry(url)
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
