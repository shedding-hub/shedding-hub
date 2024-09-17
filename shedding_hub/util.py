import pathlib
import re
import requests
import textwrap
from typing import Optional
import yaml


def normalize_str(
    value: str, *, dedent: bool = True, strip: bool = True, unwrap: bool = True
) -> str:
    """
    Normalize a string.

    Args:
        value: String to normalize.
        dedent: Remove any common leading whitespace.
        strip: Remove leading and trailing whitespace.
        unwrap: Unwrap lines separated by a single line break.
    """
    if dedent:
        value = textwrap.dedent(value)
    if unwrap:
        value = re.sub(r"(?<!\n) *\n(?!\n)", " ", value)
    if strip:
        value = value.strip()
    return value


def load_dataset(
    dataset: str,
    *,
    repo: str = "shedding-hub/shedding-hub",
    ref: Optional[str] = None,
    pr: Optional[int] = None,
    local: Optional[str] = None,
) -> dict:
    """
    Load a dataset from GitHub or a local directory.

    Args:
        dataset: Dataset identifier, e.g., :code:`woelfel2020virological`.
        repo: GitHub repository to load data from.
        ref: Git reference to load. Defaults to the most recent data on the :code:`main`
            branch of https://github.com/shedding-hub/shedding-hub and is automatically
            fetched if a :code:`pr` number is specified.
        pr: Pull request to fetch data from.
        local: Local directory to load data from.

    Returns:
        Loaded dataset.
    """
    # Check that at most one of `ref`, `pr`, and `local` is given.
    specified = {"ref": ref, "pr": pr, "local": local}
    n_specified = sum(1 if x else 0 for x in specified.values())
    if n_specified > 1:
        raise ValueError(
            f"At most one of `ref`, `pr`, or `local` may be specified; got {specified}."
        )

    # If we have a local file, just read it.
    if local:
        path = (pathlib.Path(local) / dataset / dataset).with_suffix(".yaml")
        with path.open() as fp:
            return yaml.safe_load(fp)

    # If a PR is specified, resolve it so we can get the relevant file.
    if pr:
        response = requests.get(f"https://api.github.com/repos/{repo}/pulls/{pr}")
        response.raise_for_status()
        response = response.json()
        repo = response["head"]["repo"]["full_name"]
        # Get the sha rather than just the ref because the branch may have been deleted,
        # but the commit will exist.
        ref = response["head"]["sha"]

    # Download the contents and parse the file.
    ref = ref or "main"
    response = requests.get(
        f"https://raw.githubusercontent.com/{repo}/{ref}/data/{dataset}/{dataset}.yaml"
    )
    # Backwards compatibility before change of folder structure.
    if response.status_code == 404:
        response = requests.get(
            f"https://raw.githubusercontent.com/{repo}/{ref}/data/{dataset}.yaml"
        )
    response.raise_for_status()
    return yaml.safe_load(response.text)
