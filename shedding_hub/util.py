import difflib
import pathlib
import re
import requests
import textwrap
from typing import Optional
import warnings
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
            data = yaml.safe_load(fp)
        data["dataset_id"] = dataset
        return data

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
    data = yaml.safe_load(response.text)
    data["dataset_id"] = dataset
    return data


def check_dataset(
    *,
    doi: Optional[str] = None,
    title: Optional[str] = None,
    local: Optional[str] = None,
    similarity_threshold: float = 0.6,
) -> bool:
    """
    Check whether a paper is in the curated datasets.

    Args:
        doi: DOI of the paper to check.
        title: Title of the paper to check.
        local: Local directory containing datasets. Defaults to the ``data``
            directory in the repository root.
        similarity_threshold: Minimum similarity ratio (0 to 1) for reporting
            near-matches when no exact title match is found.

    Returns:
        True if the paper is found (exact DOI or exact title match), False
        otherwise.  When a title is provided and no exact match is found, a
        warning is issued for the most similar dataset above the threshold.
    """
    if doi is None and title is None:
        raise ValueError("At least one of `doi` or `title` must be specified.")

    # Resolve data directory.
    if local:
        data_dir = pathlib.Path(local)
    else:
        data_dir = pathlib.Path(__file__).parent.parent / "data"

    if not data_dir.is_dir():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    # Collect metadata from all datasets.  Only the first few lines of each
    # file are read to avoid loading very large YAML files entirely.
    _header_keys = {"title", "doi", "url"}
    _stop_keys = ("analyte:", "analytes:", "participants:")
    datasets = []
    for yaml_path in sorted(data_dir.glob("*/*.yaml")):
        if yaml_path.name.startswith("."):
            continue
        metadata: dict = {"title": "", "doi": "", "url": ""}
        current_key = None
        with yaml_path.open(encoding="utf-8") as fp:
            for line in fp:
                # Stop once we reach a section that comes after the header.
                if line.startswith(_stop_keys):
                    break
                # Check if this line starts a new top-level key.
                matched_key = False
                for key in _header_keys:
                    if line.startswith(f"{key}:"):
                        metadata[key] = line.split(":", 1)[1].strip()
                        current_key = key
                        matched_key = True
                        break
                if matched_key:
                    continue
                # Continuation line for the current key (indented).
                if current_key and line.startswith((" ", "\t")):
                    metadata[current_key] += " " + line.strip()
                else:
                    current_key = None
        datasets.append(metadata)

    # Check for exact DOI match.
    if doi is not None:
        doi_normalized = doi.strip().lower()
        for ds in datasets:
            if ds["doi"] and ds["doi"].strip().lower() == doi_normalized:
                return True

    # Check for exact title match.
    if title is not None:
        title_normalized = title.strip().lower()
        for ds in datasets:
            if ds["title"].strip().lower() == title_normalized:
                return True

        # No exact match — look for the closest title above the threshold.
        best_match = None
        best_ratio = 0.0
        for ds in datasets:
            ratio = difflib.SequenceMatcher(
                None, title_normalized, ds["title"].strip().lower()
            ).ratio()
            if ratio > best_ratio:
                best_ratio = ratio
                best_match = ds

        if best_match and best_ratio >= similarity_threshold:
            identifier = best_match["doi"] or best_match["url"]
            warnings.warn(
                f"No exact title match found, but a similar dataset exists: "
                f'"{best_match["title"]}" ({identifier}), '
                f"similarity: {best_ratio:.2f}."
            )

    return False


class folded_str(str):
    """
    Folded string in yaml representation.
    """


class literal_str(str):
    """
    Literal string in yaml representation.
    """


def _folded_str_representer(dumper, data):
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style=">")


def _literal_str_representer(dumper, data):
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")


yaml.add_representer(folded_str, _folded_str_representer)
yaml.add_representer(literal_str, _literal_str_representer)
