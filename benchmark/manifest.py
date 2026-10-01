"""
Build `benchmark/manifest.csv`: one row per study in the benchmark analysis set.

The manifest is read from a frozen release (a git tag), not from the working
tree, so a later edit to a dataset cannot change what the benchmark is scored
against. Every column except `pmid`, `has_figure_csv` and `evidence_tier`
comes from the dataset's own YAML at that tag.

`pmid` is looked up from the DOI through NCBI E-utilities, because the
extraction agent is driven by PMID. A hit is kept only if PubMed's record
carries the same DOI. Results are cached in `benchmark/pmid_lookup.csv`, so a
rebuild needs no network; `--refresh-pmids` looks up the uncached ones.

`has_figure_csv` records whether a curator-prepared input CSV exists for the
study. The CSVs live in the private agent repository, so only the list of
study ids is kept here, in `benchmark/input_csv_studies.txt`.
`--agents-repo` rewrites that list from a tag of that repository.

`input_csv_origin` says what kind of data the input CSV carries. It follows
from `data_source` with a fixed precedence, because the same CSV holds
whichever of these the curator laid out: author-shared data, then data from a
repository, then values digitized from figures. A study with more than one of
these is flagged in `notes`. A study that has a CSV but none of those sources
is `table_supplement`: the curator transcribed the paper's tables or
supplement into the CSV. `none` means the study has no input CSV.

`evidence_tier` is filled when the evidence bundles are built.
"""

import argparse
import csv
import json
import pathlib
import re
import subprocess
import sys
import time
import unicodedata
import urllib.parse
import urllib.request

import yaml

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
BENCHMARK = REPO_ROOT / "benchmark"

# A population-scale routine-testing dataset, outside the analysis set.
EXCLUDED = {"jones2021estimating"}

COLUMNS = [
    "study_id",
    "doi",
    "pmid",
    "era",
    "pipeline_version",
    "extraction_model",
    "data_source",
    "has_figure_csv",
    "input_csv_origin",
    "evidence_tier",
    "in_analysis_A",
    "notes",
]

# Highest precedence first.
CSV_ORIGINS = [
    ("author_shared", "author_shared"),
    ("data_repository", "data_repository"),
    ("figure", "figure_digitized"),
]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PUBMED_URL = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")


def git(*args: str, cwd: pathlib.Path = REPO_ROOT) -> str:
    out = subprocess.run(["git", *args], capture_output=True, cwd=cwd, check=True)
    return out.stdout.decode("utf-8")


def release_datasets(ref: str) -> dict:
    """study_id -> parsed dataset, for every dataset at the release."""
    datasets = {}
    for line in git("ls-tree", "-r", "--name-only", ref, "--", "data").splitlines():
        parts = line.split("/")
        if len(parts) == 3 and parts[2] == parts[1] + ".yaml":
            datasets[parts[1]] = yaml.safe_load(git("show", f"{ref}:{line}"))
    return datasets


def ascii_id(name: str) -> str:
    """Study id without accents: the agent repository keeps a few with them."""
    decomposed = unicodedata.normalize("NFKD", name)
    return "".join(ch for ch in decomposed if not unicodedata.combining(ch))


def input_csv_studies(agents_repo: pathlib.Path, agents_ref: str) -> list[str]:
    """Study ids that have a `figure_data/<study>.csv` at the agent repo tag."""
    listing = git(
        "-c",
        "core.quotepath=false",
        "ls-tree",
        "--name-only",
        agents_ref,
        "figure_data/",
        cwd=agents_repo,
    )
    names = [pathlib.PurePosixPath(line).name for line in listing.splitlines()]
    return sorted(ascii_id(n[: -len(".csv")]) for n in names if n.endswith(".csv"))


def fetch_json(endpoint: str, **params: str) -> dict:
    query = urllib.parse.urlencode({**params, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}{endpoint}?{query}", timeout=30) as resp:
        payload = json.load(resp)
    time.sleep(0.4)  # NCBI allows three requests a second without a key
    return payload


def lookup_pmid(doi: str) -> str:
    """The PMID whose PubMed record carries this DOI, or '' if there is none."""
    found = fetch_json("esearch.fcgi", db="pubmed", term=f"{doi}[AID]")
    for pmid in found.get("esearchresult", {}).get("idlist", []):
        record = fetch_json("esummary.fcgi", db="pubmed", id=pmid)["result"][pmid]
        dois = [
            a["value"].lower()
            for a in record.get("articleids", [])
            if a.get("idtype") == "doi"
        ]
        if doi.lower() in dois:
            return pmid
    return ""


def read_pmid_cache(path: pathlib.Path) -> dict:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8", newline="") as fh:
        return {row["study_id"]: row for row in csv.DictReader(fh)}


def write_pmid_cache(path: pathlib.Path, cache: dict) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["study_id", "doi", "pmid", "source"], lineterminator="\n"
        )
        writer.writeheader()
        for study in sorted(cache):
            writer.writerow(cache[study])


def resolve_pmids(datasets: dict, cache: dict, refresh: bool) -> None:
    """Fill the cache for every study it does not yet cover."""
    for study, dataset in datasets.items():
        doi = (dataset.get("doi") or "").strip()
        cached = cache.get(study)
        if cached and cached["doi"] == doi:
            continue
        from_url = PUBMED_URL.search(dataset.get("url") or "")
        if from_url:
            row = {"pmid": from_url.group(1), "source": "url"}
        elif not doi:
            row = {"pmid": "", "source": "no_doi"}
        elif refresh:
            pmid = lookup_pmid(doi)
            row = {"pmid": pmid, "source": "eutils" if pmid else "not_in_pubmed"}
        else:
            raise SystemExit(
                f"{study}: no cached PMID for {doi}; rerun with --refresh-pmids"
            )
        cache[study] = {"study_id": study, "doi": doi, **row}


def manifest_row(study: str, dataset: dict, pmid: str, has_csv: bool) -> dict:
    curation = dataset["curation"]
    sources = curation["data_source"]
    era = curation["method"]
    notes = []

    present = [origin for source, origin in CSV_ORIGINS if source in sources]
    if len(present) > 1:
        notes.append("input CSV may mix: " + " + ".join(present))
    if present:
        origin = present[0]
    else:
        origin = "table_supplement" if has_csv else "none"
    if present and not has_csv:
        notes.append("no input CSV yet")
    if not pmid:
        notes.append("no PMID")

    extraction = [m for m in curation.get("models", []) if m.endswith("(extraction)")]
    return {
        "study_id": study,
        "doi": dataset.get("doi") or "",
        "pmid": pmid,
        "era": era,
        "pipeline_version": curation.get("pipeline_version", ""),
        "extraction_model": extraction[0].split(" ")[0] if extraction else "",
        "data_source": ";".join(sources),
        "has_figure_csv": "yes" if has_csv else "no",
        "input_csv_origin": origin,
        "evidence_tier": "",
        "in_analysis_A": "yes" if era == "manual" else "no",
        "notes": "; ".join(notes),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="v1.1.0", help="Frozen release tag.")
    parser.add_argument(
        "--agents-repo",
        type=pathlib.Path,
        default=None,
        help="Agent repository; rewrites the list of studies with an input CSV.",
    )
    parser.add_argument(
        "--agents-ref",
        default="pipeline-v2-benchmark",
        help="Tag of the agent repository to list input CSVs from.",
    )
    parser.add_argument(
        "--refresh-pmids",
        action="store_true",
        help="Look up PMIDs that are not cached (needs network).",
    )
    args = parser.parse_args()

    csv_list = BENCHMARK / "input_csv_studies.txt"
    if args.agents_repo is not None:
        studies = input_csv_studies(args.agents_repo, args.agents_ref)
        csv_list.write_text("\n".join(studies) + "\n", encoding="utf-8", newline="\n")
    with_csv = set(csv_list.read_text(encoding="utf-8").split())

    datasets = {
        study: dataset
        for study, dataset in release_datasets(args.release).items()
        if study not in EXCLUDED
    }

    cache_path = BENCHMARK / "pmid_lookup.csv"
    cache = read_pmid_cache(cache_path)
    resolve_pmids(datasets, cache, args.refresh_pmids)
    cache = {study: cache[study] for study in datasets}
    write_pmid_cache(cache_path, cache)

    unknown = with_csv - set(datasets) - EXCLUDED
    if unknown:
        print(f"input CSVs with no dataset in {args.release}: {sorted(unknown)}")

    rows = [
        manifest_row(study, datasets[study], cache[study]["pmid"], study in with_csv)
        for study in sorted(datasets)
    ]
    with (BENCHMARK / "manifest.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"{len(rows)} studies written to benchmark/manifest.csv ({args.release})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
