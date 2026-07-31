# A generated documentation site for the package — design

**Date:** 2026-07-31
**Status:** Draft (awaiting review)

## Problem

`shedding-hub` 0.2.0 exposes 42 public names. The documentation a user can find
describes 19 of them, and three of the things it describes do not work.

The package page on the project website, `shedding-hub.github.io/package.html`,
is hand-written HTML with a manually maintained "Function Reference". Measured
against the shipped package:

- **23 of 42 public names are absent**, including the entire fitting, catalog,
  simulation and selection surface added in 0.2.0.
- **Three documented calls raise `AttributeError`.** The page tells readers to
  call `sh.plot_clearance_curve(...)`, `sh.plot_detection_probability(...)` and
  `sh.plot_value_distribution_by_time(...)`. All three are implemented in
  `shedding_hub/viz.py` and none is exported from `__init__.py`. They were
  already in that state at the 0.1.3 release, so this is long-standing.

Nothing checks the page against the code, so it drifts silently and only a
reader discovers it. That is the actual defect; the missing entries are a
symptom.

Two supporting facts shape the fix. **Every one of the 42 public names has a
docstring** — none missing, and they are substantial: `simulate_shedding` runs
to 3,165 characters covering the dispersion rationale and the time-origin
caveats. The reference is therefore a rendering problem, not a writing one.
And **the existing long-form documentation is already markdown** —
`docs/modeling-methods.md` (244 lines) and `examples/simulating-shedding.md`
(422 lines, a jupytext notebook verified to execute) — so a markdown toolchain
reuses it without conversion.

## Decisions

**Read the Docs, with the website page reduced to a pointer.** This is the
arrangement pystan uses. It supplies hosting, search, and a version selector
without anyone maintaining them. `package.html` keeps install instructions and a
short quickstart and links out; its hand-maintained function list is deleted
rather than migrated, that list being the thing that rotted.

**MkDocs Material with mkdocstrings.** Both it and Sphinx run on Read the Docs;
the deciding factor is that the existing corpus is markdown, so it drops in
unchanged and stays editable by anyone who can write markdown. mkdocstrings
reads the Google-style docstrings (`Args:` / `Returns:` / `Raises:`) the package
already uses. Sphinx would need MyST to ingest the same files and carries a
larger configuration surface; its main advantage, intersphinx linking into
numpy and pandas, is not worth that here.

**The three orphans are exported, with tests written first.** Someone wrote
1,648 to 2,420 characters of docstring for each and put them on the website, so
they were intended to be public; the export was forgotten. They currently have
**zero test references anywhere in `tests/`**, against 63 for `plot_time_course`
and 16 for `plot_simulated_shedding`. Publishing them in a generated reference
makes them discoverable, so they get smoke tests before that happens rather
than after. `__all__` goes from 42 to 45.

**Latest-only versioning to begin with.** A per-version selector needs a git tag
per release, and the release process creates none — the existing `v1.0.x` tags
are repository releases on an independent line, having shipped pyproject
versions 0.1.0 and 0.1.2 respectively. Read the Docs will build `main`, and the
site carries a banner stating that it documents the development version and
naming the current release. Adding tags later turns the banner into a real
`stable`/`latest` split without rework.

## Structure

MkDocs defaults to `docs/` as its source, and that directory already holds
`docs/superpowers/specs/` and `docs/superpowers/plans/` — internal design
records, including candid accounts of bugs found and decisions reversed. Built
as-is they would all become public pages. The configuration excludes them
explicitly via `exclude_docs`, native since MkDocs 1.6.

```
mkdocs.yml
.readthedocs.yaml
docs/
  index.md              what the package is; install
  getting-started.md    load a dataset, summarise, visualise
  tutorial.md           the simulation walkthrough
  modeling-methods.md   already present, rendered unchanged
  reference/            one mkdocstrings stub per module
  superpowers/          excluded from the build
```

`docs/shedding_parameters.{json,csv}` are not markdown, so they are not built as
pages; they are served as static downloads, which is a small bonus rather than
a design goal.

## Single source, no copies

`docs/modeling-methods.md` is rendered where it already lives. The tutorial's
source of truth stays `examples/simulating-shedding.md` — the file the README
links to, and the one that can be executed to verify it still runs — and the
docs page includes it rather than duplicating it. A copy would drift, which is
the failure this whole design exists to correct. The jupytext front matter is
stripped at include time so it does not render.

## What stops it rotting again

A generated reference cannot disagree with the code about a signature. It can
still silently omit a module nobody added to the nav, which is the same rot in a
new place. So a test asserts that **every name in `__all__` appears in the
reference**, failing with the missing names when it does not. That test, not the
tooling, is what prevents `package.html`'s history repeating.

The docs build itself runs in CI on pull requests, so a broken include, an
unresolvable mkdocstrings target, or a nav entry pointing at a deleted file
fails the build rather than the published site.

## Testing

- Every name in `__all__` is reachable in the generated reference; the failure
  message names the omissions.
- The three newly exported functions each return a `Figure` and honour their
  `biomarker` / `specimen` filters — the smoke coverage they have never had.
- `mkdocs build --strict` succeeds, so warnings (broken internal links, missing
  nav targets, unresolved references) are failures.
- The tutorial page renders the notebook's content rather than a stub, checked
  by asserting a known phrase from it appears in the built HTML.
- The existing README doctests and the full suite stay green; exporting three
  names must not alter any existing behaviour.

## Out of scope

Rewriting `docs/modeling-methods.md`, which is current and stands as written.
Migrating the website's other pages — datasets, curation, vocabulary, models,
team — which are hand-written by intent and do not track code. Per-version
documentation and the release tagging it needs. Intersphinx linking. Publishing
the internal design records under `docs/superpowers/`.
