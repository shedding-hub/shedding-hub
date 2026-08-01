"""
The reference must not drift from the package.

A generated reference cannot disagree with the code about a signature, but it
can silently omit a module nobody added to the nav -- the same rot that left
the project website documenting 19 of 42 names, three of which did not exist.
"""

import pathlib
import re

import shedding_hub as sh

REFERENCE = pathlib.Path(__file__).parent.parent / "docs" / "reference"


def _documented_names() -> set:
    text = "\n".join(p.read_text(encoding="utf-8") for p in REFERENCE.glob("*.md"))
    # Anchored to a `::: shedding_hub.<name>` directive on its own line, rather
    # than an unanchored substring match: otherwise a name that is a strict
    # prefix of another (e.g. `calc_shedding_duration` inside
    # `calc_shedding_durations`) could be deleted from the reference and this
    # test would still pass.
    entries = set(re.findall(r"^::: shedding_hub\.(\w+)\s*$", text, re.M))
    return {name for name in sh.__all__ if name in entries}


def test_every_public_name_appears_in_the_reference():
    missing = sorted(set(sh.__all__) - _documented_names())
    assert not missing, (
        f"{len(missing)} public name(s) are exported but absent from "
        f"docs/reference/: {missing}. Add them to the matching module page."
    )


def test_every_callable_has_a_runnable_example():
    """
    Constants are exempt, and cannot not be.

    ``inspect.getdoc`` on a module-level constant returns its *container type's*
    docstring -- ``tuple.__doc__`` for MODELS, ``dict.__doc__`` for PARAM_NAMES
    and REFERENCE_EVENT_CLASSES -- so no example can ever be attached to the
    name itself. Those three carry executed examples in their defining module's
    docstring instead, and test_every_public_name_appears_in_the_reference still
    holds all 45 to appearing in the reference.
    """
    import inspect

    needs_example = [
        n
        for n in sh.__all__
        if callable(getattr(sh, n)) or inspect.isclass(getattr(sh, n))
    ]
    missing = sorted(
        n for n in needs_example if ">>>" not in (inspect.getdoc(getattr(sh, n)) or "")
    )
    assert not missing, (
        f"{len(missing)} callable(s) have no worked example: {missing}. "
        "Add an Example: block; it is executed by --doctest-modules."
    )
