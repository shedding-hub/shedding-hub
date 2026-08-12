"""
Pytest configuration for the whole repository.

The Agg backend must be selected before any test or doctest imports pyplot.
`tests/conftest.py` only governs the `tests/` tree, so module doctests under
`shedding_hub/` were relying on `tests/` being collected first and setting the
backend as a process-wide side effect. Nothing guaranteed that ordering, and
running the module doctests alone failed intermittently with a Tk error.
"""

import matplotlib

matplotlib.use("Agg")
