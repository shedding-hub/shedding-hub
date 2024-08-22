import pytest
from shedding_hub import util


@pytest.mark.parametrize(
    "value, kwargs, expected",
    [
        ("    asdf\n    jkl;", {}, "asdf jkl;"),
        ("    asdf\n    jkl;", {"dedent": False}, "asdf     jkl;"),
        ("    asdf\n    jkl;", {"dedent": False, "strip": False}, "    asdf     jkl;"),
        ("    asdf\n    jkl;", {"dedent": False, "unwrap": False}, "asdf\n    jkl;"),
        (
            "    asdf\n    jkl;",
            {"dedent": False, "unwrap": False, "strip": False},
            "    asdf\n    jkl;",
        ),
    ],
)
def test_normalize_str(value: str, kwargs: dict, expected: str) -> None:
    assert util.normalize_str(value, **kwargs) == expected

