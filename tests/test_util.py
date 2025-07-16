import hashlib
import io
import pytest
from shedding_hub import util
import yaml


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


@pytest.mark.parametrize(
    "kwargs, expected_sha1",
    [
        # This fetches the *current* dataset which may change and hence need adjustment
        # of the hash. However, this dataset is relatively stable and is unlikely to
        # need updates.
        (
            {"dataset_ID": "woelfel2020virological"},
            "d93f77da4babf4247f7775121adc749337f53f31",
        ),
        # An old version of the Woelfel dataset from a PR before folder restructuring.
        (
            {"dataset_ID": "woelfel2020", "pr": 1},
            "7bf4dabae5ba95b06a62f6102fee571b46068786",
        ),
        # The same old version of the Woelfel dataset using a commit reference.
        (
            {"dataset_ID": "woelfel2020", "ref": "534c30a"},
            "7bf4dabae5ba95b06a62f6102fee571b46068786",
        ),
        # Invalid because requesting local and pr.
        (
            {"dataset_ID": "woelfel2020virological", "local": "data", "pr": 7},
            ValueError,
        ),
        # Load from local directory.
        (
            {"dataset_ID": "woelfel2020virological", "local": "data"},
            "d93f77da4babf4247f7775121adc749337f53f31",
        ),
    ],
)
def test_load(kwargs: dict, expected_sha1: str) -> None:
    if isinstance(expected_sha1, str):
        data = util.load_dataset(**kwargs)
        assert "title" in data
        stream = io.StringIO()
        yaml.safe_dump(data, stream)
        actual_sha1 = hashlib.sha1(stream.getvalue().encode()).hexdigest()
        assert actual_sha1 == expected_sha1
    else:
        with pytest.raises(expected_sha1):
            data = util.load_dataset(**kwargs)


def test_str_representer() -> None:
    x = {"a": util.folded_str("foo\nbar\n"), "b": util.literal_str("foo\nbar\n")}
    dumped = yaml.dump(x)
    assert (
        dumped.strip()
        == """
a: >
  foo

  bar
b: |
  foo
  bar
""".strip()
    )
    y = yaml.safe_load(io.StringIO(dumped))
    assert x == y
