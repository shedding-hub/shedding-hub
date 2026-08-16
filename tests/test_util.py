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
            {"dataset": "woelfel2020virological"},
            "7a7453c9259f1043657f8d19fbfdf2f69aaf5a30",
        ),
        # An old version of the Woelfel dataset from a PR before folder restructuring.
        (
            {"dataset": "woelfel2020", "pr": 1},
            "dbf4335ebae87445c821a0772178180a596f5615",
        ),
        # The same old version of the Woelfel dataset using a commit reference.
        (
            {"dataset": "woelfel2020", "ref": "534c30a"},
            "dbf4335ebae87445c821a0772178180a596f5615",
        ),
        # Invalid because requesting local and pr.
        (
            {"dataset": "woelfel2020virological", "local": "data", "pr": 7},
            ValueError,
        ),
        # Load from local directory.
        (
            {"dataset": "woelfel2020virological", "local": "data"},
            "7a7453c9259f1043657f8d19fbfdf2f69aaf5a30",
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
    assert dumped.strip() == """
a: >
  foo

  bar
b: |
  foo
  bar
""".strip()
    y = yaml.safe_load(io.StringIO(dumped))
    assert x == y


def test_github_api_headers_uses_a_token_when_offered(monkeypatch) -> None:
    monkeypatch.delenv("GH_TOKEN", raising=False)
    monkeypatch.setenv("GITHUB_TOKEN", "sekrit")
    assert util._github_api_headers() == {"Authorization": "Bearer sekrit"}

    # GH_TOKEN is what the gh CLI exports, so it is honoured too.
    monkeypatch.delenv("GITHUB_TOKEN")
    monkeypatch.setenv("GH_TOKEN", "other")
    assert util._github_api_headers() == {"Authorization": "Bearer other"}


def test_github_api_headers_are_empty_without_a_token(monkeypatch) -> None:
    """Unauthenticated still has to work: the package is used outside CI."""
    monkeypatch.delenv("GITHUB_TOKEN", raising=False)
    monkeypatch.delenv("GH_TOKEN", raising=False)
    assert util._github_api_headers() == {}


def test_token_goes_only_to_the_api_host(monkeypatch) -> None:
    """A credential must not be sent where it is not needed.

    raw.githubusercontent.com serves public content without one, so the token
    belongs on the api.github.com request resolving the pull request and
    nowhere else.
    """
    monkeypatch.setenv("GITHUB_TOKEN", "sekrit")
    seen = []

    class _Response:
        status_code = 200
        text = "title: x\n"

        def raise_for_status(self) -> None:
            pass

        def json(self) -> dict:
            return {"head": {"repo": {"full_name": "o/r"}, "sha": "abc123"}}

    def _fake_get(url, **kwargs):
        seen.append((url, kwargs.get("headers") or {}, kwargs.get("timeout")))
        return _Response()

    monkeypatch.setattr(util.requests, "get", _fake_get)
    util.load_dataset("somestudy", pr=1)

    api = [s for s in seen if "api.github.com" in s[0]]
    raw = [s for s in seen if "raw.githubusercontent.com" in s[0]]
    assert api and raw, f"expected both hosts, saw {[s[0] for s in seen]}"
    assert api[0][1] == {"Authorization": "Bearer sekrit"}
    assert all(headers == {} for _, headers, _ in raw), "token leaked to raw host"
    assert all(timeout is not None for _, _, timeout in seen), "a call had no timeout"
