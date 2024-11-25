import argparse
import numpy as np
from typing import Mapping, Sequence
import yaml


def assert_close(x1, x2, path=()):
    template = f"values differ at path `/{'/'.join(map(str, path))}`:"
    assert (
        x1.__class__ == x2.__class__
    ), f"{template} class {x1.__class__} does not match {x2.__class__}"
    if isinstance(x1, Mapping):
        missing = set(x1) - set(x2)
        assert not missing, f"{template} keys {missing} in x1 are not in x2"
        extra = set(x2) - set(x1)
        assert not extra, f"{template} keys {extra} in x2 are not in x1"
        for key, value in x1.items():
            assert_close(value, x2[key], path + (key,))
    elif isinstance(x1, Sequence) and not isinstance(x1, str):
        assert len(x1) == len(
            x2
        ), f"{template} x1 has length {len(x1)} but x2 has length {len(x2)}"
        for i, (value1, value2) in enumerate(zip(x1, x2)):
            assert_close(value1, value2, path + (i,))
    elif isinstance(x1, float):
        assert np.allclose(x1, x2), f"{template} x1={x1} and x2={x2} are not close"
    else:
        assert x1 == x2, f"{template} x1={x1} and x2={x2} are different"


def __main__(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("filename1")
    parser.add_argument("filename2")
    args = parser.parse_args(argv)

    with open(args.filename1) as fp:
        x1 = yaml.safe_load(fp)
    with open(args.filename2) as fp:
        x2 = yaml.safe_load(fp)

    assert_close(x1, x2)


if __name__ == "__main__":
    __main__()
