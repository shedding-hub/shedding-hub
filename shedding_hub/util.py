import re
import textwrap


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
