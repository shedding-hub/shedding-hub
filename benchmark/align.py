"""Align an agent draft (S0) with a reference dataset (S2).

Three steps, each feeding the next:

1. analytes are matched on what they describe (biomarker, specimen, unit,
   reference event), never on their key, which is a free-form name;
2. participants are matched by the measurements they share, since their order
   is arbitrary;
3. within a matched participant and analyte, measurements are matched on time
   and then on value, each used at most once.

Nothing here decides which side is right. `discrepancies.py` turns the
alignment into rows of differences.

The tolerances are the proposed defaults of the analysis plan. They are not
final until the plan is tagged.
"""

from __future__ import annotations

import difflib
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linear_sum_assignment

EPS = 1e-9
CT_UNIT = "cycle threshold"
POSITIVE = {"positive", "weak positive", "strong positive"}
ANALYTE_FIELDS = (
    "biomarker",
    "specimen",
    "unit",
    "reference_event",
    "gene_target",
    "limit_of_detection",
    "limit_of_quantification",
)


@dataclass(frozen=True)
class Tolerances:
    """How far apart two values may be and still count as the same."""

    time: float  # days
    log10: float | None  # concentrations; None means use `relative`
    relative: float  # concentrations, as a share of the larger value
    ct: float  # cycles


FIGURE = Tolerances(time=0.5, log10=0.1, relative=0.005, ct=0.5)
REPORTED = Tolerances(time=0.0, log10=None, relative=0.005, ct=0.0)


def tolerances_for(input_csv_origin: str) -> Tolerances:
    """Digitized values get the wider tolerances; reported values do not."""
    return FIGURE if input_csv_origin == "figure_digitized" else REPORTED


def is_number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def value_kind(value) -> str:
    """numeric, negative, positive, inconclusive, or other."""
    if is_number(value):
        return "numeric"
    if value == "negative":
        return "negative"
    if value in POSITIVE:
        return "positive"
    if value == "inconclusive":
        return "inconclusive"
    return "other"


def values_agree(a, b, unit: str | None, tol: Tolerances) -> bool:
    if not (is_number(a) and is_number(b)):
        return type(a) is type(b) and a == b
    if unit == CT_UNIT:
        return abs(a - b) <= tol.ct + EPS
    if tol.log10 is not None and a > 0 and b > 0:
        return abs(math.log10(a) - math.log10(b)) <= tol.log10 + EPS
    return abs(a - b) <= tol.relative * max(abs(a), abs(b)) + EPS


def value_distance(a, b) -> float:
    """How far apart two values are, for choosing the nearest candidate."""
    if is_number(a) and is_number(b):
        if a > 0 and b > 0:
            return abs(math.log10(a) - math.log10(b))
        return abs(a - b)
    return 0.0 if a == b else math.inf


def times_agree(a, b, tol: Tolerances) -> bool:
    if is_number(a) and is_number(b):
        return abs(a - b) <= tol.time + EPS
    return not is_number(a) and not is_number(b) and a == b


def time_distance(a, b) -> float:
    if is_number(a) and is_number(b):
        return abs(a - b)
    return 0.0 if a == b else math.inf


# --------------------------------------------------------------------- analytes


def analytes_of(dataset) -> dict:
    """The dataset's analytes as {key: specification}, tolerating bad drafts."""
    analytes = dataset.get("analytes") if isinstance(dataset, dict) else None
    if not isinstance(analytes, dict):
        return {}
    return {k: v for k, v in analytes.items() if isinstance(v, dict)}


def _analyte_score(a0: dict, a2: dict) -> float | None:
    """How alike two analytes are, or None if they may not be paired.

    Two analytes may be paired when biomarker and specimen both agree, or when
    one of the two differs while unit and reference event agree. The second
    case lets a wrong biomarker or specimen be reported as one field error.
    """
    same = {f: a0.get(f) == a2.get(f) for f in ANALYTE_FIELDS}
    core = same["biomarker"] + same["specimen"]
    if core < 1 or (core == 1 and not (same["unit"] and same["reference_event"])):
        return None
    text = difflib.SequenceMatcher(
        None, str(a0.get("description", "")), str(a2.get("description", ""))
    ).ratio()
    return (
        100 * core
        + 10 * (same["unit"] + same["reference_event"])
        + 2 * same["gene_target"]
        + text
    )


def match_analytes(a0: dict, a2: dict):
    """Pair S0 analytes with S2 analytes.

    Returns (pairs, extra, missing): `pairs` maps an S0 key to an S2 key,
    `extra` lists unmatched S0 keys and `missing` unmatched S2 keys.
    """
    keys0, keys2 = list(a0), list(a2)
    pairs: dict = {}
    if keys0 and keys2:
        scores = np.zeros((len(keys0), len(keys2)))
        for i, k0 in enumerate(keys0):
            for j, k2 in enumerate(keys2):
                scores[i, j] = _analyte_score(a0[k0], a2[k2]) or 0.0
        for i, j in zip(*linear_sum_assignment(scores, maximize=True)):
            if scores[i, j] > 0:
                pairs[keys0[i]] = keys2[j]
    matched2 = set(pairs.values())
    return (
        pairs,
        [k for k in keys0 if k not in pairs],
        [k for k in keys2 if k not in matched2],
    )


# ----------------------------------------------------------------- measurements


def participants_of(dataset) -> list[dict]:
    """The dataset's participants, with anything that is not a mapping emptied."""
    participants = dataset.get("participants") if isinstance(dataset, dict) else None
    if not isinstance(participants, list):
        return []
    return [p if isinstance(p, dict) else {} for p in participants]


def measurements_of(participant: dict) -> list[dict]:
    measurements = participant.get("measurements")
    if not isinstance(measurements, list):
        return []
    return [m if isinstance(m, dict) else {} for m in measurements]


def by_analyte(participant: dict) -> dict:
    """{analyte key: [(index, time, value), ...]} for one participant."""
    groups: dict = {}
    for index, m in enumerate(measurements_of(participant)):
        groups.setdefault(m.get("analyte"), []).append(
            (index, m.get("time"), m.get("value"))
        )
    return groups


@dataclass
class MeasurementPairs:
    """The outcome of matching one participant's measurements for one analyte.

    Indices are positions in each participant's `measurements` list.
    """

    matched: list[tuple[int, int, bool]] = field(default_factory=list)
    shifted: list[tuple[int, int]] = field(default_factory=list)
    extra: list[int] = field(default_factory=list)
    missing: list[int] = field(default_factory=list)


def pair_measurements(m0: list, m2: list, unit, tol: Tolerances) -> MeasurementPairs:
    """Match two lists of (index, time, value).

    First, measurements whose times agree are paired, preferring pairs whose
    values agree, then the closest time, then the nearest value. `matched`
    holds these with a flag for whether the values agree. Then, among what is
    left, measurements whose values agree are paired on the closest time:
    these are the same reading at a different time (`shifted`). The rest are
    `extra` (S0 only) or `missing` (S2 only).
    """
    out = MeasurementPairs()
    used0: set = set()
    used2: set = set()

    candidates = []
    for i0, t0, v0 in m0:
        for i2, t2, v2 in m2:
            if times_agree(t0, t2, tol):
                agree = values_agree(v0, v2, unit, tol)
                candidates.append(
                    (
                        not agree,
                        time_distance(t0, t2),
                        value_distance(v0, v2),
                        i0,
                        i2,
                    )
                )
    for disagree, _, _, i0, i2 in sorted(candidates):
        if i0 not in used0 and i2 not in used2:
            used0.add(i0)
            used2.add(i2)
            out.matched.append((i0, i2, not disagree))

    candidates = []
    for i0, t0, v0 in m0:
        if i0 in used0:
            continue
        for i2, t2, v2 in m2:
            if i2 not in used2 and values_agree(v0, v2, unit, tol):
                candidates.append((time_distance(t0, t2), i0, i2))
    for _, i0, i2 in sorted(candidates):
        if i0 not in used0 and i2 not in used2:
            used0.add(i0)
            used2.add(i2)
            out.shifted.append((i0, i2))

    out.matched.sort()
    out.shifted.sort()
    out.extra = sorted(i for i, _, _ in m0 if i not in used0)
    out.missing = sorted(i for i, _, _ in m2 if i not in used2)
    return out


# ----------------------------------------------------------------- participants

_UNKNOWN_TIME = 1e15  # stands in for a non-numeric time; agrees only with itself


def _columns(participants: list[dict], key, codes: dict):
    """One analyte's measurements across participants, as arrays."""
    who, times, values, kinds = [], [], [], []
    for p_index, participant in enumerate(participants):
        for m in measurements_of(participant):
            if m.get("analyte") != key:
                continue
            time, value = m.get("time"), m.get("value")
            who.append(p_index)
            times.append(float(time) if is_number(time) else _UNKNOWN_TIME)
            if is_number(value):
                values.append(float(value))
                kinds.append(0)
            else:
                values.append(np.nan)
                kinds.append(codes.setdefault(repr(value), len(codes) + 1))
    return (
        np.array(who, dtype=np.int64),
        np.array(times, dtype=float),
        np.array(values, dtype=float),
        np.array(kinds, dtype=np.int64),
    )


def _agree_arrays(v0, k0, v2, k2, unit, tol: Tolerances):
    """`values_agree` for arrays of values and their kind codes."""
    numeric = (k0 == 0) & (k2 == 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        if unit == CT_UNIT:
            close = np.abs(v0 - v2) <= tol.ct + EPS
        else:
            close = (
                np.abs(v0 - v2)
                <= tol.relative * np.maximum(np.abs(v0), np.abs(v2)) + EPS
            )
            if tol.log10 is not None:
                positive = (v0 > 0) & (v2 > 0)
                logs = np.abs(np.log10(v0) - np.log10(v2)) <= tol.log10 + EPS
                close = np.where(positive, logs, close)
    return np.where(numeric, close, (k0 == k2) & (k0 != 0))


def _distinct_counts(who0, who2, items, shape):
    """Per participant pair, how many distinct `items` appear."""
    counts = np.zeros(shape, dtype=np.int64)
    if len(items):
        triples = np.unique(np.stack([who0, who2, items]), axis=1)
        np.add.at(counts, (triples[0], triples[1]), 1)
    return counts


def _overlap(p0: list[dict], p2: list[dict], k0, k2, unit, tol: Tolerances, codes):
    """For one analyte, measurements each participant pair shares.

    Returns two matrices: pairs agreeing on time and value, and pairs agreeing
    on time. A measurement is counted once however many partners it has, and
    the smaller of the two directions is taken, so the count never exceeds
    what one-to-one matching could give.
    """
    shape = (len(p0), len(p2))
    who0, t0, v0, c0 = _columns(p0, k0, codes)
    who2, t2, v2, c2 = _columns(p2, k2, codes)
    if not len(who0) or not len(who2):
        zeros = np.zeros(shape, dtype=np.int64)
        return zeros, zeros

    order = np.argsort(t2, kind="stable")
    sorted_t2 = t2[order]
    lo = np.searchsorted(sorted_t2, t0 - tol.time - EPS, side="left")
    hi = np.searchsorted(sorted_t2, t0 + tol.time + EPS, side="right")
    width = hi - lo
    left = np.repeat(np.arange(len(t0)), width)
    starts = np.repeat(lo - np.concatenate(([0], np.cumsum(width)[:-1])), width)
    right = order[np.arange(width.sum()) + starts]

    def count(keep):
        a, b = left[keep], right[keep]
        return np.minimum(
            _distinct_counts(who0[a], who2[b], a, shape),
            _distinct_counts(who0[a], who2[b], b, shape),
        )

    agree = _agree_arrays(v0[left], c0[left], v2[right], c2[right], unit, tol)
    return count(agree), count(np.ones(len(left), dtype=bool))


def _attributes(participant: dict) -> dict:
    attributes = participant.get("attributes")
    return attributes if isinstance(attributes, dict) else {}


def _attribute_agreement(p0: list[dict], p2: list[dict]) -> np.ndarray:
    """Share of attribute fields on which each participant pair agrees."""
    out = np.zeros((len(p0), len(p2)))
    attrs0 = [_attributes(p) for p in p0]
    attrs2 = [_attributes(p) for p in p2]
    for i, a0 in enumerate(attrs0):
        for j, a2 in enumerate(attrs2):
            fields = a0.keys() | a2.keys()
            if fields:
                same = sum(1 for f in fields if f in a0 and f in a2 and a0[f] == a2[f])
                out[i, j] = same / len(fields)
    return out


def match_participants(
    p0: list[dict], p2: list[dict], analyte_pairs: dict, a2: dict, tol: Tolerances
):
    """Pair S0 participants with S2 participants.

    The pairing maximizes the number of shared (analyte, time, value) triples,
    then shared (analyte, time) pairs, then agreement on attributes. Two
    participants may be paired only if they share at least one triple, or if
    at least half of the measurements of the larger one agree on time. The
    second case keeps a participant whose values are all wrong from being
    reported as one missing and one extra participant.

    Returns (pairs, extra, missing) with `pairs` as {S0 index: S2 index}.
    """
    shape = (len(p0), len(p2))
    pairs: dict = {}
    if all(shape):
        triples = np.zeros(shape, dtype=np.int64)
        times = np.zeros(shape, dtype=np.int64)
        codes: dict = {}
        for k0, k2 in analyte_pairs.items():
            shared, on_time = _overlap(p0, p2, k0, k2, a2[k2].get("unit"), tol, codes)
            triples += shared
            times += on_time

        def size(participants, keys):
            return np.array(
                [
                    sum(1 for m in measurements_of(p) if m.get("analyte") in keys)
                    for p in participants
                ]
            )

        n0 = size(p0, set(analyte_pairs))
        n2 = size(p2, set(analyte_pairs.values()))
        larger = np.maximum(n0[:, None], n2[None, :])
        attributes = _attribute_agreement(p0, p2)
        allowed = (
            (triples >= 1) | ((times >= 1) & (2 * times >= larger)) | (larger == 0)
        )
        scores = np.where(
            allowed, 1e8 * triples + 1e4 * times + 1e3 * attributes + 1, 0
        )
        for i, j in zip(*linear_sum_assignment(scores, maximize=True)):
            if allowed[i, j]:
                pairs[int(i)] = int(j)
    matched2 = set(pairs.values())
    return (
        pairs,
        [i for i in range(len(p0)) if i not in pairs],
        [j for j in range(len(p2)) if j not in matched2],
    )
