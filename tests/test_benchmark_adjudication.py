"""The adjudication packets are blinded, follow the plan's sampling and rebuild exactly.

The discrepancies come from drafts made by perturbing copies of released
datasets, as in test_benchmark_compare.py.
"""

import copy
import csv
import math
import re
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark"))
import adjudication as adj  # noqa: E402
import align  # noqa: E402
import discrepancies as disc  # noqa: E402
import match_findings as mf  # noqa: E402

MANIFEST = disc.read_manifest()
MANUAL, AI = "tan2021early", "antinori2022epidemiological"
# Words that would say which version a value came from, or which way a
# difference runs.
GIVES_AWAY = re.compile(
    r"\b(s0|s2|draft|reference|agent|omission|extra|missing|released?)\b", re.I
)


def load(study):
    path = REPO_ROOT / "data" / study / f"{study}.yaml"
    return yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.CSafeLoader)


def first_numeric(dataset, skip=0):
    """(participant, measurement) of a numeric reading whose value is unique."""
    seen = 0
    for i, participant in enumerate(dataset["participants"]):
        ms = participant["measurements"]
        for j, m in enumerate(ms):
            if not align.is_number(m["value"]) or j + 1 >= len(ms):
                continue
            if sum(o["value"] == m["value"] for o in ms) > 1:
                continue
            if seen == skip:
                return i, j
            seen += 1
    raise AssertionError("no numeric reading")


def manual_draft():
    """A draft of the manual-era study with one error of each common kind."""
    reference = load(MANUAL)
    draft = copy.deepcopy(reference)
    i, j = first_numeric(reference)
    ms = draft["participants"][i]["measurements"]
    ms[j]["value"] = ms[j]["value"] + 9  # wrong value
    del ms[j + 1]  # omission
    ms.append({**ms[0], "time": 999})  # extra
    draft["participants"][i]["attributes"]["age"] += 1  # wrong attribute value
    del draft["participants"][i + 1]["attributes"]["sex"]  # missing attribute
    key = next(iter(draft["analytes"]))
    draft["analytes"][key]["unit"] = "gc/mL"  # wrong unit
    draft["analytes"][key]["gene_target"] = "reworded"  # listed, not scored
    del draft["participants"][-1]  # missing participant
    return draft, reference, (i, j)


def ai_draft():
    reference = load(AI)
    draft = copy.deepcopy(reference)
    i, j = first_numeric(reference)
    draft["participants"][i]["measurements"][j]["value"] *= 3
    return draft, reference


def comparison():
    """Discrepancy rows and Study objects for the two drafts and a failed one."""
    rows, studies = [], {}
    d0, r0, _ = manual_draft()
    d1, r1 = ai_draft()
    for study, draft, reference in ((MANUAL, d0, r0), (AI, d1, r1)):
        rows += disc.compare(study, draft, reference, align.REPORTED).rows
        studies[study] = adj.Study(draft, reference, align.REPORTED)
    failed = "woelfel2020virological"
    rows += disc.no_draft(failed, load(failed), "no draft").rows
    return rows, studies


def build(seed=adj.SEED, findings=None, finding_rows=None):
    rows, studies = comparison()
    return (
        adj.build(rows, finding_rows or [], findings or {}, studies, MANIFEST, seed),
        rows,
    )


def by_type(packet, key):
    items = {item["item_id"]: item for item in packet}
    out = {}
    for k in key:
        if k["packet"] == "discrepancies":
            item = items[k["item_id"]]
            sides = {k["A_is"]: item["A"], k["B_is"]: item["B"]}
            out.setdefault(k["type"], []).append((item, sides))
    return out


def test_the_packet_does_not_say_which_version_is_which():
    (packet, _, key), _ = build()
    assert packet
    for item in packet:
        assert list(item) == adj.DISCREPANCY_COLUMNS
        for column in ("difference", "where", "A", "B"):
            assert not GIVES_AWAY.search(item[column]), (column, item[column])
    # Omissions and extras read the same way.
    found = by_type(packet, key)
    assert {i["difference"] for i, _ in found["omission"]} == {
        i["difference"] for i, _ in found["extra"]
    }


def test_the_key_recovers_each_injected_error():
    (packet, _, key), _ = build()
    found = by_type(packet, key)
    draft, reference, (i, j) = manual_draft()
    wrong = reference["participants"][i]["measurements"][j]["value"] + 9
    values = [sides["S0"] for _, sides in found["wrong value"]]
    assert any(f"value={adj.text(wrong)}" in v for v in values)
    assert all(sides["S0"] == adj.ABSENT for _, sides in found["omission"])
    assert all(sides["S2"] == adj.ABSENT for _, sides in found["extra"])
    assert all(sides["S0"] == adj.ABSENT for _, sides in found["missing"])
    assert all(sides["S0"] == adj.NOT_GIVEN for _, sides in found["missing attribute"])
    assert [sides for _, sides in found["wrong unit"]] == [
        {
            "S0": "gc/mL",
            "S2": reference["analytes"][next(iter(reference["analytes"]))]["unit"],
        }
    ]
    # An extra measurement is located on its participant's partner.
    for item, _ in found["extra"]:
        assert item["where"].startswith("participant at position")


def test_only_scored_comparable_rows_are_adjudicated():
    (packet, _, key), rows = build()
    eligible = [
        r for r in rows if r["scored"] == "yes" and r["type"] not in adj.NOT_ADJUDICATED
    ]
    assert sorted(k["ref_id"] for k in key) == sorted(r["id"] for r in eligible)
    assert {r["type"] for r in rows} - {k["type"] for k in key} == {
        "wrong gene_target",
        "no draft",
    }


def test_samples_follow_the_plan():
    era_of = {"m": "manual", "a": "ai_assisted"}

    def rows(study, n):
        return [{"study_id": study, "_ref_id": f"{study}{i:04d}"} for i in range(n)]

    def count(chosen):
        return {
            e: sum(era_of[r["study_id"]] == e for r in chosen)
            for e in ("manual", "ai_assisted")
        }

    take_all = lambda r: True  # noqa: E731
    big = rows("m", 450) + rows("a", 500)
    assert count(adj.select(big, era_of, 1, "D", take_all)) == {
        "manual": 300,
        "ai_assisted": 200,
    }
    assert count(adj.select(big, era_of, 1, "F", take_all)) == {
        "manual": 450,
        "ai_assisted": 200,
    }
    small = rows("m", 120) + rows("a", 50)
    assert count(adj.select(small, era_of, 1, "D", take_all)) == {
        "manual": 120,
        "ai_assisted": 50,
    }


def test_the_recheck_is_a_tenth_of_each_era_rounded_up():
    items = [
        {"item_id": f"D{i:04d}", "era": "manual" if i < 34 else "ai_assisted"}
        for i in range(37)
    ]
    chosen = adj.recheck_ids(items, adj.SEED, "D")
    manual = {f"D{i:04d}" for i in range(34)}
    assert len(chosen & manual) == math.ceil(3.4) == 4
    assert len(chosen - manual) == 1
    assert chosen == adj.recheck_ids(items, adj.SEED, "D")


def test_the_same_seed_gives_the_same_packets_and_another_seed_does_not():
    first, _ = build()
    again, _ = build()
    other, _ = build(seed=adj.SEED + 1)
    assert first == again
    assert [k["A_is"] for k in first[2]] != [k["A_is"] for k in other[2]] or [
        k["ref_id"] for k in first[2]
    ] != [k["ref_id"] for k in other[2]]


def test_a_and_b_are_drawn_with_even_odds():
    shares = []
    for seed in range(20):
        (_, _, key), _ = build(seed=seed)
        shares += [k["A_is"] == "S0" for k in key]
    assert 0.4 < sum(shares) / len(shares) < 0.6


def test_only_findings_that_hit_nothing_are_judged():
    rows, studies = comparison()
    draft = studies[MANUAL].s0
    wrong = next(
        r for r in rows if r["study_id"] == MANUAL and r["type"] == "wrong value"
    )
    # A participant the draft leaves untouched.
    i, j = manual_draft()[2][0] + 2, 0
    quiet = f"/participants/{i}/measurements/{j}/value"
    findings = [
        {
            "id": "hit",
            "check": "analyte.value_vs_paper",
            "severity": "error",
            "path": wrong["s0_path"],
            "message": "m1",
        },
        {
            "id": "unseen",
            "check": "analyte.value_vs_paper",
            "severity": "warning",
            "path": quiet,
            "message": "m2",
            "confidence": 0.9,
        },
        {
            "id": "gap",
            "check": "schema.missing_field",
            "severity": "info",
            "path": "/analytes",
            "message": "m3",
        },
    ]
    finding_rows = mf.match(
        MANUAL, findings, [r for r in rows if r["study_id"] == MANUAL]
    )
    stored = {(MANUAL, f["id"]): f for f in findings}
    _, packet, key = adj.build(rows, finding_rows, stored, studies, MANIFEST)
    assert [item["message"] for item in packet] == ["m2"]
    assert packet[0]["draft_content"] == adj.text(
        draft["participants"][i]["measurements"][j]["value"]
    )
    assert packet[0]["confidence"] == "0.9"
    assert [k["ref_id"] for k in key if k["packet"] == "findings"] == [
        f"{MANUAL}:unseen"
    ]


def test_the_command_line_writes_packets_that_rebuild_exactly(tmp_path):
    draft, _, _ = manual_draft()
    study_dir = tmp_path / "run" / MANUAL
    study_dir.mkdir(parents=True)
    (study_dir / "s0.yaml").write_text(yaml.safe_dump(draft), encoding="utf-8")
    (study_dir / "findings.yaml").write_text(
        yaml.safe_dump(
            {
                "findings": [
                    {
                        "id": "x",
                        "check": "attr.value_vs_paper",
                        "severity": "warning",
                        "path": "/participants/0/attributes/sex",
                        "message": "m",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    cmp_dir = tmp_path / "cmp"
    assert (
        disc.main([str(tmp_path / "run"), "--out", str(cmp_dir), "--ref", "HEAD"]) == 0
    )
    assert mf.main([str(tmp_path / "run"), "--out", str(cmp_dir)]) == 0

    outputs = {}
    for name in ("a", "b"):
        out = tmp_path / name
        argv = [
            str(tmp_path / "run"),
            "--comparison",
            str(cmp_dir),
            "--out",
            str(out),
            "--ref",
            "HEAD",
        ]
        assert adj.main(argv) == 0
        outputs[name] = {p.name: p.read_bytes() for p in sorted(out.iterdir())}
    assert outputs["a"] == outputs["b"]
    assert set(outputs["a"]) == {
        "README.txt",
        "discrepancy_packet.csv",
        "findings_packet.csv",
        "key.csv",
        "packet_info.json",
        "recheck_discrepancy_packet.csv",
        "recheck_findings_packet.csv",
    }

    def read(name):
        return list(csv.DictReader((tmp_path / "a" / name).open(encoding="utf-8")))

    key = {k["item_id"]: k for k in read("key.csv")}
    for packet, recheck in (
        ("discrepancy_packet.csv", "recheck_discrepancy_packet.csv"),
        ("findings_packet.csv", "recheck_findings_packet.csv"),
    ):
        whole = {item["item_id"]: item for item in read(packet)}
        again = read(recheck)
        assert again and all(whole[item["item_id"]] == item for item in again)
        assert {i["item_id"] for i in again} == {
            i for i in whole if key[i]["recheck"] == "yes"
        }
        assert all(item["outcome"] == "" for item in whole.values())
    assert chr(0x2014) not in outputs["a"]["README.txt"].decode("utf-8")  # em dash
