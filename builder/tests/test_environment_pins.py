"""Pins for the MMseqs2 upgrade comparison (ticket #11, ADR 0003 step 2).

The baseline environment keeps MMseqs2 13.45111 (the tool the 2021 graph was built with). The upgrade
environment must differ from it in the MMseqs2 pin only, so any change in the graph is attributable to the tool.
"""
from pathlib import Path

BUILDER = Path(__file__).resolve().parents[1]


def pins(name):
    out = {}
    for line in (BUILDER / name).read_text().splitlines():
        line = line.strip()
        if line.startswith("- ") and "=" in line:
            tool, _, version = line[2:].partition("=")
            out[tool] = version
    return out


def test_baseline_keeps_mmseqs2_13():
    assert pins("environment.yml")["mmseqs2"] == "13.45111"


def test_upgrade_environment_differs_only_in_mmseqs2():
    base, new = pins("environment.yml"), pins("environment-mmseqs18.yml")
    assert new["mmseqs2"] == "18.8cc5c"
    assert {k: v for k, v in new.items() if k != "mmseqs2"} == {k: v for k, v in base.items() if k != "mmseqs2"}
