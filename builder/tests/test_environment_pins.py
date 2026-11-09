"""Pins for the builder environments (ticket #11, ADR 0003).

environment.yml uses the current MMseqs2 (adopted 2026-11-09 after reviewing the diff report on the 2021 inputs).
environment-mmseqs13.yml reproduces the shipped 2021 graph and must differ from it in the MMseqs2 pin only.
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


def test_main_environment_uses_mmseqs2_18():
    assert pins("environment.yml")["mmseqs2"] == "18.8cc5c"


def test_reproduction_environment_differs_only_in_mmseqs2():
    main, old = pins("environment.yml"), pins("environment-mmseqs13.yml")
    assert old["mmseqs2"] == "13.45111"
    assert {k: v for k, v in old.items() if k != "mmseqs2"} == {k: v for k, v in main.items() if k != "mmseqs2"}
