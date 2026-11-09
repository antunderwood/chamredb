"""parse_card_metadata.py must give the same JSON whatever PYTHONHASHSEED (ticket #18).

Hand-built mini ontology (oracle level 1, expected values worked out by hand below) plus a byte-identity
check over several hash seeds (level 5, invariant). Found in #9 (PMID order) and #2 (phenotype inheritance
depended on term order).
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pronto")
SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "parse_card_metadata.py"

OBO = """format-version: 1.2
ontology: antibiotic_resistance

[Typedef]
id: confers_resistance_to_antibiotic
name: confers_resistance_to_antibiotic

[Term]
id: ARO:3100001
name: drugX
def: "x" []

[Term]
id: ARO:3100002
name: drugY
def: "x" []

[Term]
id: ARO:3100003
name: drugZ
def: "x" []

[Term]
id: ARO:3000001
name: child
def: "child" [PMID:3, PMID:1, PMID:2]
is_a: ARO:3000002 ! middle

[Term]
id: ARO:3000002
name: middle
def: "x" []
is_a: ARO:3000003 ! root

[Term]
id: ARO:3000003
name: root
def: "x" []
relationship: confers_resistance_to_antibiotic ARO:3100001 ! drugX

[Term]
id: ARO:3000004
name: twoparents
def: "x" []
is_a: ARO:3000006 ! p6
is_a: ARO:3000005 ! p5

[Term]
id: ARO:3000005
name: p5
def: "x" []
relationship: confers_resistance_to_antibiotic ARO:3100002 ! drugY

[Term]
id: ARO:3000006
name: p6
def: "x" []
relationship: confers_resistance_to_antibiotic ARO:3100003 ! drugZ
"""


def run(tmp_path, seed):
    src = tmp_path / "mini.obo"
    src.write_text(OBO)
    out = tmp_path / f"out{seed}.json"
    env = {**os.environ, "PYTHONHASHSEED": str(seed)}
    subprocess.run([sys.executable, str(SCRIPT), "--input", str(src), "--output", str(out)],
                   check=True, env=env, capture_output=True)
    return out.read_bytes()


def test_output_is_identical_for_all_hash_seeds(tmp_path):
    outs = {run(tmp_path, seed) for seed in range(1, 11)}
    assert len(outs) == 1


def test_hand_worked_expectations(tmp_path):
    m = json.loads(run(tmp_path, 1))
    assert m["ARO:3000001"]["PMID"] == ["1", "2", "3"]          # sorted
    assert m["ARO:3000004"]["is_a"] == "ARO:3000005:p5,ARO:3000006:p6"  # parents sorted
    # terms are processed in sorted id order and inherit from the first parent as it stands at that moment:
    # child (3000001) is handled before middle (3000002) has inherited drugX, so it gets nothing
    assert "phenotype" not in m["ARO:3000001"]
    assert m["ARO:3000002"]["phenotype"] == "confers resistance to antibiotic: drugX"
    assert m["ARO:3000003"]["phenotype"] == "confers resistance to antibiotic: drugX"
    assert m["ARO:3000004"]["phenotype"] == "confers resistance to antibiotic: drugY"  # first parent p5
