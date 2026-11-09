"""Spec #1 predicates 1 to 4 as tests (ticket #3).

Each predicate is shown to fail on a deliberately broken tiny graph (oracle level 1, by construction), then
run on the shipped 2021 graph and committed FASTAs (level 5, invariants). Needs no MMseqs2 and no network.
"""
import copy
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "builder" / "scripts" / "functions"))
import graph_invariants as gi  # noqa: E402

SHIPPED = ROOT / "chamredb" / "data"

GOOD = {
    "nodes": [{"id": "card:ARO:1"}, {"id": "ncbi:WP_1.1"}, {"id": "resfinder:dfrB1"}],
    "links": [
        {"source": "card:ARO:1", "target": "ncbi:WP_1.1", "type": "RBH", "identity": 1.0, "coverage": 1.0},
        {"source": "ncbi:WP_1.1", "target": "card:ARO:1", "type": "RBH", "identity": 1.0, "coverage": 1.0},
        {"source": "card:ARO:1", "target": "resfinder:dfrB1", "type": "OWH", "identity": 0.5, "coverage": 0.25},
    ],
}


def broken(fn):
    g = copy.deepcopy(GOOD)
    fn(g)
    return g


def test_good_graph_has_no_violations():
    assert gi.check_graph(GOOD) == {}


@pytest.mark.parametrize("name,fn", [
    ("p1_edge_endpoint_not_a_node", lambda g: g["nodes"].pop()),
    ("p1_node_id_without_database", lambda g: g["nodes"].append({"id": "nodb"})),
    ("p1_unknown_database", lambda g: g["nodes"].append({"id": "foo:1"})),
    ("p1_edge_within_one_database", lambda g: (g["nodes"].append({"id": "card:ARO:2"}),
        g["links"].append({"source": "card:ARO:1", "target": "card:ARO:2", "type": "OWH", "identity": 1, "coverage": 1}))),
    ("p2_rbh_not_reciprocal", lambda g: g["links"].pop(1)),
    ("p2_reverse_edge_is_not_rbh", lambda g: g["links"][1].update(type="OWH")),
    ("p3_identity_zero", lambda g: g["links"][0].update(identity=0)),
    ("p3_identity_above_one", lambda g: g["links"][0].update(identity=1.01)),
    ("p3_coverage_zero", lambda g: g["links"][2].update(coverage=0)),
])
def test_each_graph_predicate_fails_on_a_broken_graph(name, fn):
    assert name[:2] in {k[:2] for k in gi.check_graph(broken(fn))}, name


def test_duplicate_fasta_ids_are_detected(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(">a desc\nAC\n>b\nAC\n>a other\nGG\n")
    assert gi.duplicate_fasta_ids(f) == ["a"]
    f.write_text(">a\nAC\n>b\nAC\n")
    assert gi.duplicate_fasta_ids(f) == []


# ---- the shipped 2021 graph and committed FASTAs ----

@pytest.fixture(scope="module")
def shipped_graph():
    return json.loads((SHIPPED / "graph.json").read_text())


@pytest.mark.parametrize("predicate", ["p1", "p2", "p3"])
def test_shipped_graph_predicate(shipped_graph, predicate):
    violations = {k: v for k, v in gi.check_graph(shipped_graph).items() if k.startswith(predicate)}
    assert not violations, {k: v[:5] for k, v in violations.items()}


RESFINDER_FINDING = ("FINDING (ticket #3): the committed resfinder.nucl.fasta (cat of the per-class *.fsa) has 5 "
                     "duplicate ids (aac(6')-Ib-cr_1_DQ303918, aac(6')-Ib-cr_2_EF636461, cfr(B)_3_KR610408, "
                     "cfr_1_AM408573, cfr_2_AJ879565): genes listed in two class files. Parsed protein FASTA is unique.")


def test_coverage_above_one_is_documented_behaviour(shipped_graph):
    # Decision 2026-11-09: coverage = alignment length over query end and may exceed 1 (gapped alignments).
    # Characterisation (level 6) of the 2021 graph: 729 of 33,632 edges, max 1.198. See GLOSSARY.md.
    over = [e["coverage"] for e in shipped_graph["links"] if e["coverage"] > 1]
    assert len(over) == 729 and max(over) < 1.2


@pytest.mark.parametrize("name", [
    "card.protein.fasta", "ncbi.protein.fasta", "resfinder.protein.fasta", "ncbi.protein.raw.fasta",
    pytest.param("resfinder.nucl.fasta", marks=pytest.mark.xfail(strict=True, reason=RESFINDER_FINDING))])
def test_shipped_fasta_ids_unique_p4(name):
    assert gi.duplicate_fasta_ids(SHIPPED / "db_fastas" / name) == []
