"""Hand-checkable cases for graph_compare (ticket #2, oracle level 1: tiny graphs whose differences are known by construction)."""
import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "functions"))
from graph_compare import compare  # noqa: E402

BASE = {
    "nodes": [{"id": "card:A", "PMID": ["2", "1"]}, {"id": "ncbi:B", "name": "x"}, {"id": "resfinder:C"}],
    "links": [
        {"source": "card:A", "target": "ncbi:B", "type": "RBH", "identity": 0.95, "coverage": 1.0},
        {"source": "ncbi:B", "target": "card:A", "type": "RBH", "identity": 0.95, "coverage": 1.0},
        {"source": "card:A", "target": "resfinder:C", "type": "OWH", "identity": 0.80, "coverage": 0.9},
    ],
}


def variant(**changes):
    g = copy.deepcopy(BASE)
    for key, fn in changes.items():
        fn(g)
    return g


def test_identical_graphs_pass():
    r = compare(BASE, copy.deepcopy(BASE))
    assert r["ok"] and not r["findings"]


def test_node_order_and_pmid_order_ignored():
    g = variant(a=lambda g: g["nodes"].reverse(), b=lambda g: g["nodes"][0].update(PMID=["1", "2"]))
    assert compare(BASE, g)["ok"]


def test_missing_node_fails():
    g = variant(a=lambda g: g["nodes"].pop())
    r = compare(BASE, g)
    assert not r["ok"] and r["nodes_lost"] == ["resfinder:C"]


def test_extra_edge_and_lost_edge_fail():
    g = variant(a=lambda g: g["links"].pop(),
                b=lambda g: g["links"].append({"source": "resfinder:C", "target": "card:A", "type": "OWH", "identity": 0.5, "coverage": 0.5}))
    r = compare(BASE, g)
    assert not r["ok"]
    assert r["edges_lost"] == [("card:A", "resfinder:C")] and r["edges_gained"] == [("resfinder:C", "card:A")]


def test_type_change_fails():
    g = variant(a=lambda g: g["links"][2].update(type="RBH"))
    r = compare(BASE, g)
    assert not r["ok"] and r["type_changed"] == [(("card:A", "resfinder:C"), "OWH", "RBH")]


def test_tolerance_is_0_001_on_identity_and_coverage():
    ok = variant(a=lambda g: g["links"][0].update(identity=0.9505, coverage=0.9995))
    bad_id = variant(a=lambda g: g["links"][0].update(identity=0.952))
    bad_cov = variant(a=lambda g: g["links"][2].update(coverage=0.897))
    assert compare(BASE, ok)["ok"]
    assert not compare(BASE, bad_id)["ok"] and not compare(BASE, bad_cov)["ok"]


def test_attribute_differences_are_reported_not_failed():
    g = variant(a=lambda g: g["nodes"][1].update(name="y"))
    r = compare(BASE, g)
    assert r["ok"] and r["attribute_diffs"] == [("ncbi:B", "name", "x", "y")]
