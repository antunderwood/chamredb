"""graph_diff on a hand-built pair of small graphs (ticket #7, oracle level 1).

Every expected value below is worked out by hand from the two graphs in this file.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "functions"))
from graph_diff import diff, report  # noqa: E402


def node(i, name):
    return {"id": i, "name": name}


def edge(s, t, ty, ident, cov):
    return {"source": s, "target": t, "type": ty, "identity": ident, "coverage": cov}


OLD = {
    "nodes": [node("card:ARO:1", "a"), node("card:ARO:2", "b"), node("card:ARO:3", "c"),
              node("ncbi:N1", "x"), node("ncbi:N2", "y"), node("resfinder:r1", "r")],
    "links": [
        edge("card:ARO:1", "ncbi:N1", "RBH", 0.95, 1.0), edge("ncbi:N1", "card:ARO:1", "RBH", 0.95, 1.0),
        edge("card:ARO:2", "resfinder:r1", "OWH", 0.80, 0.9),
        edge("ncbi:N2", "resfinder:r1", "OWH", 0.70, 0.9),
        edge("card:ARO:3", "ncbi:N1", "OWH", 0.85, 0.9),
    ],
}
NEW = {
    "nodes": [node("card:ARO:1", "a"), node("card:ARO:2", "B2"), node("card:ARO:30", "c"),
              node("ncbi:N1", "x"), node("ncbi:N9", "z"), node("resfinder:r1", "r")],
    "links": [
        edge("card:ARO:1", "ncbi:N1", "RBH", 0.95, 1.0), edge("ncbi:N1", "card:ARO:1", "RBH", 0.97, 1.0),
        edge("card:ARO:2", "resfinder:r1", "RBH", 0.80, 0.9), edge("resfinder:r1", "card:ARO:2", "RBH", 0.80, 0.9),
        edge("card:ARO:30", "ncbi:N1", "OWH", 0.85, 0.9),
        edge("ncbi:N9", "resfinder:r1", "OWH", 0.60, 0.9),
    ],
}


def test_nodes_and_edges_gained_and_lost():
    d = diff(OLD, NEW)
    assert d["nodes_gained"] == ["card:ARO:30", "ncbi:N9"]
    assert d["nodes_lost"] == ["card:ARO:3", "ncbi:N2"]
    assert d["edges_gained"] == [("card:ARO:30", "ncbi:N1"), ("ncbi:N9", "resfinder:r1"), ("resfinder:r1", "card:ARO:2")]
    assert d["edges_lost"] == [("card:ARO:3", "ncbi:N1"), ("ncbi:N2", "resfinder:r1")]


def test_type_and_identity_changes():
    d = diff(OLD, NEW)
    assert d["type_changed"] == [(("card:ARO:2", "resfinder:r1"), "OWH", "RBH")]
    assert d["identity_off"] == [(("ncbi:N1", "card:ARO:1"), 0.95, 0.97)]
    assert d["coverage_off"] == []


def test_id_continuity_kept_reassigned_absent_with_successor():
    c = diff(OLD, NEW)["continuity"]
    assert c["kept"] == ["card:ARO:1", "ncbi:N1", "resfinder:r1"]
    assert c["reassigned"] == [("card:ARO:2", "b", "B2")]
    assert c["absent"] == [("card:ARO:3", ["card:ARO:30"]), ("ncbi:N2", [])]
    assert c["by_database"] == {"card": {"kept": 1, "reassigned": 1, "absent": 1},
                                "ncbi": {"kept": 1, "reassigned": 0, "absent": 1},
                                "resfinder": {"kept": 1, "reassigned": 0, "absent": 0}}


def test_successor_must_be_in_the_same_database():
    new = {"nodes": [node("ncbi:N7", "c")], "links": []}
    old = {"nodes": [node("card:ARO:3", "c")], "links": []}
    assert diff(old, new)["continuity"]["absent"] == [("card:ARO:3", [])]


def test_identical_graphs_report_no_changes():
    d = diff(OLD, OLD)
    assert d["ok"] and not d["continuity"]["absent"] and not d["continuity"]["reassigned"]


def test_report_states_the_counts():
    text = report(diff(OLD, NEW))
    for line in ("Nodes gained: 2", "Nodes lost: 2", "Edges gained: 3", "Edges lost: 2", "Edges with changed type: 1",
                 "IDs absent in the new graph: 2", "IDs reassigned (same id, different name): 1",
                 "card:ARO:3 -> card:ARO:30", "ncbi:N2 -> (no successor)", "card:ARO:2: b -> B2"):
        assert line in text, line
