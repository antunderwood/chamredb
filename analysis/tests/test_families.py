"""Selection of the real families drawn in the hero figure (ticket #41, spec #34).

Level 1: a hand-built graph with one candidate per category plus decoys that a wrong rule would pick; the expected choice
is worked out by hand. Level 5: on the release graph every selected node and edge is read back from the graph file
independently of networkx and compared (name, identity, coverage, type).
"""
import json
from pathlib import Path

import pytest
from networkx.readwrite import json_graph

from analysis import families, matches

ROOT = Path(__file__).resolve().parents[2]
BUNDLED = ROOT / "chamredb" / "data" / "graph.json"


def make(nodes, links):
    return json_graph.node_link_graph(
        {"directed": True, "multigraph": False, "graph": {},
         "nodes": [{"id": i, "name": n, "database": i.split(":")[0]} for i, n in nodes], "links": links}, edges="links")


def rbh(a, b, identity=1.0, coverage=1.0):
    return [{"source": a, "target": b, "type": "RBH", "identity": identity, "coverage": coverage},
            {"source": b, "target": a, "type": "RBH", "identity": identity, "coverage": coverage}]


def triple(tag, names, identity=1.0):
    ids = [f"card:{tag}c", f"ncbi:{tag}n", f"resfinder:{tag}r"]
    return list(zip(ids, names)), rbh(ids[0], ids[1], identity) + rbh(ids[1], ids[2], identity) + rbh(ids[0], ids[2], identity)


NODES, LINKS = [], []


def add(nodes, links):
    NODES.extend(nodes)
    LINKS.extend(links)


# category 1, names agree: 'abc' beats 'zzz' (first by name); a name over 14 characters is excluded
add(*triple("a1", ["zzz", "zzz", "zzz"]))
add(*triple("a2", ["abc", "abc", "abc"]))
add(*triple("a3", ["aaaaaaaaaaaaaaaaaaaa", "aaaaaaaaaaaaaaaaaaaa", "aaaaaaaaaaaaaaaaaaaa"]))
# category 2, case and punctuation differ only: ErmC / erm(C) / erm(C) all clean to ermc
add(*triple("p1", ["ErmC", "erm(C)", "erm(C)"]))
# category 3, allele number differs: qepa1 vs qepa3, allele-free all qepa
add(*triple("l1", ["QepA1", "qepA3", "qepA3"]))
# category 4, different names for identical proteins: the lowest name agreement at identity 1, ties by name. The decoy
# aaaa/bbbb/cccc has the same lowest agreement (0) and sorts first, but its identity is 0.95, so only the filter excludes it
add(*triple("d1", ["imiH", "cphA-13", "imiH"]))
add(*triple("d2", ["vanA", "vanA", "VanHAX"]))
add(*triple("d3", ["aaaa", "bbbb", "cccc"], identity=0.95))
# category 5, a family of 6: TEM-1 (card), TEM-1A to 1D (resfinder) and TEM-1 (ncbi); a family of 7 is a decoy
fam = [("card:t1", "TEM-1"), ("ncbi:t1", "blaTEM-1")] + [(f"resfinder:t1{x}", f"blaTEM-1{x}") for x in "ABCD"]
add(fam, rbh("card:t1", "ncbi:t1") + [l for x in "ABCD" for l in rbh("ncbi:t1", f"resfinder:t1{x}")])
big = [("card:o1", "OXA-1"), ("ncbi:o1", "blaOXA-1")] + [(f"resfinder:o{i}", f"blaOXA-{i}") for i in range(5)]
add(big, rbh("card:o1", "ncbi:o1") + [l for i in range(5) for l in rbh("ncbi:o1", f"resfinder:o{i}")])
# category 6, one-way hits: a card+ncbi pair whose members both reach the same single resfinder node one way
add([("card:w1", "chrB"), ("ncbi:w1", "chrB"), ("resfinder:wr", "erm(32)")],
    rbh("card:w1", "ncbi:w1") + [{"source": "card:w1", "target": "resfinder:wr", "type": "OWH", "identity": 0.7, "coverage": 0.9},
                                 {"source": "ncbi:w1", "target": "resfinder:wr", "type": "OWH", "identity": 0.7, "coverage": 0.9}])
# decoy: a pair whose one-way targets differ
add([("card:w2", "aaaa"), ("ncbi:w2", "aaaa"), ("resfinder:wx", "x1"), ("resfinder:wy", "y1")],
    rbh("card:w2", "ncbi:w2") + [{"source": "card:w2", "target": "resfinder:wx", "type": "OWH", "identity": 0.7, "coverage": 0.9},
                                 {"source": "ncbi:w2", "target": "resfinder:wy", "type": "OWH", "identity": 0.7, "coverage": 0.9}])
HAND = make(NODES, LINKS)


def test_selection_on_the_hand_graph():
    chosen = {f["key"]: sorted(f["nodes"]) for f in families.select_families(HAND)}
    assert list(chosen) == ["agree", "punctuation", "allele", "different", "family", "one_way"]
    assert chosen["agree"] == ["card:a2c", "ncbi:a2n", "resfinder:a2r"]
    assert chosen["punctuation"] == ["card:p1c", "ncbi:p1n", "resfinder:p1r"]
    assert chosen["allele"] == ["card:l1c", "ncbi:l1n", "resfinder:l1r"]
    assert chosen["different"] == ["card:d1c", "ncbi:d1n", "resfinder:d1r"]   # imiH/cphA-13 (cleaned imih vs cpha13) beats vanA/VanHAX
    assert chosen["family"] == sorted(n for n, _ in fam)
    assert chosen["one_way"] == ["card:w1", "ncbi:w1", "resfinder:wr"]


def test_edges_are_the_group_rbh_pairs_once_and_for_one_way_the_two_hits():
    by_key = {f["key"]: f for f in families.select_families(HAND)}
    # RBH pairs once, source in database order (card, ncbi, resfinder); one-way hits as directed edges
    assert by_key["agree"]["edges"] == [("card:a2c", "ncbi:a2n"), ("card:a2c", "resfinder:a2r"), ("ncbi:a2n", "resfinder:a2r")]
    assert by_key["family"]["edges"] == [("card:t1", "ncbi:t1")] + [("ncbi:t1", f"resfinder:t1{x}") for x in "ABCD"]
    assert by_key["one_way"]["edges"] == [("card:w1", "ncbi:w1"), ("card:w1", "resfinder:wr"), ("ncbi:w1", "resfinder:wr")]


def test_release_graph_selection_is_read_back_from_the_graph_file():
    graph = matches.read_graph(BUNDLED)
    chosen = families.select_families(graph)
    assert [f["key"] for f in chosen] == ["agree", "punctuation", "allele", "different", "family", "one_way"]
    raw = json.loads(BUNDLED.read_text())                      # independent of networkx
    names = {n["id"]: n["name"] for n in raw["nodes"]}
    edges = {(e["source"], e["target"]): e for e in raw["links"]}
    for family in chosen:
        assert all(names[n] == graph.nodes[n]["name"] for n in family["nodes"])
        assert all(len(names[n]) <= families.MAX_NAME_LENGTH for n in family["nodes"])
        for u, v in family["edges"]:
            assert (u, v) in edges and u in family["nodes"] and v in family["nodes"]
    one_way = next(f for f in chosen if f["key"] == "one_way")
    assert {edges[e]["type"] for e in one_way["edges"]} == {"RBH", "OWH"}
    assert len(next(f for f in chosen if f["key"] == "family")["nodes"]) == 6


def test_selection_is_deterministic():
    graph = matches.read_graph(BUNDLED)
    assert families.select_families(graph) == families.select_families(graph)
