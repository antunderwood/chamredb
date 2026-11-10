"""Equivalence groups: proteins joined by reciprocal best hits (ticket #49, spec #34).

Level 1: a hand-built graph whose groups and database profiles are worked out by hand in this file. Level 4: the groups of
the bundled release graph checked against a second implementation (union-find) written independently of networkx.
One-way hits never join groups; a node without an RBH is a group of its own.
"""
from pathlib import Path

import pytest
from networkx.readwrite import json_graph

from analysis import groups, matches, run

ROOT = Path(__file__).resolve().parents[2]
BUNDLED = ROOT / "chamredb" / "data" / "graph.json"


def pair(a, b):
    return [{"source": a, "target": b, "type": "RBH", "identity": 1.0, "coverage": 1.0},
            {"source": b, "target": a, "type": "RBH", "identity": 1.0, "coverage": 1.0}]


# 24 nodes in 8 groups, worked out by hand:
#   G1 c1 - n1 - r1 (a chain: c1 and r1 are not joined directly)      -> card+ncbi+resfinder, size 3
#   G2 c2 - n2                                                         -> card+ncbi, size 2
#   G3 a path of 12: c10 n10 c11 n11 c12 n12 c13 n13 c14 n14 c15 n15   -> card+ncbi, size 12 (6 + 6)
#   G4 n20 alone; G5 c20 alone (one-way hit c20 -> n20 does not join them)
#   G6 c30 - r30 (card+resfinder); G7 r40 alone; G8 n50 - r50 (ncbi+resfinder)
PATH = [f"{db}:{db[0]}{i}" for i in range(10, 16) for db in ("card", "ncbi")]
NODES = ["card:c1", "ncbi:n1", "resfinder:r1", "card:c2", "ncbi:n2"] + PATH + [
    "ncbi:n20", "card:c20", "card:c30", "resfinder:r30", "resfinder:r40", "ncbi:n50", "resfinder:r50"]
LINKS = (pair("card:c1", "ncbi:n1") + pair("ncbi:n1", "resfinder:r1") + pair("card:c2", "ncbi:n2")
         + [l for a, b in zip(PATH, PATH[1:]) for l in pair(a, b)]
         + [{"source": "card:c20", "target": "ncbi:n20", "type": "OWH", "identity": 0.8, "coverage": 0.9}]
         + pair("card:c30", "resfinder:r30") + pair("ncbi:n50", "resfinder:r50"))
HAND = json_graph.node_link_graph(
    {"directed": True, "multigraph": False, "graph": {}, "nodes": [{"id": n} for n in NODES], "links": LINKS}, edges="links")


def test_groups_and_profiles_on_the_hand_graph():
    found = groups.equivalence_groups(HAND)
    assert len(NODES) == 24 and sum(len(g) for g in found) == 24 and len(found) == 8
    profile = groups.profiles(found)
    assert profile["groups"] == {"card": 1, "card+ncbi": 2, "card+ncbi+resfinder": 1, "card+resfinder": 1,
                                 "ncbi": 1, "ncbi+resfinder": 1, "resfinder": 1}
    assert profile["proteins"] == {"card": 1, "card+ncbi": 14, "card+ncbi+resfinder": 3, "card+resfinder": 2,
                                   "ncbi": 1, "ncbi+resfinder": 2, "resfinder": 1}
    assert profile["size_distribution"] == {1: 3, 2: 3, 3: 1, 12: 1}
    assert profile["largest_group_size"] == 12 and profile["largest_group_by_database"] == {"card": 6, "ncbi": 6}


def test_one_way_hits_do_not_join_groups_and_a_node_without_rbh_is_its_own_group():
    found = groups.equivalence_groups(HAND)
    by_node = {n: g for g in found for n in g}
    assert by_node["card:c20"] != by_node["ncbi:n20"]
    assert by_node["card:c20"] == frozenset({"card:c20"}) and by_node["ncbi:n20"] == frozenset({"ncbi:n20"})
    assert by_node["card:c1"] == by_node["resfinder:r1"]  # joined through n1


def union_find_groups(graph):
    """Independent implementation: union-find over the undirected RBH pairs."""
    parent = {n: n for n in graph.nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u, v, d in graph.edges(data=True):
        if d["type"] == "RBH":
            parent[find(u)] = find(v)
    out = {}
    for n in graph.nodes:
        out.setdefault(find(n), set()).add(n)
    return {frozenset(g) for g in out.values()}


def test_groups_equal_a_union_find_on_the_hand_graph_and_the_release_graph():
    assert set(groups.equivalence_groups(HAND)) == union_find_groups(HAND)
    release = matches.read_graph(BUNDLED)
    found = groups.equivalence_groups(release)
    assert set(found) == union_find_groups(release)
    assert sum(len(g) for g in found) == release.number_of_nodes()
    profile = groups.profiles(found)
    assert sum(profile["groups"].values()) == len(found) and sum(profile["proteins"].values()) == release.number_of_nodes()


def test_group_statistics_as_written_to_statistics_json():
    g = run.group_statistics(HAND)
    assert g["total"] == 8 and g["by_databases"]["card+ncbi"] == 2 and g["proteins_by_databases"]["card+ncbi"] == 14
    assert g["largest_group_size"] == 12 and g["largest_group_by_database"] == {"card": 6, "ncbi": 6}
    assert g["size_distribution"] == {"1": 3, "2": 3, "3": 1, "12": 1}   # JSON keys are strings
