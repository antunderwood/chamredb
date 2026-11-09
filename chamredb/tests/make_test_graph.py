"""Build chamredb/tests/test_data/graph_2021_subset.json: the part of the 2021 graph that the unit tests query.

The unit tests used to read the bundled graph, so every data refresh would change their expected outputs. They now
read this fixed subset (every queried source node, its out-neighbours and the edges between them, with all node
attributes), so the expected TSV and text outputs stay valid and a refresh cannot change them.
Provenance only; the fixture is committed. Usage (from the repository root):
    git show 359d3fe:chamredb/data/graph.json > /tmp/graph_2021.json   # the 2021 graph as shipped
    python3 -c "import sys; sys.path.insert(0, 'chamredb/tests'); import make_test_graph as m; m.main('/tmp/graph_2021.json')"
"""
import json
import os
import sys

import networkx as nx
from networkx.readwrite import json_graph

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from chamredb.functions import graph_functions, utility_functions  # noqa: E402

DATA = os.path.join(HERE, "test_data")
OUT = os.path.join(DATA, "graph_2021_subset.json")


def queried_ids():
    ids = [{"id": "WP_012695489.1", "database": "ncbi"}, {"id": "ARO:3002735", "database": "card"}]  # second: CLI smoke test
    ids += [{"id": i, "database": "card"} for i in utility_functions.parse_id_file(os.path.join(DATA, "card_ids.txt"))]
    ids += utility_functions.parse_hamronization_json_file(os.path.join(DATA, "multiple_dbs_summary.json"))
    return ids


def main(graph_path):
    graph = json_graph.node_link_graph(json.load(open(graph_path)), edges="links")
    keep = set()
    for q in queried_ids():
        node_id, _ = graph_functions.get_graph_node(q["id"], q["database"], graph)
        if node_id:
            keep.add(node_id)
            keep.update(graph.successors(node_id))
    sub = graph.subgraph(keep).copy()
    with open(OUT, "w") as fh:
        fh.write(json.dumps(json_graph.node_link_data(sub, edges="links"), sort_keys=True, indent=2))
    print(f"{sub.number_of_nodes()} nodes, {sub.number_of_edges()} edges")


if __name__ == "__main__":
    main(sys.argv[1])
