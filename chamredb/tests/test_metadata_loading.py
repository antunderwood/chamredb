"""Database metadata that goes with the graph (ticket #33).

The metadata JSON files name ids that are missing from the graph in the 'could not find matches' message. They
must come from the same release as the graph: pinned by sha256 next to it, and every graph node must be in them.
"""
import hashlib
import json
import os

import pytest

from chamredb import graph_loader
from chamredb.functions import graph_functions, graph_visualisation_functions

GRAPH = {
    "directed": True, "multigraph": False, "graph": {},
    "nodes": [{"id": "card:ARO:1", "database": "card", "name": "a"}, {"id": "ncbi:WP_1", "database": "ncbi", "name": "b"}],
    "links": [{"source": "card:ARO:1", "target": "ncbi:WP_1", "type": "RBH", "identity": 1.0, "coverage": 1.0},
              {"source": "ncbi:WP_1", "target": "card:ARO:1", "type": "RBH", "identity": 1.0, "coverage": 1.0}],
}


def make_pin(tmp_path, card_meta):
    meta_dir = tmp_path / "db_metadata"
    meta_dir.mkdir()
    (meta_dir / "card.metadata.json").write_text(json.dumps(card_meta))
    sha = hashlib.sha256((meta_dir / "card.metadata.json").read_bytes()).hexdigest()
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps({"data_release": "t", "sha256": "0" * 64, "url": None,
                               "metadata": {"card": {"file": "db_metadata/card.metadata.json", "sha256": sha}}}))
    return pin


def test_metadata_is_read_from_the_pinned_file(tmp_path):
    pin = make_pin(tmp_path, {"ARO:9": {"name": "newgene"}})
    assert graph_loader.resolve_metadata("card", pin_file=pin) == {"ARO:9": {"name": "newgene"}}


def test_metadata_with_the_wrong_sha256_is_rejected(tmp_path):
    pin = make_pin(tmp_path, {"ARO:9": {"name": "newgene"}})
    (tmp_path / "db_metadata" / "card.metadata.json").write_text(json.dumps({"ARO:9": {"name": "changed"}}))
    with pytest.raises(ValueError, match="sha256"):
        graph_loader.resolve_metadata("card", pin_file=pin)


def test_a_custom_graph_has_no_pinned_metadata():
    assert graph_loader.resolve_metadata("card", graph_path="whatever.json") is None


def test_missing_ids_are_named_from_the_metadata_of_the_same_release(tmp_path, capsys):
    id_data = [{"id": "ARO:1", "database": "card"}, {"id": "ARO:9", "database": "card"}]
    from networkx.readwrite import json_graph
    graph = json_graph.node_link_graph(GRAPH, edges="links")
    graph_visualisation_functions.write_multiple_node_info(
        id_data, graph, str(tmp_path / "o.tsv"), metadata=lambda db: {"ARO:9": {"name": "newgene"}})
    assert "ARO:9 (newgene) in 1" in capsys.readouterr().out


def test_every_node_of_the_bundled_graph_is_in_the_bundled_metadata():
    graph = graph_functions.read_graph()
    missing = {}
    for node_id in graph.nodes:
        database, identifier = node_id.split(":", 1)
        meta = graph_loader.resolve_metadata(database)
        if identifier not in meta:
            missing.setdefault(database, []).append(identifier)
    assert not missing, {db: (len(v), v[:3]) for db, v in missing.items()}
