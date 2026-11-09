"""Graph loading from a pinned release, a cache and --graph PATH (ticket #6, ADR 0004).

No network: the 'release asset' is a file:// URL to a tiny graph made here (level 1, by construction).
"""
import hashlib
import json
import os
import subprocess
import sys
import tempfile

import pytest

from chamredb import graph_loader
from chamredb.functions import graph_functions, graph_visualisation_functions, utility_functions

TINY = {
    "directed": True, "multigraph": False, "graph": {},
    "nodes": [{"id": "card:ARO:1", "database": "card", "name": "a"}, {"id": "ncbi:WP_1", "database": "ncbi", "name": "b"}],
    "links": [{"source": "card:ARO:1", "target": "ncbi:WP_1", "type": "RBH", "identity": 1.0, "coverage": 1.0},
              {"source": "ncbi:WP_1", "target": "card:ARO:1", "type": "RBH", "identity": 1.0, "coverage": 1.0}],
}


def make_release(tmp_path, graph=TINY, data_release="test-release"):
    asset = tmp_path / "asset" / "graph.json"
    asset.parent.mkdir()
    asset.write_text(json.dumps(graph))
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps({"data_release": data_release, "url": asset.as_uri(),
                               "sha256": hashlib.sha256(asset.read_bytes()).hexdigest()}))
    return asset, pin


def test_read_graph_from_explicit_path(tmp_path):
    asset, _ = make_release(tmp_path)
    assert graph_functions.read_graph(asset).number_of_nodes() == 2


def test_download_verifies_caches_and_then_works_offline(tmp_path):
    asset, pin = make_release(tmp_path)
    cache = tmp_path / "cache"
    path, release = graph_loader.resolve_graph(pin_file=pin, cache_dir=cache)
    assert release == "test-release" and path.read_text() == asset.read_text() and cache in path.parents
    asset.unlink()  # source gone: the cached copy is used, no download needed
    path2, _ = graph_loader.resolve_graph(pin_file=pin, cache_dir=cache)
    assert path2 == path


def test_download_with_wrong_sha256_is_rejected_and_not_cached(tmp_path):
    asset, pin = make_release(tmp_path)
    p = json.loads(pin.read_text())
    p["sha256"] = "0" * 64
    pin.write_text(json.dumps(p))
    cache = tmp_path / "cache"
    with pytest.raises(ValueError, match="sha256"):
        graph_loader.resolve_graph(pin_file=pin, cache_dir=cache)
    assert not list(cache.rglob("graph.json"))


def test_corrupted_cache_is_not_trusted(tmp_path):
    asset, pin = make_release(tmp_path)
    cache = tmp_path / "cache"
    path, _ = graph_loader.resolve_graph(pin_file=pin, cache_dir=cache)
    path.write_text("corrupted")
    path2, _ = graph_loader.resolve_graph(pin_file=pin, cache_dir=cache)  # re-downloaded
    assert path2.read_text() == asset.read_text()


def test_default_is_the_bundled_graph_and_its_pin_matches():
    path, release = graph_loader.resolve_graph()
    pin = json.loads((graph_loader.PACKAGE_DATA / "graph_release.json").read_text())
    assert release == pin["data_release"] and pin["url"] is None
    assert hashlib.sha256(path.read_bytes()).hexdigest() == pin["sha256"]  # editing graph.json requires a new pin


def run_cli(*args, **env):
    e = {**os.environ, "COLUMNS": "200", "PYTHONIOENCODING": "utf-8", **env}
    return subprocess.run([sys.executable, "-m", "chamredb.run_chamredb", "query", *args],
                          capture_output=True, text=True, env=e, check=True).stdout


def test_cli_graph_flag_uses_the_given_graph_and_reports_its_release(tmp_path):
    asset, _ = make_release(tmp_path)
    out = run_cli("--graph", str(asset), "-d", "card", "-i", "ARO:1")
    assert "WP_1" in out and "data release: unknown (custom graph)" in out


def test_cli_prints_data_release_of_the_default_graph():
    pin = json.loads((graph_loader.PACKAGE_DATA / "graph_release.json").read_text())
    out = run_cli("-d", "card", "-i", "ARO:3002735")
    assert f"data release: {pin['data_release']}" in out and "QnrB20" in out


def test_tsv_gets_data_release_as_extra_last_column_and_otherwise_unchanged():
    graph = graph_functions.read_graph(os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data", "graph_2021_subset.json"))
    id_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data", "card_ids.txt")
    id_data = [{"id": i, "database": "card"} for i in utility_functions.parse_id_file(id_file)]
    with tempfile.NamedTemporaryFile(mode="r") as plain, tempfile.NamedTemporaryFile(mode="r") as tagged:
        graph_visualisation_functions.write_multiple_node_info(id_data, graph, plain.name)
        graph_visualisation_functions.write_multiple_node_info(id_data, graph, tagged.name, data_release="2021-10-20")
        a, b = plain.read().splitlines(), tagged.read().splitlines()
    assert len(a) == len(b) and a[0] + "\tdata_release" == b[0]
    assert all(y == x + "\t2021-10-20" for x, y in zip(a[1:], b[1:]))
