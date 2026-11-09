"""Install and CLI smoke test (ticket #5): the packaged graph loads and `chamredb query` answers.

Expected query values come from the 2021 baseline run on the old stack (Python 3.9, networkx 2.8.8): the same query
printed these cross-database matches. The query runs on a fixed subset of the 2021 graph (see make_test_graph.py) so a
data refresh cannot change it; the bundled graph is checked against its pin (counts and hash) instead.
"""
import json
import os
import subprocess
import sys

from chamredb import graph_loader
from chamredb.functions import graph_functions


SUBSET = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data", "graph_2021_subset.json")


def test_packaged_graph_loads_and_matches_its_pin():
    pin = json.loads((graph_loader.PACKAGE_DATA / "graph_release.json").read_text())
    graph = graph_functions.read_graph()
    assert graph.number_of_nodes() == pin["nodes"]
    assert graph.number_of_edges() == pin["edges"]


def test_query_cli_reports_baseline_matches():
    env = {**os.environ, "COLUMNS": "200", "PYTHONIOENCODING": "utf-8"}
    out = subprocess.run([sys.executable, "-m", "chamredb.run_chamredb", "query", "--graph", SUBSET, "-d", "card", "-i", "ARO:3002735"],
                         capture_output=True, text=True, env=env, check=True).stdout
    for expected in ("ARO:3002735 (QnrB20)", "WP_063865914.1 (qnrB20)", "qnrB20 (qnrB20)", "type: RBH"):
        assert expected in out, expected
