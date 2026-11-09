"""Install and CLI smoke test (ticket #5): the packaged graph loads and `chamredb query` answers.

Expected values come from the 2021 baseline run on the old stack (Python 3.9, networkx 2.8.8): the
same query printed these cross-database matches. The graph is read from the installed package data.
"""
import os
import subprocess
import sys

from chamredb.functions import graph_functions


def test_packaged_graph_loads():
    graph = graph_functions.read_graph()
    assert graph.number_of_nodes() == 11718
    assert graph.number_of_edges() == 33632


def test_query_cli_reports_baseline_matches():
    env = {**os.environ, "COLUMNS": "200", "PYTHONIOENCODING": "utf-8"}
    out = subprocess.run([sys.executable, "-m", "chamredb.run_chamredb", "query", "-d", "card", "-i", "ARO:3002735"],
                         capture_output=True, text=True, env=env, check=True).stdout
    for expected in ("ARO:3002735 (QnrB20)", "WP_063865914.1 (qnrB20)", "qnrB20 (qnrB20)", "type: RBH"):
        assert expected in out, expected
