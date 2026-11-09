import pytest
import os
import json
import networkx as nx
from networkx.readwrite import json_graph
import pandas as pd
from chamredb.functions import graph_functions

@pytest.fixture(scope="module") 
def expected_rbh_graph_json():
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data", "card_vs_ncbi.rbh.graph.json")) as json_file:
        graph_json = json.loads(json_file.read())
    return graph_json

@pytest.fixture(scope="module") 
def expected_search_graph_json():
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data", "card_vs_ncbi.search.graph.json")) as json_file:
        graph_json = json.loads(json_file.read())
    return graph_json

def test_adding_rbh_hits(expected_rbh_graph_json):
    G = nx.DiGraph()
    rbh_data = pd.read_csv(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data", "card_vs_ncbi.rbh.sorted.tsv"),
        sep="\t"
    )
    graph_functions.add_rbh_hits_to_graph("card","ncbi",rbh_data,G)
    graph_json = json.loads(json.dumps(json_graph.node_link_data(G, edges="links"), sort_keys=True))
    assert expected_rbh_graph_json == graph_json
    
def test_adding_search_hits(expected_search_graph_json):
    G = nx.DiGraph()
    search_data = pd.read_csv(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data", "card_vs_ncbi.search.sorted.tsv"),
        sep="\t"
    )
    graph_functions.add_search_hits_to_graph("card","ncbi",search_data,G)
    graph_json = json.loads(json.dumps(json_graph.node_link_data(G, edges="links"), sort_keys=True))
    assert expected_search_graph_json == graph_json
