from functions import hit_parsing
from functions import graph_building

import os
import json
import networkx as nx
from networkx.readwrite import json_graph
import itertools

G = nx.DiGraph()

print("Building graph...")

for permutation in itertools.permutations(['card', 'ncbi', 'resfinder'], r=2):
    sourceDB = permutation[0]
    targetDB = permutation[1]
    print("    - {} vs {}".format(sourceDB, targetDB))
    hit_parsing.filter_and_sort_rbhs(sourceDB,targetDB)
    hit_parsing.filter_and_sort_non_rbhs(sourceDB,targetDB)
    rbh_data = graph_building.get_rbh_data(sourceDB,targetDB)
    graph_building.add_rbh_hits_to_graph(sourceDB,targetDB,rbh_data,G)
    search_data = graph_building.get_search_data(sourceDB,targetDB)
    graph_building.add_search_hits_to_graph(sourceDB,targetDB,search_data,G)


out_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..", "graph.json"
            )
print("Saving graph to {}".format(out_path))
with open(out_path, "w") as out_file:
    out_file.write(
        json.dumps(json_graph.node_link_data(G), sort_keys=True, indent=2)
    )


