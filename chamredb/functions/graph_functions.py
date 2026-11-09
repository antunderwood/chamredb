import pandas as pd
from chamredb import graph_loader
from networkx.readwrite import json_graph
import os
import json
import itertools

def add_rbh_hits_to_graph(sourceDB,targetDB,rbh_data,graph,debug=False,metadata=None):
    # metadata: optional (source metadata dict, target metadata dict); default reads the bundled files
    sourceDB_metadata,targetDB_metadata = metadata or __load_metadata(sourceDB,targetDB)
    for _, row in rbh_data.iterrows():
        source_node_attributes = sourceDB_metadata[row['query']]
        source_node_attributes['database'] = sourceDB
        target_node_attributes = targetDB_metadata[row['target']]
        target_node_attributes['database'] = targetDB
        __add_connection(
            graph,
            row['query'],
            source_node_attributes,
            row['target'],
            target_node_attributes,
            {
                'type': 'RBH',
                'identity': row['fident'],
                'coverage': row['alnlen']/row['qend']
            },
            debug = debug
        )

def add_search_hits_to_graph(sourceDB,targetDB,search_data,graph,debug=False,metadata=None):
    sourceDB_metadata,targetDB_metadata = metadata or __load_metadata(sourceDB,targetDB)

    # OWH == One Way Hit
    for _, row in search_data.iterrows():
        source_node_attributes = sourceDB_metadata[row['query']]
        source_node_attributes['database'] = sourceDB
        target_node_attributes = targetDB_metadata[row['target']]
        target_node_attributes['database'] = targetDB
        __add_connection(
            graph,
            row['query'],
            source_node_attributes,
            row['target'],
            target_node_attributes,
            {
                'type': 'OWH',
                'identity': row['fident'],
                'coverage': row['alnlen']/row['qend']
            },
            debug=debug
        )

def get_graph_node(identifier,database,graph):
    """
    gets node from graph first trying the id in the graph which is 'database:id', else searching name, else returning None
    """
    full_id = f'{database}:{identifier}'
    if full_id in graph.nodes:
        return(full_id, graph.nodes[full_id])
    else:
        for node_id, node_data in graph.nodes(data=True):
            if ((
                    'alternative_id' in node_data
                    and identifier in node_data['alternative_id'].values()
                ) or (
                    'duplicate_allele_ids' in node_data
                    and identifier in itertools.chain(*node_data['duplicate_allele_ids'].values())
                )
                or
                node_data['name'] == identifier
            ) and node_data['database'] == database:
                return(node_id, graph.nodes[node_id])
        return(None, None)

def read_graph(graph_path=None):
    """Load the graph from graph_path, or from the pinned default release (see chamredb.graph_loader)."""
    json_path, _ = graph_loader.resolve_graph(graph_path)
    with open(json_path) as json_file:
        G = json_graph.node_link_graph(json.load(json_file), edges="links")
    return(G)

# private methods
def __load_metadata(sourceDB,targetDB):
    sourceDB_metadata_path = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "..", "data", "db_metadata",
                    f"{sourceDB}.metadata.json"
                )
    with open(sourceDB_metadata_path) as sourceDB_metadata_file:
        sourceDB_metadata = json.load(sourceDB_metadata_file)

    targetDB_metadata_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..", "data", "db_metadata",
                f"{targetDB}.metadata.json"
            )
    with open(targetDB_metadata_path) as targetDB_metadata_file:
        targetDB_metadata = json.load(targetDB_metadata_file)
    
    return(sourceDB_metadata,targetDB_metadata)

def __add_connection(graph,
                    source_node_id,
                    source_node_attributes,
                    target_node_id,
                    target_node_attributes,
                    edge_attributes,
                    debug=False):
    source_composite_key = f'{source_node_attributes["database"]}:{source_node_id}'
    target_composite_key = f'{target_node_attributes["database"]}:{target_node_id}'
    if not graph.has_node(source_composite_key):
        graph.add_node(source_composite_key, **source_node_attributes)
    if not graph.has_node(target_composite_key):
        graph.add_node(target_composite_key, **target_node_attributes)
    if graph.has_edge(source_composite_key, target_composite_key):
        if debug:
            print(f"Warning! Graph already has edge {source_composite_key} => {target_composite_key}")
    else:
        graph.add_edge(source_composite_key, target_composite_key, **edge_attributes)
    # add reciprocal edge if RBH
    if edge_attributes['type'] == 'RBH':
        if graph.has_edge(target_composite_key, source_composite_key):
            if debug:
                print(f"Warning! Graph already has edge {target_composite_key} => {source_composite_key}")
        else:
            graph.add_edge(target_composite_key, source_composite_key, **edge_attributes)
