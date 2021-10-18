import pandas as pd
from networkx.readwrite import json_graph
import os
import json
import itertools


def get_rbh_data(sourceDB,targetDB):
    """
    get_rbh_data get pre-calculated sorted reciprocal best hit data

    Args:
        sourceDB (string): source database
        targetDB (string): target database

    Returns:
        pandas dataframe: dataframe of mmseqs2 data when RBHs were found by comparing the source and target databases
    """
    rbh_data = pd.read_csv(
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "data", "mmseqs2_results",
                    f"{sourceDB}_vs_{targetDB}.rbh.sorted.tsv"
                ), sep="\t")
    return rbh_data


def add_rbh_hits_to_graph(sourceDB,targetDB,rbh_data,graph,debug=False):
    sourceDB_metadata,targetDB_metadata = __load_metadata(sourceDB,targetDB)
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

def get_search_data(sourceDB,targetDB):
    """
    get_search_data get pre-calculated sorted search data

    Args:
        sourceDB (string): source database
        targetDB (string): target database

    Returns:
        pandas dataframe: dataframe of mmseqs2 data when matches were found found by comparing the source and target databases
    """
    search_data = pd.read_csv(
                    os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "data", "mmseqs2_results",
                        f"{sourceDB}_vs_{targetDB}.search.sorted.tsv"
                    ), sep="\t")
    return search_data


def add_search_hits_to_graph(sourceDB,targetDB,search_data,graph,debug=False):
    sourceDB_metadata,targetDB_metadata = __load_metadata(sourceDB,targetDB)

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

# private methods
def __load_metadata(sourceDB,targetDB):
    sourceDB_metadata_path = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "data", "db_metadata",
                    f"{sourceDB}.metadata.json"
                )
    with open(sourceDB_metadata_path) as sourceDB_metadata_file:
        sourceDB_metadata = json.load(sourceDB_metadata_file)

    targetDB_metadata_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..", "..", "data", "db_metadata",
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
