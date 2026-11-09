import sys
import json
import os
import itertools
from rich.progress import track
from chamredb import graph_loader
from chamredb.functions import graph_functions

def single_node_info_text(identifier, database, graph, coverage_threshold = 0.9, identity_threshold = 0.9):
    """
    single_node_info_text Given an id and its corresponding database find matches and write the match metadata to a file

    Args:
        identifier (string): The identifier in the database
        database (string): One of ncbi,card or resfinder
        graph (networkx DiGraph): A graph of the AMR database matches
        coverage_threshold (float, optional): The coverage value below which a match will not be reported. Defaults to 0.9.
        identity_threshold (float, optional): The identity value below which a match will not be reported. Defaults to 0.9.

    Returns:
        [string]: A string formatted for the printing using the rich library
    """
    colours = ["blue_violet", "yellow1", "orange_red1","chartreuse1" ]
    source_node_id, source_node_data = graph_functions.get_graph_node(identifier, database, graph)
    if not source_node_id:
        print(f"Could not find a match for {identifier} in the {database} database")
        sys.exit(1)

    info_text = f"\n:dna: [cyan bold]{identifier}[/cyan bold] [magenta bold]({source_node_data['name']})[/magenta bold]\n"

    source_node_metadata = __node_metadata(source_node_data)
    for key in source_node_metadata:
        info_text += f"[bright_cyan]:page_facing_up: {key}: {source_node_metadata[key]}[/bright_cyan]\n"
    
    node_targets = __node_target_info(graph,source_node_id)
    database_colours = {}
    for database in node_targets:
        database_colours[database] = colours.pop()

    # print edge (hit) data
    for database in node_targets:
        colour = database_colours[database]
        info_text += f"  :file_cabinet: [{colour}]{database}[/{colour}]\n"
        # return just first match
        target_ids = list(node_targets[database])
        target_id = target_ids[0]
        edge_data = node_targets[database][target_id]['edge_data']
        if edge_data['coverage'] < coverage_threshold or edge_data['identity'] < identity_threshold:
            info_text += f"    No significant matches with coverage >= {coverage_threshold} and identity >= {identity_threshold}\n"
            continue
        target_node_data = node_targets[database][target_id]['node']
        if edge_data['type'] == 'RBH':
            info_text += f"    :left_right_arrow: [{colour}]{target_id} ({target_node_data['name']})[/{colour}]\n"
        else:
            info_text += f"      :right_arrow: [{colour}]{target_id} ({target_node_data['name']})[/{colour}]\n"
        for key in edge_data:
            info_text += f"        :link: [white]{key}:[/white] [grey66]{edge_data[key]}[/grey66]\n"
        other_data = {key:target_node_data[key] for key in target_node_data if key != 'database' and key != 'name'}
        for key in other_data:
            info_text += f"        :page_facing_up: [white]{key}:[/white] [grey66]{other_data[key]}[/grey66]\n"
    info_text += f"{'='*80}\n\n"
    return info_text


def write_multiple_node_info(id_data,graph,out_filepath,coverage_threshold=0.9,identity_threshold=0.9,data_release=None,metadata=None):
    """
    write_multiple_node_info Given a list of dict objects containing ids and their corresponding databases find
    matches and write the match metadata to a file

    Args:
        id_data (list):  id_data: a list of dicts with keys id, database and optionally file
        graph (networkx DIGraph): A graph of the AMR database matches
        out_filepath (string): path to the outfile that will contain the id matches
        coverage_threshold (float, optional): The covergae value below which a match will not be reported. Defaults to 0.9.
        identity_threshold (float, optional): The identity value below which a match will not be reported. Defaults to 0.9.
        data_release (string, optional): if given, written as an extra last column `data_release`
        metadata (callable, optional): database name -> metadata dict (or None), used to name ids missing from the graph;
            defaults to the metadata pinned with the bundled graph
    """
    suffix = f'\t{data_release}' if data_release else ''
    all_databases = set([graph.nodes[node_name]['database'] for node_name in graph.nodes])
    id_databases = set([ id_info['database'] for id_info in id_data])
    # if only one database in ids_and_databases then the header databases should not include this database
    if len(id_databases) == 1:
        header_databases = sorted(all_databases - id_databases)
    else:
        header_databases = sorted(all_databases)
    
    target_node_field_titles = ['id', 'name', 'match_type', 'match_identity', 'match_coverage', 'metadata'] 
    
    if 'file' in id_data[0]:
        multiple_samples = True
    else:
        multiple_samples = False
    
    with open(out_filepath, 'w') as out:
        header = __multiple_ids_header(header_databases, target_node_field_titles, multiple_samples)
        out.write(header.rstrip('\n') + ('\tdata_release' if data_release else '') + '\n')
        missing_ids = {}
        for id_info in track(id_data, description=f'Finding metadata'):
            source_id = id_info['id']
            source_database = id_info['database']
            # get a source node and it's id based on the database and an id (which could be a name)
            source_node_id, source_node = graph_functions.get_graph_node(source_id, source_database, graph)
            if not source_node_id:
                if source_database not in missing_ids:
                    missing_ids[source_database] = {}
                if source_id not in missing_ids[source_database]:
                    missing_ids[source_database][source_id] = 0
                missing_ids[source_database][source_id] += 1
                continue
            # get the name
            name = source_node['name']
            # make the metadata from the source node
            metadata_string = ', '.join([f'{key}: {value}' for key, value in __node_metadata(source_node).items()])
            # find all target nodes for the source node
            source_node_target_info = __node_target_info(graph, source_node_id)
            # make a dict of items about the target nodes
            target_node_info = {}
            for target_database in source_node_target_info:
                best_match = __best_target_node_match(source_node_target_info[target_database],coverage_threshold,identity_threshold)
                if best_match:
                    target_node_info[target_database] = {
                        'id': best_match['id'],
                        'name': best_match['node']['name'],
                        'match_type': best_match['edge_data']['type'],
                        'match_identity': best_match['edge_data']['identity'],
                        'match_coverage': best_match['edge_data']['coverage'],
                        'metadata': ', '.join([f'{key}: {value}' for key, value in __node_metadata(best_match['node']).items()])
                    }
                else:
                    target_node_info[target_database] = {
                        'id': '',
                        'name': '',
                        'match_type':'',
                        'match_identity': '',
                        'match_coverage': '',
                        'metadata': ''
                    }
            
            # write out target info
            target_node_info_list = []
            for target_database in header_databases:
                for field in target_node_field_titles:
                    if target_database == source_database or target_database not in target_node_info:
                        target_node_info_list.append('-')
                    else:
                        target_node_info_list.append(str(target_node_info[target_database][field]))
            target_node_info_string = "\t".join(target_node_info_list)
            if multiple_samples:
                source_file = id_info['file']
                out.write(f'{source_file}\t{source_id}\t{source_database}\t{name}\t{metadata_string}\t{target_node_info_string}{suffix}\n')
            else:
                out.write(f'{source_id}\t{source_database}\t{name}\t{metadata_string}\t{target_node_info_string}{suffix}\n')
        # print missing ids message
        __missing_ids_message(missing_ids, multiple_samples, id_data, metadata)


# private methods
def __node_metadata(node):
    """
    return metadata for a node as a dict
    """
    metadata = {}
    non_metadata_keys = ['name', 'database', 'alternative_id', 'duplicate_allele_ids']
    for key in node:
        if key not in non_metadata_keys :
            metadata[key] = node[key]
    return metadata

def __node_target_info(graph,node_id):
    """
    Returns dictionary containing information about the edges connected to a source node and the target nodes associated with the edge.
    The dict groups the edges by database of the target node and the database are the top level keys 
    """
    node_targets = {}
    # gather edge (hit) data
    for edge in graph.edges(node_id):
        database,id,target_node,edge_data = __edge_target_info(graph,edge)
        if database not in node_targets:
            node_targets[database] = {}
        node_targets[database][id] = {
            'node': target_node,
            'edge_data': edge_data
        }
    return node_targets

def __best_target_node_match(node_targets,coverage_threshold,identity_threshold):
    """
    Given a dict of target nodes (key is id) return best node based on identity and coverage
    """
    best_match_identity = 0
    best_match_coverage = 0
    best_match = None
    for target_node_id in node_targets:
        coverage = node_targets[target_node_id]['edge_data']['coverage']
        identity = node_targets[target_node_id]['edge_data']['identity']
        if  coverage >= coverage_threshold and coverage > best_match_coverage:
            if identity >= identity_threshold and identity > best_match_identity:
                best_match = node_targets[target_node_id]
                best_match['id'] = target_node_id
                best_match_coverage = coverage
                best_match_identity = identity
    return best_match

def __edge_target_info(graph,edge):
    target_node = graph.nodes[edge[1]]
    edge_data = graph.get_edge_data(edge[0], edge[1])
    database = edge[1].split(":")[0]
    id = ''.join(edge[1].split(":")[1:])
    return database,id,target_node,edge_data


def __multiple_ids_header(databases, field_titles, multiple_samples):
    database_header_titles = []
    for db in sorted(databases):
        for title in field_titles:
            database_header_titles.append(f'{db}: {title}')
    if multiple_samples:
        sample_header_string = "sample\tid\tdatabase\tname\tmetadata"
    else:
        sample_header_string = "id\tdatabase\tname\tmetadata"
    database_header_title_string = '\t'.join(database_header_titles)
    header = f"{sample_header_string}\t{database_header_title_string}\n"
    return header

def __missing_ids_message(missing_ids, multiple_samples, id_data, metadata_for=None):
    # read in database metadata (pinned with the graph; None when unknown, then ids are listed without names)
    metadata_for = metadata_for or (lambda database: graph_loader.resolve_metadata(database))
    metadata = {}
    databases = missing_ids.keys()
    for database in databases:
        metadata[database] = metadata_for(database)
    # print message about missing ids
    if multiple_samples:
        number_of_samples = len(set([ id_info['file'] for id_info in id_data]))
    if len(missing_ids) > 0:
        print(f"WARNING: Could not find a matches for some ids")
        for database in missing_ids:
            print(f"Missing Ids in {database}:")
            for missing_id in missing_ids[database]:
                num_missing_ids = missing_ids[database][missing_id]
                if multiple_samples:
                    proportion = f' ({round(num_missing_ids/number_of_samples*100,1)}%)'
                else:
                    proportion = ""
                database_name = __get_name_for_id(missing_id, metadata[database]) if metadata[database] else None
                if database_name:
                    database_name_string = f' ({database_name})'
                else:
                    database_name_string = ''
                print(f"\t{missing_id}{database_name_string} in {num_missing_ids}{proportion} samples")

def __get_name_for_id(id_to_match, metadata):
    """get name of AMR determinant from metadata based on id or alternative id

    Args:
        id_to_match (string): id or alternative id of AMR determinant to find in database
        metadata (json): The consistently encoded database metadata in JSON format
    
    Returns:
        name(string) or None: Name of AMR determinant or None if not found
    """
    for database_id in metadata:
        if id_to_match == database_id:
            return metadata[database_id]['name']
        elif(
                'alternative_id' in metadata[database_id]
                    and id_to_match in metadata[database_id]['alternative_id'].values()
            ) or (
                    'duplicate_allele_ids' in metadata[database_id]
                    and id_to_match in itertools.chain(*metadata[database_id]['duplicate_allele_ids'].values())
            ):
            return metadata[database_id]['name']
    return None