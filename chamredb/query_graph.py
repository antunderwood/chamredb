from chamredb.functions import graph_visualisation_functions
from chamredb.functions import graph_functions
from chamredb import graph_loader
from rich import print as rprint

def query_graph_single_id(id,database,coverage_threshold,identity_threshold,graph_path=None):
    """
    query_graph_single_id query the graph with AMR metadata with a single id from one database

    Args:
        id (string): the identifier used in the database
        database (string): database name
        coverage_threshold (float): coverage threshold below which a match will not be reported
        identity_threshold (float): identity threshold below which a match will not be reported
        graph_path (string, optional): graph file to use instead of the pinned default release

    """
    path, data_release = graph_loader.resolve_graph(graph_path)
    graph = graph_functions.read_graph(path)
    rprint(f'data release: {data_release}')
    single_id_info_text = graph_visualisation_functions.single_node_info_text(id,database,graph,coverage_threshold,identity_threshold)
    rprint(single_id_info_text)

def query_graph_multiple_ids(id_data,outfile_path,coverage_threshold,identity_threshold,graph_path=None):
    """
    query_graph_multiple_ids query the graph with AMR metdata with a list of identifiers from one or more databases

    Args:
        id_data (list):  id_data: a list of dicts with keys id, database and optionally file
        out_filepath (string): path to the outfile that will contain the id matches
        coverage_threshold (float): The covergae value below which a match will not be reported
        identity_threshold (float): The identity value below which a match will not be reported
        graph_path (string, optional): graph file to use instead of the pinned default release
    """
    path, data_release = graph_loader.resolve_graph(graph_path)
    graph = graph_functions.read_graph(path)
    graph_visualisation_functions.write_multiple_node_info(id_data,graph,outfile_path,coverage_threshold,identity_threshold,data_release=data_release,metadata=lambda database: graph_loader.resolve_metadata(database, graph_path))
        