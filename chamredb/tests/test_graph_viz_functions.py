import pytest
import os
import textwrap
import tempfile
from chamredb.functions import graph_visualisation_functions
from chamredb.functions import graph_functions
from chamredb.functions import utility_functions

@pytest.fixture(scope="module") 
def expected_print_node_info_output():
    return textwrap.dedent(
        """
        :dna: [cyan bold]WP_012695489.1[/cyan bold] [magenta bold](qnrB2)[/magenta bold]
        [bright_cyan]:page_facing_up: phenotype: confers resistance to subclass QUINOLONE[/bright_cyan]
        [bright_cyan]:page_facing_up: product: quinolone resistance pentapeptide repeat protein QnrB2[/bright_cyan]
          :file_cabinet: [chartreuse1]card[/chartreuse1]
            :left_right_arrow: [chartreuse1]ARO3002735 (QnrB20)[/chartreuse1]
                :link: [white]coverage:[/white] [grey66]1.0[/grey66]
                :link: [white]identity:[/white] [grey66]0.995[/grey66]
                :link: [white]type:[/white] [grey66]RBH[/grey66]
                :page_facing_up: [white]PMID:[/white] [grey66]['18993034'][/grey66]
                :page_facing_up: [white]additional_phenotype:[/white] [grey66]confers resistance to drug class: fluoroquinolone antibiotic[/grey66]
                :page_facing_up: [white]is_a:[/white] [grey66]ARO:3000419:quinolone resistance protein (qnr)[/grey66]
                :page_facing_up: [white]phenotype:[/white] [grey66]confers resistance to antibiotic: ciprofloxacin,gatifloxacin,levofloxacin,moxifloxacin,nalidixic acid,norfloxacin,sparfloxacin[/grey66]
          :file_cabinet: [orange_red1]resfinder[/orange_red1]
            :left_right_arrow: [orange_red1]qnrB2 (qnrB2)[/orange_red1]
                :link: [white]coverage:[/white] [grey66]1.0[/grey66]
                :link: [white]identity:[/white] [grey66]1.0[/grey66]
                :link: [white]type:[/white] [grey66]RBH[/grey66]
                :page_facing_up: [white]phenotype:[/white] [grey66]confers resistance to Ciprofloxacin[/grey66]
        ================================================================================

        """
    )

@pytest.fixture(scope="module") 
def chamredb_graph():
    # fixed subset of the 2021 graph (chamredb/tests/make_test_graph.py): data refreshes cannot change these expected outputs
    return graph_functions.read_graph(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data', 'graph_2021_subset.json'))


def test_single_node_info_text(expected_print_node_info_output, chamredb_graph):
    assert expected_print_node_info_output == graph_visualisation_functions.single_node_info_text(
        'WP_012695489.1',
        'ncbi',
        chamredb_graph
    )

def test_write_multiple_node_info_from_ids(chamredb_graph):
    # get id_data from hamronization summary file
    id_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_data', 'card_ids.txt')
    ids= utility_functions.parse_id_file(id_file)
    id_data = [{'id': id, 'database': 'card'} for id in ids]
    # get expected output filepath
    expected_multiple_ids_output_file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_data', 'expected_multiple_ids.info.tsv')
    with open(expected_multiple_ids_output_file_path) as expected_multiple_ids_output_file:
        expected_multiple_dbs_output = expected_multiple_ids_output_file.read()
    with tempfile.NamedTemporaryFile(mode = "w") as temp_output:
        # run chamredb over id_data
        graph_visualisation_functions.write_multiple_node_info(
            id_data,
            chamredb_graph,
            temp_output.name
        )
        with open(temp_output.name) as multiple_ids_output_file:
            print(multiple_ids_output_file)
            multiple_ids_output = multiple_ids_output_file.read()
            assert expected_multiple_dbs_output == multiple_ids_output

def test_write_multiple_node_info_from_hamronization(chamredb_graph):
    # get id_data from hamronization summary file
    json_summary_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_data', 'multiple_dbs_summary.json')
    id_data = utility_functions.parse_hamronization_json_file(json_summary_file)
    # get expected output filepath
    expected_multiple_dbs_output_file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_data', 'expected_multiple_dbs.info.tsv')
    with open(expected_multiple_dbs_output_file_path) as expected_multiple_dbs_output_file:
        expected_multiple_dbs_output = expected_multiple_dbs_output_file.read()
    with tempfile.NamedTemporaryFile(mode = "w") as temp_output:
        # run chamredb over id_data
        graph_visualisation_functions.write_multiple_node_info(
            id_data,
            chamredb_graph,
            temp_output.name
        )
        with open(temp_output.name) as multiple_dbs_output_file:
            multiple_dbs_output = multiple_dbs_output_file.read()
            print(expected_multiple_dbs_output)
            print("============================================================")
            print(multiple_dbs_output)
            assert expected_multiple_dbs_output == multiple_dbs_output



