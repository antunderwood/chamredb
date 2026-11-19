import argparse
import os
import sys
import textwrap
from chamredb.query_graph import query_graph_single_id, query_graph_multiple_ids
from chamredb.functions.utility_functions import parse_id_file
from chamredb.functions.utility_functions import parse_hamronization_json_file

def is_valid_file(parser, arg):
    """A function to check that a filepath param passed as a string exists

    Args:
        parser (argparse): An argparse parser object
        arg (string): The file path to test

    Returns:
        string or parser.error: The filepath if it exists else a parser.error
    """
    if not os.path.isfile(arg):
        parser.error('The file {} does not exist!'.format(arg))
    else:
        # File exists so return the filename
        return arg


def is_valid_dir(parser, arg):
    """A function to check that a directory param passed as a filepath exists


    Args:
        parser (argparse): An argparse parser oblject
        arg (string): The path of the directory to test

    Returns:
        string or parser.error: The directory filepath if it exists else a parser.error
    """
    if not os.path.isdir(arg):
        parser.error('The directory {} does not exist!'.format(arg))
    else:
        # File exists so return the directory
        return arg

def parse_arguments():
    description = """
    A package to consolidate information about ARGs across multiple ARG databases
    """
    # parse all arguments
    parser = argparse.ArgumentParser(description=description,formatter_class=argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(
        help='The following commands are available. Type chamredb <COMMAND> -h for more help on a specific commands',
        dest='command'
    )
    
    subparsers.required = True

    # query command
    query = subparsers.add_parser('query', help='Query databases for matches and associated metadata')
    query.add_argument('-d', '--database', help='which database are the gene(s) in', choices=['card', 'ncbi', 'resfinder'], required='-i' in sys.argv or '--id' in sys.argv or '-f' in sys.argv or '--id_file' in sys.argv)
    query.add_argument('-ct', '--coverage_threshold', help='coverage threshold below which a match will not be reported', type=float, default=0.9)
    query.add_argument('-it', '--identity_threshold', help='identity threshold below which a match will not be reported', type=float, default=0.9)
    query.add_argument('--graph', help='Path to a graph JSON file to use instead of the default pinned release (offline use)', default=None)
    id_source = query.add_mutually_exclusive_group(required=True)
    id_source.add_argument('-i', '--id',  help='The id of a ARG in the specified database')
    id_source.add_argument('-f', '--id_file', help='Path to a file containing ids of ARGs in the specified database', type=lambda x: is_valid_file(parser, x))
    id_source.add_argument('-j', '--hamronization_json_file', help='Path to a hamronization summaty in JSON format', type=lambda x: is_valid_file(parser, x))
    query.add_argument('-o', '--outfile_path', help='Path to file where query results will be written', required='-f' in sys.argv or '--id_file' in sys.argv or '-j' in sys.argv or '--hamronization_json_file' in sys.argv)

    variants = subparsers.add_parser(
        'variants',
        help='Compare resistance mutations (SNPs, indels, promoter and rRNA changes) between CARD, NCBI and PointFinder',
        description='Compare resistance mutations between CARD, NCBI and PointFinder. Positions are reported in the '
                    'numbering of each locus\'s anchor database (CARD, else NCBI, else PointFinder).')
    variants.add_argument('-d', '--database', help='database the mutation(s) come from', choices=['card', 'ncbi', 'pointfinder'])
    variants.add_argument('-t', '--taxon', help='organism or organism group, e.g. Escherichia, "Mycobacterium tuberculosis"')
    variants.add_argument('--variants', help='Path to a variants JSON file to use instead of the bundled one (offline use)', default=None)
    v_source = variants.add_mutually_exclusive_group(required=True)
    v_source.add_argument('-i', '--id', help='a mutation: gene and change in any common form (gyrA_S83L, "gyrA p.Ser83Leu", "23S A2059G"), a record id, or several joined by commas for a combination')
    v_source.add_argument('-l', '--locus', help='show all mutations of a locus side by side, e.g. Escherichia:gyrA')
    v_source.add_argument('-f', '--id_file', help='file with one mutation per line (with -d and optionally -t)', type=lambda x: is_valid_file(parser, x))
    v_source.add_argument('-j', '--hamronization_json_file', help='hAMRonization summary JSON; mutation rows of AMRFinderPlus, RGI and PointFinder are used', type=lambda x: is_valid_file(parser, x))
    variants.add_argument('-o', '--outfile_path', help='TSV output for -f and -j', required='-f' in sys.argv or '--id_file' in sys.argv or '-j' in sys.argv or '--hamronization_json_file' in sys.argv)
    return(parser.parse_args())


def print_ascii_header():
    header = textwrap.dedent(
    """
    ====================================================
          _              __  __ _____      _____  _     
         | |       /\   |  \/  |  __ \    |  __ \| |    
      ___| |__    /  \  | \  / | |__) |___| |  | | |__  
     / __| '_ \  / /\ \ | |\/| |  _  // _ \ |  | | '_ \ 
    | (__| | | |/ ____ \| |  | | | \ \  __/ |__| | |_) |
     \___|_| |_/_/    \_\_|  |_|_|  \_\___|_____/|_.__/ 
    ====================================================
    """
    )
    print(header)

def main():
    """The main function to parse the arguments and then direct the command as specified by the params
    """
    options = parse_arguments()
    print_ascii_header()
    if options.command == 'query':
        if options.id:
            query_graph_single_id(options.id, options.database, options.coverage_threshold, options.identity_threshold, options.graph)
        elif options.id_file:
            ids = parse_id_file(options.id_file)
            id_data = [{'id': id, 'database': options.database} for id in ids]
            query_graph_multiple_ids(id_data, options.outfile_path, options.coverage_threshold, options.identity_threshold, options.graph)
        elif options.hamronization_json_file:
            id_data = parse_hamronization_json_file(options.hamronization_json_file)
            query_graph_multiple_ids(id_data, options.outfile_path, options.coverage_threshold, options.identity_threshold, options.graph)
            print(f'\nResults available at {options.outfile_path}')
    elif options.command == 'variants':
        from chamredb import query_variants
        if options.id:
            sys.exit(query_variants.query_single(options.id, options.database, options.taxon, options.variants))
        elif options.locus:
            sys.exit(query_variants.query_locus(options.locus, options.variants))
        elif options.id_file:
            items = [{'id': line, 'database': options.database, 'taxon': options.taxon} for line in parse_id_file(options.id_file) if line.strip()]
            query_variants.query_multiple(items, options.outfile_path, options.variants)
            print(f'\nResults available at {options.outfile_path}')
        elif options.hamronization_json_file:
            items = query_variants.parse_hamronization_mutations(options.hamronization_json_file)
            for item in items:
                item['taxon'] = options.taxon
            query_variants.query_multiple(items, options.outfile_path, options.variants)
            print(f'\n{len(items)} mutation rows; results available at {options.outfile_path}')

if __name__ == "__main__":
    main()
