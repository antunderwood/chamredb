
import pandas as pd
import os

def filter_and_sort_rbhs(sourceDB,targetDB):
    unsorted_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "data", "mmseqs2_results",
        f"{sourceDB}_vs_{targetDB}.rbh.tsv"
    )
    sorted_file = os.path.join(os.path.dirname(unsorted_file), f'{os.path.splitext(os.path.basename(unsorted_file))[0]}.sorted.tsv')

    headers = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits']
    unsorted_data = pd.read_csv(unsorted_file, sep="\t", names = headers)

    # sort by query and get best match using fident and alnlen
    sorted_data = (unsorted_data
        .sort_values(
            by=['query', 'fident', 'alnlen'],
            ascending=[True, False, False]
        )
    )
    # remove duplicates
    sorted_data = sorted_data.drop_duplicates(subset=['query'])
    # now resort based on query
    sorted_data = (sorted_data
        .sort_values(
            by=['query'],
            key = lambda x: x.str.lower()
        )
    )

    sorted_data.to_csv(sorted_file, index=False, sep="\t")

def filter_and_sort_non_rbhs(sourceDB,targetDB):
    # get rbh ids
    rbh_ids = list(pd.read_csv(
                    os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "data", "mmseqs2_results",
                        f"{sourceDB}_vs_{targetDB}.rbh.sorted.tsv"
                    ), sep="\t"
                )['query'])


    unsorted_file = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..", "..", "data", "mmseqs2_results",
                f"{sourceDB}_vs_{targetDB}.search.tsv"
            )
    sorted_file = os.path.join(os.path.dirname(unsorted_file), f'{os.path.splitext(os.path.basename(unsorted_file))[0]}.sorted.tsv')

    # get matches from search that do not have rbhs and name them non-rbhs 
    headers = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits']

    non_rbh_results = (
        pd.read_csv(unsorted_file, sep="\t", names=headers)
        .query("query != @rbh_ids")
        .sort_values(
            by=['query', 'fident', 'alnlen'],
            ascending=[True, False, False]
        )
    )

    non_rbh_results.to_csv(
        sorted_file,
        sep="\t",
        index=False
    )
