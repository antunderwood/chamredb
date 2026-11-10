# project CharmeDb
# Appendicies
## Data Download
### CARD
#### Metadata
ARO ontology
```
wget https://github.com/arpcard/aro/raw/master/src/ontology/aro.obo -O chamredb/data/db_metadata/card.metadata.obo
```
#### Sequence data
##### Protein
```
wget https://card.mcmaster.ca/latest/data -O card_db/data
cd card_db; tar -xvf ./data protein_fasta_protein_homolog_model.fasta
```

##### Nucleotide
```
tar -xvf data ./nucleotide_fasta_protein_homolog_model.fasta
```

### NCBI

#### Metadata
```
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt -O db_metadata/ncbi.metadata.tsv
```

#### Sequence data
##### Protein
```
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt -O db_fastas/ncbi.protein.raw.fasta
```
Nucleotide sequences
```
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS -O db_fastas/ncbi.nucl.raw.fasta
```

### Resfinder

#### Metadata
```
wget https://bitbucket.org/genomicepidemiology/resfinder_db/raw/master/phenotypes.txt -O db_metadata/resfinder.metadata.tsv
```

#### Sequences
```
git clone https://bitbucket.org/genomicepidemiology/resfinder_db
cat resfinder_db/*.fsa > db_fastas/resfinder.nucl.fasta
```

## Data Parsing

### CARD
Write out protein file with just the CARD ARO ids as fasta headers
[parse_card_proteins.py](/scripts/parse_card_proteins.py)
This produced the output file [card.protein.fasta](/chamredb/data/db_fastas/card.protein.fasta) with 2979 proteins

Process the [ARO obo data](/chamredb/data/db_metadata/card.metadata.obo) to write out `confers_resistance..` and `is_a` metadata to a JSON file where the ARO ids are keys
[parse_card_metadata.py](/scripts/parse_card_metadata.py)
This produced the [JSON output file](/chamredb/data/db_metadata/card.metadata.json)

### NCBI
Only include those in the core AMR/AMR (acquired) set and remove duplicates with the script. Use the `refseq_protein_accession` as the fasta header
[parse_ncbi_proteins.py](/scripts/parse_ncbi_proteins.py)

This produced the output file [ncbi.protein.fasta](/chamredb/data/db_fastas/ncbi.protein.fasta) with 4542 protein sequences

Metadata JSON file was created using `refseq_protein_accession` as keys and drug `subclass` + `product` as the metadata values
[parse_ncbi_metadata.py](/scripts/parse_ncbi_metadata.py)
This produced the [JSON output file](/chamredb/data/db_metadata/ncbi.metadata.json)

### Resfinder
Found best ORF from nucleotide [sequences](/chamredb/data/db_fastas/resfinder.nucl.fasta) and removed duplicated gene/alleles with the script
[resfinder_orfs_to_proteins.py](/scripts/resfinder_orfs_to_proteins.py)

This produced the output file [resfinder.protein.fasta](../charmed//chamredb/data/db_fastas/resfinder.protein.fasta) with 2543 protein sequences

Metadata json was created using the [resfinder.metadata.tsv](/chamredb/data/db_metadata/resfinder.metadata.tsv) as input so that the gene name was the key and `phenotype` was the metadata value
This produced the [JSON output file](/chamredb/data/db_metadata/resfinder.metadata.json)

## Analyse for Reciprocal Best Hits (RBHs)
Used the [MMseqs2](https://github.com/soedinglab/MMseqs2) software that allows very fast protein clustering - see this [publication](https://www.nature.com/articles/s41467-018-04964-5)
Performance characteristics in this [comparison publication](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07132-6)

Commands run were
* For each query database (DB1) vs target database (DB2)
  * reciprocal blast hit analysis
    ```
    mmseqs easy-rbh db_fastas/<DB1>.protein.fasta  db_fastas/<DB2>.protein.fasta  results/mmseqs_<DB1>_vs_<DB2>.rbh.tsv tmp -s 7.5
    ```
  * run a search for top 3 hits so that proteins that did not return a RBH have matches to be reported
    ```
    # make database for query DB1
    mmseqs createdb db_fastas/<DB1>.protein.fasta mmseqs_DBs/<DB1>.protein
    mmseqs createindex mmseqs_DBs/<DB1>.protein tmp
    # make database for target DB2
    mmseqs createdb db_fastas/<DB1>.protein.fasta mmseqs_DBs/<DB1>.protein
    mmseqs createindex mmseqs_DBs/<DB1>.protein tmp

    # search DB2 with DB1 to get top 3 hits
    mmseqs search mmseqs_DBs/<DB1>.protein mmseqs_DBs/<DB2>.protein mmseqs_search_DBs/<DB1>_vs_<DB2>_search tmp -s 7.5 --max-accept 3
    mmseqs convertalis  mmseqs_DBs/<DB1>.protein mmseqs_DBs/<DB2>.protein mmseqs_search_DBs/<DB1>_vs_<DB2>_search results/mmseqs_<DB1>_vs_<DB2>.search.tsv
    ```
These commands can be found in [find_rbhs_and_search_matches.sh](scripts/find_rbhs_and_search_matches.sh)

### Summary of RBH analysis

For each pair of databases: how many proteins of the first database are in the graph (a protein is in the graph if it has at least one match in another database), and how many of them have a reciprocal best hit (RBH) in the second. Generated from `docs/statistics.json` by `python -m analysis.readme`; the percentages are of the proteins in the graph, not of all proteins parsed from the database.

<!-- table:rbh -->
| Comparison | Proteins of the first database in the graph | With a reciprocal best hit |
| --- | --- | --- |
| CARD vs NCBI | 5,971 | 5,662 (94.8%) |
| CARD vs ResFinder | 5,971 | 2,458 (41.2%) |
| NCBI vs CARD | 7,299 | 5,682 (77.8%) |
| NCBI vs ResFinder | 7,299 | 2,509 (34.4%) |
| ResFinder vs CARD | 2,589 | 2,461 (95.1%) |
| ResFinder vs NCBI | 2,589 | 2,504 (96.7%) |
<!-- /table -->

### Summary of non-RBH searches

After the RBH search, a second MMseqs2 search (sensitivity 7.0, up to three accepted hits) finds one-way hits for proteins without a reciprocal best hit. Combining the two gives the coverage of one database by another.

<!-- table:search -->
| Comparison | With a one-way hit only | Covered by a reciprocal best hit or a one-way hit |
| --- | --- | --- |
| CARD vs NCBI | 307 | 5,969 (100.0%) |
| CARD vs ResFinder | 3,414 | 5,872 (98.3%) |
| NCBI vs CARD | 1,613 | 7,295 (99.9%) |
| NCBI vs ResFinder | 4,683 | 7,192 (98.5%) |
| ResFinder vs CARD | 124 | 2,585 (99.8%) |
| ResFinder vs NCBI | 84 | 2,588 (100.0%) |
<!-- /table -->

## Building a Networkx graph
Steps include
* Sort and filter RBHs to 
  * sort by query then percent_id match and alignment length
  * remove duplicates (this will keep best match)
  * see the `filter_and_sort_rbhs` function in [hit_functions.py](/chamredb/functions/hit_functions.py)
* Filter the non-RBH search matches to find those that do not already have a RBH match
  * get search matches from search that do not have a corresponding accession in the RBH hits
  * sort by query then percent_id match and alignment length
  * see the `filter_and_sort_non_rbhs` function in [hit_functions.py](/chamredb/functions/hit_functions.py)
* Make a directed graph (DiGraph) in networkx where
  * the nodes have attributes pulled from the metadata
  * the edges link nodes which have attributes
    * `type` either `RBH` or `OWH` (one way hit) for a search match
    * `coverage` alignment length/query length
    * `identity` fraction of identical residues in the hit (0 to 1)
  * See functions `add_rbh_hits_to_graph` and `add_search_hits_to_graph` in [graph_functions.py](/chamredb/functions/graph_functions.py)

These are all combined in the [build_graph_from_mmseqs_data.py](/chamredb/build_graph_from_mmseqs_data.py) script
