# CharmeDb

![CharmeDb](docs/images/CharmeDb.png)  
(pronounced 'charmed' `/tʃɑː(r)md/`)
  
Previously known as Project ![mAMRite](docs/images/mAMRite_small.png)  

(Abandoned for obvious trademark issues and the fact that the joke may be lost on non-Brits)  

## Contributors
[Adam Witney](https://gitlab.com/awitney)  
[Alex Manuele](https://gitlab.com/alexmanuele)  
[Inês Mendes](https://gitlab.com/cimendes)  
[Thanh Le Viet](https://gitlab.com/thanhleviet)  
[Trestan Pillonel](https://gitlab.com/tpillone)  
[Varun Shamanna](https://gitlab.com/varunshamanna4)  

## Introduction

This project originated from the dilemma a scientist faces when choosing a database that stores antimicrobial resistance determinants. Multiple databases exist with comparative strengths and weaknesses. This project builds on the concepts of the [haAMRonization](https://github.com/pha4ge/hAMRonization) project aiming to aggeregate and combine the information contained within the metadata associated with each project. The problem is exacerbated by the fact that the equivalent antimicrobial resistance genes (ARGs) can be named differently in each database.

The hypothesis for the project is as follows:  

* given a match in one database
* find the matches in other databases
* aggregate the combined descriptive information pertaining to antimicrobial resistance contained in the union of the metadata
* report this to user for them to make intelligent informed choices

## Methodology

* Download sequences and associated metadata of ARGs from 3 databases
  * [CARD](https://card.mcmaster.ca/) ([Manuscript](http://www.ncbi.nlm.nih.gov/pubmed/31665441))
  * [NCBI AMR Reference Gene Catalog](https://www.ncbi.nlm.nih.gov/pathogens/refgene/)
  * [Resfinder 4](https://bitbucket.org/genomicepidemiology/resfinder/src/4.0/) ([Manuscipt](https://academic.oup.com/jac/article/75/12/3491/5890997))  
   Details can be found in the [appendices](/docs/appendices.md#data-download)
* Parse the data to
  * extract the protein sequences and write into fasta format with the gene identifiers as the record ids.
  * extract the associated metadata and convert to a consistent `JSON` format  
   Details can be found in the [appendices](/docs/appendices.md#data-parsing)
* Find best matches of each gene from one source database against the other two target databases
  * Where a reciprocal best hit (RBH) exists, report this.  
     Details can be found in the [appendices](/docs/appendices.md#analyse-for-reciprocal-best-hits-rbhs).  
     A summary of the results can be found [here](/docs/appendices.md#summary-of-rbh-analysis)
  * If a RBH does not exist, report the best match as long as thresholds for coverage and indentity are met.
    A summary of the results can be found [here](/docs/appendices.md#summary-of-non-rbh-searches)
  
  For this purpose the [MMseqs2](https://pubmed.ncbi.nlm.nih.gov/29035372/) search tool was used that in its most sensitive mode is 100x faster than blastp and almost as sensitive. In a [comparative manuscript](https://pubmed.ncbi.nlm.nih.gov/33099302/) demonstrated that even in the worst cases MMseqs2 would not miss more than 10% of the RBH produced by blastp. MMseqs2 also contains a convenient wrapper to perform the all-by-all search necessary to find RBHs.

* From the outputs of the MMseqs2 searches the RBHs or best matches of each gene from one database against the other two databases can be parsed to produce a `Directed Graph`. This network was constructed using the [networkx](https://networkx.org/) python package.  
  Details of the method can be found [here](docs/appendices.md#building-a-networkx-graph)  
  In this graph
  * the nodes represent a protein from one database
    * Node attributes contain the phenotype from the JSON metadata
  * the edges link nodes and represent the matches and attributes include
    * type, either RBH or OWH (one way hit)
    * coverage, (alignment length/query length)
    * identity, (percent identity of match)  
  See the image below for a pictoral example using made up data  


\
\
![network diagram](docs/images/chamredb_network.png)

## Current data release

The graph shipped with this version of the package is the release of **2026-11-09**. It was built from

| Database | Release |
| --- | --- |
| CARD | 4.0.2 (sequences and ontology from the same release) |
| NCBI AMRFinderPlus | database 4.2, 2026-10-29.1 |
| ResFinder | db 2.6.0, commit `eecf0aa` |

using MMseqs2 18.8cc5c (reciprocal best hits `easy-rbh -s 7.5`, one-way hits `search -s 7.0 --max-accept 3`). The graph is built by `builder/` (see the [builder tests](builder/tests/README.md)); the exact input files, with source URL and sha256, are listed in the `data_manifest.json` of the release.

| | Count |
| --- | --- |
| Nodes (proteins with at least one match) | 15,859 (NCBI 7,299, CARD 5,971, ResFinder 2,589) |
| Directed edges | 51,971 (RBH 21,428, one-way hits 30,543) |
| Reciprocal best hit pairs | 10,714 (CARD-NCBI 5,685, NCBI-ResFinder 2,538, CARD-ResFinder 2,491) |

The version of the data a result came from is printed by `chamredb query` (`data release:`) and written to the last column (`data_release`) of the TSV output. The earlier 2021 graph (11,718 nodes, 33,632 edges) can be rebuilt with `builder/reproduce_2021/`.

## Assessing the graph

_N.B. The assessment in this section was made on an earlier version of the graph, before the 2026 data refresh. Its counts, figures and CSV files have **not** been regenerated for the current release._

In order to look at the 19132 matches within the database and assess the effectiveness of the methodology the database names for matches were compared with the [Normalized Levensthein algorithm](https://devopedia.org/levenshtein-distance#qst-ans-3).
Before calculating the name similarity between the source and target of a match the name was cleaned using the following steps

1. Removing species names from database names (exclusively in the CARD database) e.g `Staphylococcus aureus mupB conferring resistance to mupirocin`
2. Coversion to lower case
3. Removing the bla prefix
4. Removing parentheses
5. Removing hyphens

_N.B blaPAO-N and blaPDC matches are the source of 562 low name similarities so were skipped_

The resulting data is plotted below showing
A: distribution of levenshtein smilarities between the database names of the best matches
B: distribution sequence identities for the best matches
C: plot of the levenshtein smilarity versus the sequence identity for each match  

![analysis plots](scripts/seq_id_and_name_sim_analysis.png)  

The red line shows the correlation including 95% confidence intervals.

Based on this regression the expected name similarity for a sequence identity of 0.95 can be calculated (0.69)

```
linear_fit = np.polyfit(
    distance_dataframe['sequence_identity'],
    distance_dataframe['name_similarity'],
    1
)
np.polyval(linear_fit, 0.95)
0.6898250936869554
```
To examine data matches where the sequence identity is > 0.95 **BUT** the name similarity is less than the predicted 0.69 was created and explored, a [CSV file](scripts/low_name_similarity.csv) was created.

In this data, many of the differences in the names are due to matches with the same gene family but different alleles e.g `blaADC-125` in the `ncbi` database and `blaADC-25` in the `resfinder 4` database.  
Therefore data calculating name smilarities ignoring alleles was created.


A second series of plots explores this data by plotting the distributions of the name smilarities.
In the top panel the violin plots show the distribution of the name similarity differences for those matches where the sequence identity is greater than 0.95. In A: this name similarities are based on complete cleaned locus names and in B: they are based on names where the alleles are ignored. 
The lower panel of the figure contains violin plots showing the distribution of the **difference** between the **observed name similarity** and that **predicted** by a linear regression model fitting name smilarity to sequence identity. The right hand 2 plots are data where the name simialrities were calculated excluding alleles.
![analysis plots](scripts/name_similarity_distribution_analysis.png)

The data was filtered for those matches where the sequence identity is > 0.95 but the name similarity is less than the predicted value of 0.86 based on the linear regression model. 

To examine these anomalous results a [CSV sheet](scripts/low_name_similarity_without_alleles.csv) was created.  

Exploring this data some of these are clearly related genes but the databases have different nomneclature e.g  
`vanA` in `card` and `VanHAX` in `resfinder 4.0` or  
`catA15` in `ncbi` or `Clostridium butyricum catB` in `card`

_**N.B** The species names are removed in the name cleaning function._  

In other cases the names are completely different, e.g  
`gimA` in `card` and `mgt`  in `ncbi` but the sequences are 99.5% identical. `gimA` is a macrolide glycosyltransferase and may confer resistance to spiramycin. `mgt` in the `ncbi` database stands for `macrolide-inactivating glycosyltransferase`. Clearly the genes are likely to have the same function but have been assigned different names in the two databases.

## Querying the graph

```
usage: chamredb query
  [-h]
  [-d {card,ncbi,resfinder}]
  [-ct COVERAGE_THRESHOLD] [-it IDENTITY_THRESHOLD]
  [--graph GRAPH]
  (-i ID | -f ID_FILE | -j HAMRONIZATION_JSON_FILE)
  [-o OUTFILE_PATH]
```

`--graph` uses a graph JSON file instead of the one shipped with the package (for example a newer release downloaded from the releases page), which also works offline.  

The graph can be queried in one of 3 ways  

### 1. Querying an individual

Requires specifying the identifier `-i` and database `-d`  

```
chamredb query -d ncbi -i WP_012695489.1 
```

Alternatively the gene name can be used

```
chamredb query -d ncbi -i qnrB2
```

The output reports the matches and metadata from the other databases  
![qnrB2](/docs/images/qnrB2.png)

Another example where the matches are one way hits not RBHs

```
chamredb query -d resfinder -i "aac(3)-IIIb"
```

![aac(3)-IIIb](/docs/images/aac(3)-IIIb.png)
In these outputs ↔ means a RBH, and ➡ a search hit

### 2. Providing a list of identifiers from a single database

Requires specifying the database `-d`, the text file containing the ids `-f`, and a path for the tsv output file `-o`  

```
chamredb query -d card -f docs/data/card_ids.txt  -o docs/data/card_ids.tsv
```

This will produce a [TSV file](/docs/data/card_ids.tsv) containing the matches and associated metadata with one row per id in the text file

### 3. Use hAMRonization summary output

Use the [hAMRonization softare](https://github.com/pha4ge/hAMRonization) to convert the outputs from antimicrobial resistance gene detection tools into a unified format. Concatenate and summarize AMR detection reports into a single summary JSON file using the `hamronize summarize` command from this package. The JSON output from this step can be used to query ChamreDb.  
Use `-j` to specify the json summary file and `-o` the path for the TSV output  
**Please Note**  
Only outputs using data derived from AMR detection tools that have searched either the `CARD`, `NCBI` or `Resfinder 4` databases can be used.

```
chamredb query -j docs/data/hamronize_summary.json -o docs/data/hamronize_summary.tsv
```

This will produce a [TSV file](/docs/data/hamronize_summary.tsv) containing the matches and associated metadata with one row per id in the text file
