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
[Finlay Maguire](https://github.com/fmaguire)

## Introduction

This project originated from the dilemma a scientist faces when choosing a database that stores antimicrobial resistance determinants. Multiple databases exist with comparative strengths and weaknesses. This project builds on the concepts of the [hAMRonization](https://github.com/pha4ge/hAMRonization) project aiming to aggregate and combine the information contained within the metadata associated with each project. The problem is exacerbated by the fact that the equivalent antimicrobial resistance genes (ARGs) can be named differently in each database.

The hypothesis for the project is as follows:  

* given a match in one database
* find the matches in other databases
* aggregate the combined descriptive information pertaining to antimicrobial resistance contained in the union of the metadata
* report this to user for them to make intelligent informed choices

## At a glance

![Real families: the same protein under different names in CARD, NCBI and ResFinder](docs/images/figure_real_families.png)

The same protein often carries a different name in each database (panels B to D), and many proteins exist in only some of them. Of the <!-- stat:groups.total -->7,589<!-- /stat --> groups of proteins linked by reciprocal best hits, only <!-- stat:groups.all_three -->2,406<!-- /stat --> (<!-- stat:groups.all_three_share -->32%<!-- /stat -->) contain all three databases.

![Overlap of the equivalence groups between CARD, NCBI and ResFinder](docs/images/figure_overlap_upset.png)

## Methodology

* Download sequences and associated metadata of ARGs from 3 databases
  * [CARD](https://card.mcmaster.ca/) ([Manuscript](http://www.ncbi.nlm.nih.gov/pubmed/31665441))
  * [NCBI AMR Reference Gene Catalog](https://www.ncbi.nlm.nih.gov/pathogens/refgene/)
  * [Resfinder 4](https://bitbucket.org/genomicepidemiology/resfinder/src/4.0/) ([Manuscript](https://academic.oup.com/jac/article/75/12/3491/5890997))  
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
    * coverage, (alignment length/query length, a fraction that can exceed 1)
    * identity, (fraction of identical residues in the match, 0 to 1)  

![Three proteins, one from each database, with the metadata stored on each node and the type, coverage and identity stored on each edge. CARD and NCBI are reciprocal best hits (RBH, one edge each way); NCBI also has a one-way hit (OWH) to ResFinder.](docs/images/figure_graph_example.png)

The example is read from the release graph (CARD ARO:3000579, NCBI WP_029728367.1, ResFinder blaB-11). Only a subset of the stored fields is shown, long lines are wrapped, and both RBH edges carry the same values. CARD also has a one-way hit to ResFinder (identity 0.976), which is not drawn. This pair is a clean illustration; many edges have lower identity or coverage.

![Pipeline from the three databases to the released graph](docs/images/figure_pipeline.png)

The numbers in this figure are read from the release's manifest and statistics, so they match the graph described below.

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
| Nodes (proteins with at least one match) | <!-- stat:graph.nodes -->15,859<!-- /stat --> |
| Directed edges | <!-- stat:graph.edges -->51,971<!-- /stat --> (reciprocal best hits <!-- stat:graph.edges_rbh -->21,428<!-- /stat -->, one-way hits <!-- stat:graph.edges_owh -->30,543<!-- /stat -->) |
| Reciprocal best hit pairs | <!-- stat:graph.rbh_pairs -->10,714<!-- /stat --> |
| Equivalence groups (proteins joined by reciprocal best hits) | <!-- stat:groups.total -->7,589<!-- /stat --> |

The version of the data a result came from is printed by `chamredb query` (`data release:`) and written to the last column (`data_release`) of the TSV output. The earlier 2021 graph (11,718 nodes, 33,632 edges) can be rebuilt with `builder/reproduce_2021/`. How the graph grew from 2021 to this release, step by step, is in [docs/release-2026-11-09.md](docs/release-2026-11-09.md).

## Assessing the graph

How well do database names agree where the sequences match? For every protein the best match in each other database was found (the match with the highest coverage and identity, ties broken by edge order), and the two names were compared with the [normalised Levenshtein similarity](https://devopedia.org/levenshtein-distance#qst-ans-3), where 1 means identical names. Before the comparison each name was cleaned:

1. Removing species names from database names (CARD only), for example `Staphylococcus aureus mupB conferring resistance to mupirocin` becomes `mupB`
2. Conversion to lower case
3. Removing a leading `bla`
4. Removing parentheses
5. Removing hyphens

For the analysis that ignores alleles a trailing allele number is also removed (`qnrB20` becomes `qnrb`).

<!-- stat:analysis.matches -->30,152<!-- /stat --> best matches were analysed. A further <!-- stat:analysis.skipped -->1,349<!-- /stat --> were skipped because their target is the ResFinder entry `blaPAO`, a single entry that is the best ResFinder match of every one of the PDC alleles in CARD and NCBI, so comparing its name says nothing (the blaPDC alleles themselves are kept). Of the analysed matches, <!-- stat:analysis.high_identity -->24,637<!-- /stat --> have a sequence identity of at least 0.95.

![Name similarity against sequence identity](docs/images/figure_name_similarity_vs_identity.png)

Names agree more where sequences agree more, but not always. A linear regression of name similarity on sequence identity predicts a similarity of <!-- stat:analysis.full.predicted -->0.84<!-- /stat --> at identity 0.95 when names are compared in full, and <!-- stat:analysis.free.predicted -->0.93<!-- /stat --> when alleles are ignored (each analysis uses its own regression). Among the matches with identity of at least 0.95, <!-- stat:analysis.full.low_share -->16.9%<!-- /stat --> (<!-- stat:analysis.full.low_rows -->4,153<!-- /stat -->) have a name similarity at or below the full-name prediction, and <!-- stat:analysis.free.low_share -->2.4%<!-- /stat --> (<!-- stat:analysis.free.low_rows -->597<!-- /stat -->) when alleles are ignored. Most of the disagreement is therefore about alleles: the same gene family, a different allele number.

![Distributions of name similarity and of the difference from the prediction](docs/images/figure_name_distributions.png)

The tables of these matches are [docs/analysis/low_name_similarity.csv](docs/analysis/low_name_similarity.csv) and [docs/analysis/low_name_similarity_without_alleles.csv](docs/analysis/low_name_similarity_without_alleles.csv).

Some disagreements are not about alleles. These examples are proteins whose sequences are identical (identity of at least 0.999) and whose names share nothing once alleles are ignored; the caption explains how they were chosen and that they are illustrative, not representative.

![Examples of identical proteins with unrelated names](docs/images/figure_name_disagreements.png)

All figures and tables in this section are made from the release graph by `python -m analysis.run --graph chamredb/data/graph.json --out docs --manifest docs/data_manifest.json`, with the data behind every figure in `docs/images/source_data_*`.

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
![chamredb query -d ncbi -i qnrB2](docs/images/figure_query_qnrB2.svg)

Another example where the matches are one way hits not RBHs

```
chamredb query -d resfinder -i "aac(3)-IIIb"
```

![chamredb query -d resfinder -i aac(3)-IIIb](docs/images/figure_query_aac3IIIb.svg)
In these outputs ↔ means a RBH, and ➡ a search hit

### 2. Providing a list of identifiers from a single database

Requires specifying the database `-d`, the text file containing the ids `-f`, and a path for the tsv output file `-o`  

```
chamredb query -d card -f docs/data/card_ids.txt  -o docs/data/card_ids.tsv
```

This will produce a [TSV file](/docs/data/card_ids.tsv) containing the matches and associated metadata with one row per id in the text file

### 3. Use hAMRonization summary output

Use the [hAMRonization software](https://github.com/pha4ge/hAMRonization) to convert the outputs from antimicrobial resistance gene detection tools into a unified format. Concatenate and summarize AMR detection reports into a single summary JSON file using the `hamronize summarize` command from this package. The JSON output from this step can be used to query ChamreDb.  
Use `-j` to specify the json summary file and `-o` the path for the TSV output  
**Please Note**  
Only outputs using data derived from AMR detection tools that have searched either the `CARD`, `NCBI` or `Resfinder 4` databases can be used.

```
chamredb query -j docs/data/hamronize_summary.json -o docs/data/hamronize_summary.tsv
```

This will produce a [TSV file](/docs/data/hamronize_summary.tsv) containing the matches and associated metadata with one row per id in the text file

## Comparing resistance mutations

Some resistance comes from mutations in genes every isolate carries (gyrA, rpoB, 23S rRNA, promoters) rather than from acquired genes. CARD lists these in its variant models, NCBI in the AMRFinderPlus point-mutation files, and ResFinder in [PointFinder](https://bitbucket.org/genomicepidemiology/pointfinder_db). `chamredb variants` compares the three.

Data release <!-- vstat:release -->2026-11-19<!-- /vstat --> has <!-- vstat:records -->6,936<!-- /vstat --> mutations: CARD <!-- vstat:card -->3,324<!-- /vstat -->, NCBI <!-- vstat:ncbi -->1,636<!-- /vstat --> and PointFinder <!-- vstat:pointfinder -->1,976<!-- /vstat -->, <!-- vstat:combinations -->461<!-- /vstat --> of them combinations of mutations, in <!-- vstat:groups -->76<!-- /vstat --> organism groups. It was built from the CARD and NCBI releases of the graph above and pointfinder_db commit `3bd7d5e` (inputs in [docs/variants/data_manifest.json](docs/variants/data_manifest.json)).

### Method

Organisms are grouped by the organisms the tools can call mutations for: the AMRFinderPlus `--organism` labels (for example `Escherichia`, which covers *E. coli* and *Shigella*, and `Klebsiella_pneumoniae`), then the PointFinder species (`mycobacterium_tuberculosis`), then the CARD species for organisms neither tool covers. The mapping is in [builder/variants/taxa.tsv](builder/variants/taxa.tsv).

Within an organism group, references from different databases are the same locus when their gene names agree or their sequences are near-identical (identity and coverage of at least 0.9). Each locus has one anchor sequence, CARD's if CARD has the gene, else NCBI's, else PointFinder's. Each database's positions are converted to the anchor by a global alignment of its reference sequence to the anchor, and positions are reported in the anchor's numbering. Records keep the database's own notation too. Promoter positions are offsets from the gene start in all three databases and are compared directly.

The original residue each database gives is checked against its own reference sequence. It disagrees for <!-- vstat:self_ncbi -->0.2%<!-- /vstat --> of NCBI's mutations, <!-- vstat:self_pointfinder -->0.6%<!-- /vstat --> of PointFinder's and <!-- vstat:self_card -->6.8%<!-- /vstat --> of CARD's (<!-- vstat:self_card_n -->227<!-- /vstat --> residues). Most of the CARD mismatches are positions in another numbering, often *E. coli* numbering; they are listed in the [mutation report](docs/variants/mutation_report.md).

Each mutation is compared with each other database and given one class:

| Class | Meaning |
| --- | --- |
| `same_mutation` | same change, same position, same reference sequence |
| `same_mutation_different_reference` | same change on different reference sequences that agree once aligned |
| `possible_same_mutation_offset` | same change up to 10 positions away where a reference check fails or the database positions agree; probably a numbering difference, left for review |
| `same_mutation_only_in_combination` | the other database lists it only together with other mutations (PointFinder lists *E. coli* parC S80I only with a gyrA S83 change) |
| `same_position_different_change` | a different change at the same position |
| `different_mutation` | the other database has mutations in the gene but not this one (the nearest is reported) |
| `same_mutation_related_organism` | listed by the other database only for an organism of the same family |
| `locus_not_covered` / `organism_not_covered` | the other database has no mutations in this gene for this organism, or none for the organism |

Combinations are compared as sets: `same_combination`, `partial_combination`, `components_listed_separately` or `combination_not_listed`.

### Results

![Mutations of each organism by the databases that list them](docs/images/figure_variant_coverage.png)

<!-- vstat:groups_two -->20<!-- /vstat --> of the <!-- vstat:groups -->76<!-- /vstat --> organism groups have mutations in more than one database and <!-- vstat:groups_three -->9<!-- /vstat --> in all three. <!-- vstat:card_tb_share -->69%<!-- /vstat --> of CARD's mutations are in *M. tuberculosis*, which AMRFinderPlus does not cover. Of the <!-- vstat:single -->5,183<!-- /vstat --> distinct single mutations, <!-- vstat:single_three -->98<!-- /vstat --> are in all three databases and <!-- vstat:single_one -->4,065<!-- /vstat --> in one only. In the organisms both databases cover, <!-- vstat:share_ncbi_card -->21%<!-- /vstat --> of the mutations are listed by both CARD and NCBI, <!-- vstat:share_card_pointfinder -->18%<!-- /vstat --> by both CARD and PointFinder, and <!-- vstat:share_ncbi_pointfinder -->23%<!-- /vstat --> by both NCBI and PointFinder.

![Match classes of each database's mutations against each other database](docs/images/figure_variant_classes.png)

The figure leaves out mutations of organisms the other database does not cover. Most matches need the conversion between reference sequences, because the databases seldom share one. <!-- vstat:offset_pairs -->66<!-- /vstat --> pairs of mutations may be the same mutation numbered differently, and at <!-- vstat:systematic_offsets -->5<!-- /vstat --> genes the same shift occurs three or more times (CARD numbers *S. aureus* gyrB one higher than NCBI on an identical sequence). They are in [docs/variants/offset_review.tsv](docs/variants/offset_review.tsv). The [mutation report](docs/variants/mutation_report.md) has the excluded rows, the reference checks and all class counts, and [docs/variants/variant_comparison.tsv](docs/variants/variant_comparison.tsv) has one row per mutation with a status for each database.

### Querying mutations

```
usage: chamredb variants
  [-d {card,ncbi,pointfinder}] [-t TAXON] [--variants VARIANTS]
  (-i ID | -l LOCUS | -f ID_FILE | -j HAMRONIZATION_JSON_FILE)
  [-o OUTFILE_PATH]
```

Mutations can be given as the database writes them or in HGVS-like form: `gyrA_S83L`, `"gyrA p.Ser83Leu"`, `"23S A2059G"`, `"fabG1 c.-15C>T"`. Mutations joined by commas are looked up as a combination. `-t` selects the organism (`-t Escherichia`, `-t "Mycobacterium tuberculosis"`).

```
chamredb variants -d ncbi -i gyrA_S83L -t Escherichia
```

prints the mutation in anchor numbering (`p.Ser83Leu` on CARD AAC75291.1) and how CARD and PointFinder list it.

```
chamredb variants -l Escherichia:gyrA
```

prints all mutations of a locus, one column per database. `-f` (one mutation per line, with `-d`) and `-j` (the mutation rows of AMRFinderPlus, RGI and PointFinder in a hAMRonization summary) write a TSV with one row per mutation, the status and notation in each database, and a `data_release` column.

### Limitations

- This compares the database catalogues; it does not say whether a mutation causes resistance in a given isolate.
- PointFinder's rows for any premature stop codon in a gene (*M. tuberculosis* katG, ethA, pncA, gidB) have no position and are left out, so CARD's nonsense mutations in those genes are reported as not listed by PointFinder.
- A mutation written for any residue (CARD `S315Var`) matches the specific changes at that position but is counted as its own mutation.
- Some rRNA genes are stored as the reverse complement of the rRNA (four NCBI references, twelve CARD ones). The stored strand is detected from conserved rRNA stretches, anchor positions are always on the sense strand, and positions are read along the sense strand unless the stated residues fit the antisense strand at three more sites (only NCBI *S. pneumoniae* 23S); see the table in the [report](docs/variants/mutation_report.md). With one or two mutations per reference this choice rests on little evidence, so NCBI *Salmonella* 16S and *C. difficile* 23S, whose residues fit only the stored antisense strand, are read on the sense strand and left for review.
- A numbering offset is only suspected when a reference check fails, so a shifted position that happens to land on the same residue is reported as a different mutation (CARD *S. aureus* gyrB S128L against NCBI S127L).
- PointFinder's `klebsiella` folder is treated as `Klebsiella_pneumoniae`. CARD knockout models and PointFinder *Plasmodium falciparum* are left out.
- Drug classes are compared at a coarse level ([builder/variants/drug_classes.tsv](builder/variants/drug_classes.tsv)).
