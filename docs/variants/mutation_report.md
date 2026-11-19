# Mutation comparison report (data release 2026-11-19)

Built by `builder/scripts/compare_mutations.py`. Positions are in the numbering of each locus's anchor database (CARD, else NCBI, else PointFinder); rRNA positions are on the sense strand of the anchor.

## Inputs

| Database | Role | Release |
| --- | --- | --- |
| CARD | mutations | 4.0.2 |
| NCBI | mutations | AMRFinderPlus database 4.2, 2026-10-29.1 |
| PointFinder | mutations | 3bd7d5ef2dbea88e7d40e19126dc587f5abb8e4b |

## Invariants

All hold.

## Rows and records

| Database | Upstream rows | Parsed | Excluded | Records | Single | Combination | Organism groups |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CARD | 3328 | 3324 | 4 | 3324 | 3231 | 93 | 64 |
| NCBI | 1637 | 1636 | 1 | 1636 | 1636 | 0 | 31 |
| PointFinder | 1338 | 1286 | 52 | 1976 | 1608 | 368 | 10 |

Records are distinct changes: rows listing the same change twice (CARD once per evidence source, PointFinder once per phenotype) are merged, and rows listing several alternative residues are split.

### Excluded rows

| Database | Reason | Rows |
| --- | --- | --- |
| CARD | CARD insertion between non-adjacent residues: 'N330_R335insTIPY' | 1 |
| CARD | CARD mutation with a residue that is not an amino acid: 'Q10Z' | 2 |
| CARD | CARD mutation with a residue that is not an amino acid: 'Q142X' | 1 |
| NCBI | protein WP_004188318.1 not in AMRProt | 1 |
| PointFinder | PointFinder required mutation not understood: 'rpoB_1_AP020316.1_ins475' | 1 |
| PointFinder | any premature stop codon in the gene (no position) | 4 |
| PointFinder | no sequence for gyrA_2-CP035927.1 | 1 |
| PointFinder | not a bacterium | 46 |

## Reference checks

Self check: the residue each database states as the original one, against its own reference sequence. Anchor check: the same residue against the anchor sequence at the converted position.

| Database | Self ok | Self mismatch | Mismatch share | Unchecked | Anchor ok | Anchor mismatch | Not placed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CARD | 3135 | 227 | 6.8% | 64 | 3133 | 214 | 15 |
| NCBI | 1633 | 3 | 0.2% | 0 | 1598 | 8 | 1 |
| PointFinder | 2308 | 15 | 0.6% | 80 | 2213 | 44 | 1 |

### Self-check mismatches explained by one constant shift

A shift that makes at least three mismatching residues agree suggests the database numbers this gene in another convention (often *E. coli* numbering). Not corrected: listed for review.

| Organism | Gene | Changes | Mismatches | Shift | Explained |
| --- | --- | --- | --- | --- | --- |
| Streptococcus pneumoniae R6 | PBP2b | 3 | 3 | +6 | 3 |
| Staphylococcus aureus subsp. aureus MRSA252 | gyrB | 7 | 6 | -1 | 6 |
| Bartonella bacilliformis KC583 | gyrA | 4 | 4 | -40 | 3 |
| Escherichia coli CFT073 | EFTu | 6 | 5 | +15 | 5 |
| Escherichia coli str. K-12 substr. MG1655 | EFTu | 4 | 4 | +1 | 4 |
| Mycolicibacterium smegmatis MC2 155 | 16rrsA | 3 | 3 | +1 | 3 |
| Helicobacter pylori | 23S | 14 | 7 | +4 | 7 |
| Streptococcus pneumoniae | 23S | 4 | 3 | -63 | 3 |
| Clostridium ljungdahlii DSM 13528 | gyrB | 20 | 20 | +2 | 12 |
| Thermus thermophilus HB27 | 23S | 6 | 4 | +69 | 4 |
| Bacillus subtilis subsp. subtilis str. 168 | rpsE | 12 | 11 | +19 | 5 |

### rRNA references not stored or numbered along the sense strand

The stored strand is found from conserved rRNA stretches. Anchor positions of rRNA are always on the sense strand (the rRNA itself). Positions are read along the sense strand unless the stated original residues fit the antisense strand at 3 more sites; a reference with unknown strand is read as stored.

| Locus | Database | Reference | Stored | Positions read along | Sites | Fit sense | Fit antisense |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Borreliella burgdorferi:16S rRNA:rna | CARD | ARO:3003503 | unknown | unknown | 1 |  |  |
| Brachyspira hyodysenteriae:23S rRNA:rna | CARD | ARO:3004133 | antisense | sense | 1 | 0 | 0 |
| Chlamydomonas reinhardtii:23S rRNA:rna | CARD | ARO:3004132 | antisense | sense | 2 | 1 | 0 |
| Clostridioides_difficile:23S rRNA:rna | NCBI | NZ_CP019870.1@23S_ribosomal_RNA@23S:4042221-4045120 | antisense | sense | 1 | 0 | 1 |
| Helicobacter_pylori:16S rRNA:rna | CARD | ARO:3003510 | antisense | sense | 5 | 3 | 2 |
| Moraxella catarrhalis:23S rRNA:rna | CARD | ARO:3004138 | antisense | sense | 1 | 1 | 0 |
| Mycolicibacterium smegmatis:16S rRNA:rna | CARD | ARO:3003539 | antisense | sense | 3 | 0 | 0 |
| Mycolicibacterium smegmatis:16S rRNA:rna | CARD | ARO:3003543 | antisense | sense | 2 | 0 | 0 |
| Mycolicibacterium smegmatis:16S rRNA:rna | CARD | ARO:3003544 | antisense | sense | 1 | 0 | 1 |
| Mycolicibacterium smegmatis:16S rRNA:rna | CARD | ARO:3003546 | antisense | sense | 1 | 0 | 0 |
| Mycoplasma hominis:23S rRNA:rna | CARD | ARO:3004176 | antisense | sense | 2 | 0 | 0 |
| Mycoplasmopsis fermentans:23S rRNA:rna | CARD | ARO:3004177 | antisense | sense | 1 | 0 | 1 |
| Neisseria_gonorrhoeae:16S rRNA:rna | NCBI | NZ_AP019853.2@16S_ribosomal_RNA@16S:1155759-1157299 | antisense | sense | 1 | 1 | 1 |
| Salmonella:16S rRNA:rna | CARD | ARO:3003512 | antisense | sense | 2 | 0 | 1 |
| Salmonella:16S rRNA:rna | NCBI | NC_003197.2@16S_ribosomal_RNA@16S_rrsD:3570470-3572005 | antisense | sense | 2 | 0 | 2 |
| Streptococcus_pneumoniae:23S rRNA:rna | CARD | ARO:3004181 | antisense | sense | 1 | 1 | 0 |
| Streptococcus_pneumoniae:23S rRNA:rna | NCBI | NZ_CP018347.1@23S_ribosomal_RNA@23S:32347-35250 | antisense | antisense | 4 | 1 | 4 |

## Loci

486 loci: 31 with references from all three databases, 141 from two, 314 from one. 20 pairs of differently named references were joined because their sequences are near-identical.

| Group | Molecule | Reference | Joined with | Identity | Coverage |
| --- | --- | --- | --- | --- | --- |
| Enterococcus_faecium | protein | card|ARO:3003761|protein | ncbi|WP_002296175.1|protein | 0.998 | 1.0 |
| Escherichia | protein | card|ARO:3003368|protein | ncbi|WP_000031783.1|protein | 1.0 | 1.0 |
| Escherichia | protein | card|ARO:3007423|protein | ncbi|WP_000642196.1|protein | 1.0 | 1.0 |
| Escherichia | protein | ncbi|WP_094970923.1|protein | pointfinder|escherichia_coli/pmrA_1_CP047010.1|protein | 0.9955 | 1.0 |
| Escherichia | rna | card|ARO:3003333|rna | ncbi|NC_000913.3@16S_ribosomal_RNA@16S_rrsB:4166659-4168201|rna | 0.9961 | 1.0 |
| Escherichia | rna | card|ARO:3003372|rna | ncbi|NC_000913.3@16S_ribosomal_RNA@16S_rrsB:4166659-4168201|rna | 0.9961 | 0.9994 |
| Escherichia | rna | card|ARO:3003376|rna | ncbi|NC_000913.3@16S_ribosomal_RNA@16S_rrsB:4166659-4168201|rna | 1.0 | 1.0 |
| Haemophilus_influenzae | protein | card|ARO:3004446|protein | ncbi|WP_005693446.1|protein | 1.0 | 1.0 |
| Helicobacter_pylori | protein | card|ARO:3007060|protein | ncbi|WP_000913670.1|protein | 0.9803 | 1.0 |
| Klebsiella_pneumoniae | protein | card|ARO:3007421|protein | ncbi|WP_016532398.1|protein | 1.0 | 1.0 |
| Neisseria_gonorrhoeae | protein | card|ARO:3000464|protein | ncbi|WP_106333382.1|protein | 0.9302 | 1.0 |
| Neisseria_gonorrhoeae | protein | card|ARO:3000464|protein | pointfinder|neisseria_gonorrhoeae/porB_1_AE004969.1|protein | 0.9684 | 1.0 |
| Neisseria_gonorrhoeae | protein | card|ARO:3004832|protein | ncbi|WP_003693546.1|protein | 0.9983 | 1.0 |
| Neisseria_meningitidis | protein | card|ARO:3003937|protein | ncbi|WP_002218783.1|protein | 1.0 | 1.0 |
| Pseudomonas_aeruginosa | protein | card|ARO:3003582|protein | ncbi|WP_003102236.1|protein | 1.0 | 1.0 |
| Pseudomonas_aeruginosa | protein | card|ARO:3004059|protein | ncbi|WP_003094842.1|protein | 1.0 | 1.0 |
| Pseudomonas_aeruginosa | protein | card|ARO:3004060|protein | ncbi|WP_003094842.1|protein | 1.0 | 1.0 |
| Staphylococcus_aureus | protein | card|ARO:3003314|protein | pointfinder|staphylococcus_aureus/grlB_1_CP026964.1|protein | 0.9985 | 1.0 |
| Staphylococcus_aureus | protein | ncbi|WP_001289586.1|protein | pointfinder|staphylococcus_aureus/grlA_1_CP026964.1|protein | 0.9938 | 1.0 |
| mycobacterium_tuberculosis | rna | card|ARO:3007536|rna | pointfinder|mycobacterium_tuberculosis/rrs_1_CP054014.1|rna | 1.0 | 1.0 |

## Organism coverage

| Organism group | CARD | NCBI | PointFinder |
| --- | --- | --- | --- |
| Acinetobacter_baumannii | 8 | 53 | 0 |
| Burkholderia_cepacia | 3 | 21 | 0 |
| Campylobacter | 3 | 51 | 36 |
| Clostridioides_difficile | 21 | 71 | 0 |
| Enterobacter_cloacae | 2 | 12 | 0 |
| Enterococcus_faecalis | 19 | 26 | 14 |
| Enterococcus_faecium | 24 | 72 | 44 |
| Escherichia | 181 | 322 | 229 |
| Haemophilus_influenzae | 5 | 41 | 0 |
| Helicobacter_pylori | 173 | 13 | 2 |
| Klebsiella_pneumoniae | 17 | 127 | 99 |
| Neisseria_gonorrhoeae | 48 | 89 | 48 |
| Neisseria_meningitidis | 10 | 46 | 0 |
| Pseudomonas_aeruginosa | 53 | 234 | 0 |
| Salmonella | 28 | 96 | 184 |
| Staphylococcus_aureus | 115 | 259 | 140 |
| Streptococcus_pneumoniae | 20 | 23 | 0 |
| Streptococcus_pyogenes | 4 | 3 | 0 |
| Vibrio_vulnificus | 1 | 3 | 0 |
| mycobacterium_tuberculosis | 2307 | 0 | 1180 |

Only organism groups with mutations in at least two databases are shown (20 of 76).

## Match classes (single mutations)

| Class | card->ncbi | card->pointfinder | ncbi->card | ncbi->pointfinder | pointfinder->card | pointfinder->ncbi |
| --- | --- | --- | --- | --- | --- | --- |
| same_mutation | 4 | 27 | 4 | 130 | 22 | 130 |
| same_mutation_different_reference | 394 | 720 | 381 | 147 | 653 | 151 |
| possible_same_mutation_offset | 34 | 7 | 30 | 4 | 18 | 6 |
| same_mutation_only_in_combination | 0 | 16 | 14 | 109 | 9 | 0 |
| same_position_different_change | 13 | 447 | 67 | 41 | 169 | 6 |
| different_mutation | 100 | 1117 | 316 | 177 | 385 | 101 |
| same_mutation_related_organism | 11 | 33 | 43 | 31 | 10 | 0 |
| locus_not_covered | 158 | 495 | 718 | 442 | 342 | 34 |
| organism_not_covered | 2517 | 369 | 63 | 555 | 0 | 1180 |

## Match classes (combinations)

| Class | card->ncbi | card->pointfinder | pointfinder->card | pointfinder->ncbi |
| --- | --- | --- | --- | --- |
| same_combination | 0 | 7 | 7 | 0 |
| partial_combination | 0 | 1 | 50 | 0 |
| components_listed_separately | 13 | 17 | 14 | 329 |
| combination_not_listed | 15 | 44 | 271 | 39 |
| locus_not_covered | 4 | 17 | 26 | 0 |
| organism_not_covered | 61 | 7 | 0 | 0 |

## Agreement where both databases cover the organism

| Pair | Mutations (clusters) in shared organism groups | Listed by both | Share |
| --- | --- | --- | --- |
| CARD and NCBI | 1859 | 382 | 20.5% |
| CARD and PointFinder | 3706 | 663 | 17.9% |
| NCBI and PointFinder | 1194 | 269 | 22.5% |

## Mutations per organism group and databases listing them

| Group | card | card+ncbi | card+ncbi+pointfinder | card+pointfinder | ncbi | ncbi+pointfinder | pointfinder |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Acinetobacter_baumannii | 2 | 6 | 0 | 0 | 47 | 0 | 0 |
| Burkholderia_cepacia | 0 | 3 | 0 | 0 | 18 | 0 | 0 |
| Campylobacter | 0 | 0 | 3 | 0 | 25 | 15 | 4 |
| Clostridioides_difficile | 11 | 9 | 0 | 0 | 62 | 0 | 0 |
| Enterobacter_cloacae | 2 | 0 | 0 | 0 | 12 | 0 | 0 |
| Enterococcus_faecalis | 12 | 7 | 0 | 0 | 7 | 12 | 2 |
| Enterococcus_faecium | 10 | 14 | 0 | 0 | 38 | 20 | 24 |
| Escherichia | 26 | 96 | 46 | 0 | 155 | 25 | 5 |
| Haemophilus_influenzae | 0 | 5 | 0 | 0 | 36 | 0 | 0 |
| Helicobacter_pylori | 148 | 6 | 0 | 2 | 7 | 0 | 0 |
| Klebsiella_pneumoniae | 3 | 13 | 1 | 0 | 106 | 7 | 91 |
| Neisseria_gonorrhoeae | 11 | 14 | 19 | 1 | 45 | 11 | 9 |
| Neisseria_meningitidis | 8 | 2 | 0 | 0 | 44 | 0 | 0 |
| Pseudomonas_aeruginosa | 31 | 19 | 0 | 0 | 215 | 0 | 0 |
| Salmonella | 4 | 14 | 6 | 1 | 40 | 36 | 5 |
| Staphylococcus_aureus | 20 | 64 | 23 | 0 | 127 | 45 | 3 |
| Streptococcus_pneumoniae | 10 | 10 | 0 | 0 | 13 | 0 | 0 |
| Streptococcus_pyogenes | 2 | 2 | 0 | 0 | 1 | 0 | 0 |
| Vibrio_vulnificus | 1 | 0 | 0 | 0 | 3 | 0 | 0 |
| mycobacterium_tuberculosis | 1648 | 0 | 0 | 561 | 0 | 0 | 619 |

## Drug classes of matched mutations

| Pair | same_classes | shared_class | no_shared_class | unknown |
| --- | --- | --- | --- | --- |
| card->ncbi | 215 | 137 | 14 | 32 |
| card->pointfinder | 597 | 52 | 19 | 79 |
| ncbi->card | 210 | 130 | 14 | 31 |
| ncbi->pointfinder | 253 | 10 | 0 | 14 |
| pointfinder->card | 550 | 44 | 3 | 78 |
| pointfinder->ncbi | 257 | 10 | 0 | 14 |

## Systematic offsets between databases

The same change at a constant distance on the anchor, seen at least three times: probably a numbering convention, not different mutations. Review in `offset_review.tsv`.

| Locus | Database | Other | Offset | Count |
| --- | --- | --- | --- | --- |
| Escherichia:EFTu:protein | CARD | NCBI | +1 | 4 |
| Escherichia:EFTu:protein | CARD | NCBI | +15 | 6 |
| Escherichia:EFTu:protein | NCBI | CARD | -15 | 6 |
| Neisseria_gonorrhoeae:23S rRNA:rna | PointFinder | CARD | -6 | 4 |
| Neisseria_gonorrhoeae:23S rRNA:rna | PointFinder | CARD | -2 | 4 |
| Staphylococcus_aureus:gyrB:protein | CARD | NCBI | -1 | 6 |
| Staphylococcus_aureus:gyrB:protein | NCBI | CARD | +1 | 6 |
| Staphylococcus_aureus:mprF:protein | CARD | NCBI | -4 | 3 |
| Staphylococcus_aureus:mprF:protein | NCBI | CARD | +4 | 3 |
| Streptococcus_pneumoniae:PBP2b:protein | CARD | NCBI | +6 | 3 |
| Streptococcus_pneumoniae:PBP2b:protein | NCBI | CARD | -6 | 3 |

