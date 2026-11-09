# Scoping: refresh chamredb data and modernise the project

- Date: 2026-11-07
- Requested by: @antunderwood
- Workflow: Maintenance (full), with a Snakemake pipeline merged in (ADR 0002)
- Next command: /bio-spec

## Request

Bring the graph of Source database matches (CARD, NCBI, ResFinder RBH and OWH edges) up to date with current releases, and modernise the code and build, keeping the RBH plus OWH method and existing outputs.

## Why this workflow

Code and behaviour exist and must be preserved, so Maintenance. Runner-up: Pipeline, since the build must be rerunnable; that part is covered by merging the builder (ADR 0002).

## Premise check

Re-checked 2026-11-07: commit 6c2645e moved the build code to `chamredb-builder` (read locally); NCBI upstream is 2026-10-29.1 against 2021-03-01.1 shipped; ARO repo last commit 2026-10-07; ResFinder db last commit 2026-01-26. The data is about 5.5 years stale. Confirmed.

## Prior art

| Candidate | Version | Licence | Maintained? | Fit | Decision |
| --- | --- | --- | --- | --- | --- |
| hAMRonization | 1.3.1 | LGPLv3 | yes | defines the JSON input chamredb consumes | use as input format, test against current output |
| argNorm | 1.1.0 | MIT | yes | maps tool outputs to ARO; different goal (normalise, not cross-database equivalence) | not suitable as replacement, note in README |
| mmseqs2 (bioconda) | 18.8cc5c (2025-07-27; 17.b804f was the previous release) | GPLv3 | yes | the search engine | upgrade in step 2 of ADR 0003 |
| ncbi-amrfinderplus, resfinder (bioconda) | 4.2.7, 4.7.2 | public domain, Apache-2.0 | yes | sources of data, not replacements | not used |

Looked up via PyPI and the anaconda.org API on 2026-11-07. A GitHub keyword search for cross-database AMR mapping tools returned nothing useful. nf-core, bio.tools and recent papers were not searched, so "no close competitor" is unverified. NCBI's new `hierarchy_node` and `synonyms` catalogue columns may overlap with chamredb's purpose (not yet assessed).

## Findings

### Facts

- Shipped graph: 11,718 nodes, 33,632 edges (`chamredb/data/graph.json`); README still says 19,132 matches.
- Shipped inputs: CARD ARO 2021-04-23, NCBI 2021-03-01.1 (file headers); ResFinder snapshot date not recorded.
- Builder (Snakemake, conda, MMseqs2 pinned 13.45111, `--cores 1`) uses `latest` URLs and a `git clone` of ResFinder master, records no versions, has no tests.
- CARD latest download now includes `card.json`, `aro_index.tsv`, variant, knockout and overexpression models; the builder uses only the homolog model and the OBO.
- AMRFinderPlus database format is 4.2.0; ResFinder `phenotypes.txt` columns differ from the 2021 layout.
- hAMRonization input accession handling (`extract_database_and_id`) is tool-specific string parsing.
- Checked 2026-11-07 on downloaded current files (scratchpad): CARD release (`card.json` 4.0.2, 2026-08-11) has 6,059 homolog-model proteins against 2,979 shipped; the builder takes the OBO from ARO GitHub master and sequences from the CARD release, so the two can differ in release.
- NCBI current catalogue: 11,554 rows, 10,155 core AMR, 8,760 unique protein accessions among them (4,542 shipped); `db_version` 2026-10-29.1. Columns renamed (`genbank_strand_orientation` to `genbank_strand`, `blacklist` to `blacklisted_taxa`), new `synonyms`, `hierarchy_node`, `db_version`; the six columns the parsers use are present. `AMRProt.fa` (10,160 sequences) uses pipe-delimited headers starting with the accession; the parser has not been run on it.
- ResFinder db `VERSION` 2.6.0: `all.fsa` (3,212 sequences) duplicates the per-class files, so the builder's `cat *.fsa` would double every sequence (shipped 2021 file had 3,160). `phenotypes.txt` still has `Gene_accession no.` and `Phenotype`.
- CI is GitLab CI on Python 3.7; setup.py with unpinned dependencies; published as chAMReDb 0.2.1 on PyPI.

### Risks

| Risk | Likelihood | Impact | Mitigation |
| --- | --- | --- | --- |
| 2021 graph cannot be reproduced (unrecorded ResFinder date, MMseqs2 build) | medium | high | ADR 0003 step 1 first; record the mismatch as a finding |
| Upstream schema changes break parsers (NCBI 4.2.0 columns, ResFinder phenotypes, CARD models) | high | high | parser tests on current files; fail loudly on missing columns |
| MMseqs2 upgrade changes hits, mixing method and data effects | high | medium | upgrade separately, diff report |
| Accession or ID changes break hAMRonization lookups in the query path | medium | high | test against current hAMRonization 1.3.1 output; map removed IDs explicitly |
| ResFinder `cat *.fsa` doubles sequences via `all.fsa`, inflating hits and RBH ties | certain | high | use `all.fsa` only (or per-class files), assert unique sequence ids |
| CARD OBO (ARO master) and sequences (release tarball) from different releases | high | high | take both from one CARD release; record the version |
| MMseqs2 13.45111 unavailable on macOS arm64 | medium | low | run step 1 on Linux CI |
| New CARD models (variant, knockout) tempt scope creep | medium | medium | out of scope this round |

### Open questions

| Question | Owner | Blocking? |
| --- | --- | --- |
| Which ResFinder commit matches the 2021 graph? Try git history around 2021-10 | @antunderwood | no |
| Should removed IDs stay queryable (tombstones) after the refresh? | @antunderwood | no |
| Release naming and cadence for scheduled builds | @antunderwood | no |

## Scale and constraints

- Build runs on a GitHub Actions hosted runner, target up to about 2 hours; runtime to be measured in step 1 (single-core in 2021, not recorded). Executor: local.
- Consumers: mainly the maintainer and original contributors; no known external pipelines. Additive-schema rule still applies.

## Cost of being wrong

Wrong cross-database equivalences feed AMR interpretation by researchers choosing databases; errors are quiet. Oracle strength: reproduction of the 2021 graph, then reported diffs. A human scientific review (@antunderwood) is required for each data upgrade diff report.

## Pre-mortem

User confirmed all three below as the failures of most concern.

1. Refreshed graph silently loses or reassigns IDs that users have in saved hAMRonization outputs: add an ID-continuity report to each release.
2. Parser accepts a changed upstream format and produces plausible but wrong metadata: strict schema checks and row-count checks against upstream totals.
3. Nobody can tell which data a result came from: data release recorded per database and exposed in query output (additive field).

## Decisions

- ADR 0002: data build pipeline lives in this repo and stays Snakemake
- ADR 0003: prove reproduction of the 2021 graph before refreshing data or tools
- ADR 0004: distribute the graph as a versioned release asset; keep output schemas additive

## Out of scope

- New databases (MEGARes, ARG-ANNOT, PointFinder) and new CARD model types
- Changing the matching method (RBH plus one-way hits via MMseqs2)
- Nextflow port, bioconda packaging, breaking changes to CLI or output columns

## Module-design survey (ticket #3, 2026-11-09)

- **Duplicated graph logic.** `chamredb/functions/graph_functions.py` (`add_rbh_hits_to_graph`, `add_search_hits_to_graph`, `__add_connection`, `__load_metadata`) is copied into `builder/scripts/functions/graph_building.py`. The package copy is exercised by `chamredb/tests`; the builder copy produces the shipped graph. They can drift. Candidate: one shared function module (decide in #5, not changed here).
- **I/O mixed with science in the builder.** `hit_parsing`, `graph_building` and the parse scripts read and write fixed relative paths (`data/...`) inside the functions, so they cannot be unit-tested without a directory layout. The fixture build in `builder/tests` is the seam today. Candidate: pure functions taking dataframes (as `add_*_hits_to_graph` already do), paths owned by the scripts.
- **Parsers are scripts, not modules.** `parse_*.py` have logic at import time and no functions, so the strict-parser work in #4 needs a small refactor first.
- **Conventions owned in several places.** Node id `database:id`, identity as a fraction and coverage as alnlen over query end are created in `graph_building` and re-derived in the package's visualisation code. New module `builder/scripts/functions/graph_invariants.py` is now the one place that states what a valid graph is.
- **Graph loading.** `read_graph()` reads the packaged `graph.json` through networkx `node_link_graph`; this breaks on networkx 3.x because the default edge key changed from `links` to `edges` (see findings below). Relevant to #5 and #6.

## Findings from the invariant tests (ticket #3)

- **Coverage above 1.** 729 of 33,632 shipped edges (2.2%; 715 OWH, 14 RBH) have coverage above 1 (max 1.198): coverage is alignment length over query end, which exceeds 1 for gapped alignments. The spec assumption 0 < coverage <= 1 is false. **Decision (2026-11-09, @antunderwood): keep the definition, document it in GLOSSARY.md, relax the invariant to coverage > 0.** Predicate 3 in spec #1 is amended accordingly.
- **ResFinder duplicate ids.** The committed `resfinder.nucl.fasta` (cat of per-class `*.fsa`) has 5 duplicate ids (genes in two class files); the parsed protein FASTA is unique. Predicate 4 holds for the parsed output, not for the raw concatenation.
- **Predicate 8 and networkx.** `chamredb/tests` pass (7) with networkx 2.8.8, pandas 1.5.3 and Python 3.9, and fail (2 failed, 3 errors) with networkx 3.6.1, because `node_link_data` and `node_link_graph` now default to the `edges` key. The package needs a pinned bound or an explicit `edges="links"` argument (additive-schema rule: keep `links` in graph.json).
