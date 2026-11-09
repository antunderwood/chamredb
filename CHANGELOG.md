# Changelog

## 0.3.0 (2026-11-09)

Data refresh and modernisation (spec: GitHub issue #1).

### Data

- New graph built from CARD 4.0.2 (sequences and ontology from the same release), NCBI AMRFinderPlus database 4.2 (2026-10-29.1) and ResFinder db 2.6.0 (commit `eecf0aa`) with MMseqs2 18.8cc5c: 15,859 nodes and 51,971 edges (was 11,718 and 33,632, built 2021-10-20).
- Each upgrade was reviewed as a diff against the previous step (MMseqs2, then CARD, NCBI, ResFinder); reports are in the pull requests #28 to #31. Summary of the changes seen: all reciprocal best hits unchanged by the MMseqs2 upgrade, 2.4% of one-way hits shifted; CARD +1,491 nodes, 108 names changed, 17 ids absent; NCBI +2,633 nodes, 23 accessions renamed, 29 NCBI and 2 ResFinder ids absent; ResFinder step 29 ResFinder and 13 CARD ids absent, 51 edges changed type.
- Ids removed by a refresh are not kept as tombstones (ADR 0005); each release lists them in its diff and ID-continuity report.

### Added

- `chamredb query --graph PATH`: use another graph file (offline use).
- `data release:` line in single-id output and an extra last TSV column `data_release` (existing columns unchanged).
- `builder/`: the data build (Snakemake, formerly the separate `chamredb-builder`) with strict parsers, a data manifest, graph invariants, a graph diff and ID-continuity report, a monthly release workflow, and a fixture-based test of the whole build.

### Changed

- Python 3.10 or newer (was 3.7); dependencies bounded (`networkx>=3.4`, `pandas>=2.0`, `rich>=13`); `pyproject.toml` replaces `setup.py`; unused `pronto`, `seaborn` and `numpy` runtime dependencies removed.
- `chamredb query` works with current networkx (explicit `edges="links"`; `graph.json` keeps `links`).
- ResFinder: the phenotype of a gene now comes from the allele whose sequence is used (9 of 2,601 genes changed).
- The unit tests read a fixed subset of the 2021 graph, so data refreshes cannot change their expected outputs.

### Known

- The assessment in the README ("Assessing the graph") describes an earlier graph and was not regenerated.
- Edge `coverage` is alignment length over the query end position and can exceed 1 (2.6% of edges in the 2026-11-09 graph, 2.2% in the 2021 graph).
- Metadata used only to name ids that are missing from the graph is still the 2021 version.
