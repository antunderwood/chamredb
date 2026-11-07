# GLOSSARY

Agents and humans use these words with exactly these meanings. Resolve ambiguity here first; add terms as they are agreed.

## Project

- **Project goal:** given an antimicrobial resistance gene (ARG) in one database, find its equivalents in the other databases and aggregate their metadata.
- **Organism / assembly / annotation:** none (see .bio-skills.yaml).

## Terms

| Term | Meaning in this project | Not to be confused with |
| --- | --- | --- |
| Source database | One of CARD, NCBI (AMRFinderPlus Reference Gene Catalog) or ResFinder | Target database |
| RBH | Reciprocal best hit: each protein is the other's best MMseqs2 hit | One-way hit |
| OWH | One-way hit: best hit from source to target with no reciprocal hit, kept if coverage and identity thresholds are met | RBH |
| Coverage | Alignment length divided by query length | Depth |
| Identity | MMseqs2 fraction identical (`fident`) over the alignment | Name similarity |
| Name similarity | Normalised Levenshtein similarity of cleaned database names | Sequence identity |
| Node | One protein from one database, id `database:id` | Gene family |
| Edge | Directed match between nodes, with type, coverage and identity | |
| Graph | The networkx DiGraph of all nodes and edges, shipped as `chamredb/data/graph.json` | |
| Snapshot | The retrieved upstream release of each database that a graph was built from | Package version |
