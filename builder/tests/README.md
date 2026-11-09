# Builder characterisation harness

`make test` builds the graph from a tiny fixture and compares it with `expected/fixture_graph.json`.
It needs MMseqs2 18.8cc5c (Linux or macOS; the 2021 reproduction still uses 13.45111, Linux amd64 only). `make test-docker` runs the same in a
container; it takes about 5 minutes under emulation on arm64, so the 60 s target in `.bio-skills.yaml`
is only checkable on native Linux CI (the workflow prints the measured time).

## What the fixture is

- `fixtures/data`: about 20 proteins per database taken from the public files committed in `chamredb/data`,
  chosen by `make_fixtures.py`. Two caveats: the CARD raw FASTA headers are reconstructed (only parsed CARD
  FASTAs are committed, so the accession is synthetic and the parsers read only the ARO field), and decoys
  are chosen using the shipped `graph.json` (proteins with no edges there).
- It is **not** a reproduction of the shipped 2021 `graph.json`, which was built from later downloads than
  the committed files (see ticket #2).
- Designed cases: mutual best hits across all three databases (dfrB1, 2, 3, 6, tet(37)); a name mismatch
  (CARD dfrB7 is dfrB8 elsewhere); one-way hits only (dfrB9 to ResFinder, ResFinder dfrB7); CARD-NCBI only
  pairs (mecI, vanU-G); a partial-coverage pair (CARD ARO:3002718 is the NCBI protein plus one residue, so
  coverage is 214/215, derived from the sequences in the test); NCBI non-core rows and shared names;
  ResFinder multi-allele genes.
- Not designed: card ARO:3001205 and ncbi WP_004201167.1 (ble) turned out to be identical and form an RBH
  pair, and card ARO:3005036 hits that protein one way. The ble pair exists only because of which accession
  survives de-duplication, so treat it as characterisation.
- The one-way-hit target sets (`--max-accept 3` keeps three) are characterisation too: near ties
  (0.853, 0.851, 0.841) are tool-version sensitive.

## What the oracle is

`expected/` holds the frozen output, the tool versions (`VERSIONS.txt`) and fixture checksums
(`SHA256SUMS`). Frozen output proves "unchanged", not "correct". Update it only with a written reason
why the old values were wrong, never to make a test pass. The designed cases and the intermediate-file
checks (unique ids, one NCBI protein per name, nodes are parsed proteins) are the hand-checkable part.

## Mutation check (the harness must go red)

Run in a scratch copy of a finished build, rebuild the affected stage, and run the tests with
`FIXTURE_BUILD_DIR` pointing at the copy. Last run 2026-11-07:

| Mutation | Fails |
| --- | --- |
| `graph_building.py`: coverage replaced by `1.0` | `test_matches_frozen_graph`, `test_coverage_is_alignment_length_over_query_end` |
| `Snakefile`: `--max-accept 3` to `1` | `test_matches_frozen_graph` |
| `parse_ncbi_proteins.py`: remove `.drop_duplicates(subset=['allele'])` | the two NCBI parser tests |
| `hit_parsing.py`: per-query sort order reversed | **survives**: order is not observable on this fixture |

## Not covered (follow-up tickets)

- General invariants for any graph, including 0 < coverage <= 1 (ticket #3).
- Reverse-strand ResFinder genes (ORF extraction by orfipy runs only on the fixture genes), and ResFinder
  alleles with different phenotypes: sequence and metadata can come from different alleles.
- Renamed NCBI columns and the ResFinder `all.fsa` duplication (tickets #4, #13, #14): the unmodified
  2021 parsers cannot pass them, so they are not in this fixture.
- The download rules and `rule all`: the harness requests explicit targets with the fixture files in place
  of the downloads.

## Determinism of the CARD metadata (fixed in ticket #18)

`parse_card_metadata.py` used to take `PMID` lists, the `is_a` parent list and the term order for phenotype
inheritance from sets, so output changed with `PYTHONHASHSEED` (3 to 8 CARD terms differed between seeds on the
2021 ontology). It now sorts `PMID`, sorts parents and visits terms in sorted id order.
`test_card_metadata_determinism.py` checks byte identity over 10 seeds on a hand-built ontology; the fixture
graph compares `PMID` exactly (expected file stores them sorted). Phenotype inheritance still follows the
first parent only and depends on visit order by design of the 2021 code, now a fixed order.

## Notes

- `FIXTURE_BUILD_DIR=<dir>` reuses a finished build. It has no staleness check, so a build from before
  a code change will pass.
- `builder/environment.yml` was pinned for the harness (Python 3.9, Snakemake 7.32.4, pandas 1.5.3,
  networkx 2.8.8, Biopython 1.81, pronto 2.4.3, MMseqs2 18.8cc5c, pytest). The `defaults` channel was dropped,
  `pip` dependencies moved to conda, and `git` and `orfipy` are not pinned in the file (orfipy 0.0.4 is
  recorded in `VERSIONS.txt` and asserted). Python 3.9 applies to the builder only; the package target is
  Python 3.10 or newer (ADR 0004).
