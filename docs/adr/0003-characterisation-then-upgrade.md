# 0003. Prove reproduction of the 2021 graph before refreshing data or tools

- Status: accepted (2026-11-10; premise corrected 2026-11-09)
- Date: 2026-11-07
- Deciders: @antunderwood
- Kind: scientific

## Context

The shipped `graph.json` (11,718 nodes, 33,632 edges) was built from CARD ARO 2021-04-23, NCBI Reference Gene Catalog 2021-03-01.1, a ResFinder snapshot of unknown date, and MMseqs2 13.45111. Upstream data and MMseqs2 (18.8cc5c on bioconda) have both moved. Changing everything at once makes any difference in the new graph unattributable.

## Options considered

1. Existing unit fixtures only: fast, cannot show that the pipeline reproduces the shipped graph
2. Two-step characterisation: reproduce 2021 first, then change one thing at a time

## Decision

Step 1: rebuild from inputs reconstructed from archived releases nearest the build date 2021-10-20 (CARD 3.1.4 or nearest, NCBI 2021-09-30.1 or nearest, ResFinder commit by date), with MMseqs2 13.45111 (Linux CI, since the bioconda build may not exist for macOS arm64) and require the node set and edge set (identical node set and edge set and edge type; identity and coverage within 0.001) to match the shipped `graph.json`. Step 2: upgrade MMseqs2, then each database, as separate changes, each with a written diff report (edges gained, lost, changed type).

## Consequences

The shipped graph was built in builder CI on 2021-10-20 from `latest` downloads, not from the April 2021 FASTAs committed in `chamredb/data` (1,599 CARD and 168 NCBI nodes have no committed protein; found in #9, see #2). The reconstructed 2021 inputs become a test fixture or release asset. If step 1 cannot be matched, the mismatch is itself a finding, recorded before any upgrade. Each data release records source URL, release and retrieval date per database.

## Result

Step 1 matched (CI run 42593395624, native amd64, MMseqs2 13.45111, pull request 20): 11,718 nodes and 33,632 edges, an identical node set, edge set and edge types, and no edge with identity or coverage off by more than 0.001. ResFinder commit 7562716 matched on the first try and the build took 41.6 minutes. The only differences were 5 CARD node attributes (phenotype or additional_phenotype missing or moved on ARO:3003907, 3004146, 3005087, 3005088 and 3005121), consistent with a slightly different ontology or CARD release than the 2021-10-20 build; they were reported and did not fail the comparison. Step 2 then ran as four reviewed changes (MMseqs2 18.8cc5c, CARD 4.0.2, NCBI 2026-10-29.1, ResFinder 2.6.0), summarised in `docs/release-2026-11-09.md`.
