# 0003. Prove reproduction of the 2021 graph before refreshing data or tools

- Status: proposed (premise corrected 2026-11-09)
- Date: 2026-11-07
- Deciders: @antunderwood
- Kind: scientific

## Context

The shipped `graph.json` (11,718 nodes, 33,632 edges) was built from CARD ARO 2021-04-23, NCBI Reference Gene Catalog 2021-03-01.1, a ResFinder snapshot of unknown date, and MMseqs2 13.45111. Upstream data and MMseqs2 (17.b804f on bioconda) have both moved. Changing everything at once makes any difference in the new graph unattributable.

## Options considered

1. Existing unit fixtures only: fast, cannot show that the pipeline reproduces the shipped graph
2. Two-step characterisation: reproduce 2021 first, then change one thing at a time

## Decision

Step 1: rebuild from inputs reconstructed from archived releases nearest the build date 2021-10-20 (CARD 3.1.4 or nearest, NCBI 2021-09-30.1 or nearest, ResFinder commit by date), with MMseqs2 13.45111 (Linux CI, since the bioconda build may not exist for macOS arm64) and require the node set and edge set (identical node set and edge set and edge type; identity and coverage within 0.001) to match the shipped `graph.json`. Step 2: upgrade MMseqs2, then each database, as separate changes, each with a written diff report (edges gained, lost, changed type).

## Consequences

The shipped graph was built in builder CI on 2021-10-20 from `latest` downloads, not from the April 2021 FASTAs committed in `chamredb/data` (1,599 CARD and 168 NCBI nodes have no committed protein; found in #9, see #2). The reconstructed 2021 inputs become a test fixture or release asset. If step 1 cannot be matched, the mismatch is itself a finding, recorded before any upgrade. Each data release records source URL, release and retrieval date per database.
