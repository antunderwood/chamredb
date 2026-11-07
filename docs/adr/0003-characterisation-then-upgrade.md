# 0003. Prove reproduction of the 2021 graph before refreshing data or tools

- Status: proposed
- Date: 2026-11-07
- Deciders: @antunderwood
- Kind: scientific

## Context

The shipped `graph.json` (11,718 nodes, 33,632 edges) was built from CARD ARO 2021-04-23, NCBI Reference Gene Catalog 2021-03-01.1, a ResFinder snapshot of unknown date, and MMseqs2 13.45111. Upstream data and MMseqs2 (17.b804f on bioconda) have both moved. Changing everything at once makes any difference in the new graph unattributable.

## Options considered

1. Existing unit fixtures only: fast, cannot show that the pipeline reproduces the shipped graph
2. Two-step characterisation: reproduce 2021 first, then change one thing at a time

## Decision

Step 1: rebuild from the committed 2021 inputs with MMseqs2 13.45111 (Linux CI, since the bioconda build may not exist for macOS arm64) and require the node set and edge set (identical node set and edge set and edge type; identity and coverage within 0.001) to match the shipped `graph.json`. Step 2: upgrade MMseqs2, then each database, as separate changes, each with a written diff report (edges gained, lost, changed type).

## Consequences

The frozen 2021 inputs become a test fixture or release asset. If step 1 cannot be matched, the mismatch is itself a finding, recorded before any upgrade. Each data release records source URL, release and retrieval date per database.
