# 0001. Data policy is open

- Status: accepted
- Date: 2026-11-07
- Deciders: @antunderwood
- Kind: data-policy

## Context

All inputs are public, non-human reference databases: CARD (ARO ontology and sequences), NCBI AMRFinderPlus Reference Gene Catalog and sequences, and ResFinder (resfinder_db). They are already committed under `chamredb/data/`. No patient, sample or controlled-access data are in scope.

## Options considered

1. restricted: safest default, but blocks reading and committing the public reference data that this project exists to process
2. open: agents may read, print and commit the public databases, with provenance (source URL, release, retrieval date) recorded

## Decision

`data_policy: open`. Every data file added or refreshed must record its source URL and release or version.

## Consequences

Agents may inspect and commit database files. If any non-public or human-derived data are ever added (for example clinical isolate hAMRonization outputs), this ADR must be superseded and the policy set to `restricted` for those paths.
