# 0005. No tombstones, dated graph releases checked monthly, GitLab as a mirror

- Status: accepted
- Date: 2026-11-09
- Deciders: @antunderwood
- Kind: engineering

## Context

Refreshing the data (spec #1) will remove or rename some IDs. Users keep saved hAMRonization outputs that reference IDs. The graph will be built and released by scheduled GitHub Actions (ADR 0004). The repository also has a GitLab `origin` and a `.gitlab-ci.yml` pinned to Python 3.7, which cannot run the package any more (Python 3.10 or newer, ADR 0004).

## Options considered

1. Tombstones: ship a removed-ids table and make `chamredb query` explain removed ids
2. Report only: removed IDs leave the graph; each release carries an ID-continuity report (`builder/scripts/diff_graphs.py`)
3. Release naming: dated tags, or semantic version tags; build monthly, manually, or on upstream change
4. GitLab: retire the CI file and keep a mirror, retire GitLab entirely, or leave as is

## Decision

- No tombstones. Removed IDs are not queryable. Every release reports ID continuity (kept, reassigned, absent with same-name successors) so users can see what moved. Revisit if users report lost saved outputs.
- Release tags are dated, `graph-YYYY.MM.DD`. A scheduled workflow checks monthly and builds and publishes only when an upstream release (CARD, NCBI or ResFinder) has changed; it can also be run by hand.
- GitHub is the source of truth. `.gitlab-ci.yml` is deleted. The GitLab `origin` stays as a mirror; pushing to it is the maintainer's action.

## Consequences

A user with a saved output that names a removed ID gets "not found" and must read the release's continuity report. Schemas stay unchanged (ADR 0004). The release workflow (#8) needs an upstream-change check. A tombstone table can be added later without a breaking change.
