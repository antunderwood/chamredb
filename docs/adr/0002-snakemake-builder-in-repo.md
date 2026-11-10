# 0002. Data build pipeline lives in this repo and stays Snakemake

- Status: accepted (2026-11-10)
- Date: 2026-11-07
- Deciders: @antunderwood
- Kind: engineering

## Context

Since commit 6c2645e (2021-10-20) the graph build lived in the separate `chamredb-builder` repo (Snakemake, conda, GitLab CI producing a `graph.json` artefact). That split left this repo with no way to rebuild or refresh its own data, and no record of which upstream releases produced the shipped graph.

## Options considered

1. Keep two repos: data and code drift apart, as they have for 5 years
2. Merge builder into this repo, keep Snakemake: smallest change, existing rules keep working
3. Merge and port to Nextflow DSL2: matches personal preference, but is a rewrite with no scientific gain

## Decision

Merge `chamredb-builder` into this repo with its history (under `builder/`), keep Snakemake, and add `snakemake-pipeline` to `project_types`. Local execution is the default executor.

## Consequences

One repo owns data, build and query. Pipeline tests and a test profile are now required. The old builder repo is to be archived after the merge (a manual step for the maintainer, not done yet).
