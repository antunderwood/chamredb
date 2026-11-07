# 0004. Distribute the graph as a versioned release asset; keep output schemas additive

- Status: proposed
- Date: 2026-11-07
- Deciders: @antunderwood
- Kind: engineering

## Context

The 9 MB `graph.json` is shipped inside the Python package, so every data refresh needs a package release. Existing users rely on the TSV output columns, the hAMRonization JSON input and the CLI flags.

## Options considered

1. Keep the graph in the package: simplest for users, couples data and code releases
2. Versioned GitHub release asset, downloaded and cached on first use, built by a scheduled GitHub Actions workflow

## Decision

Option 2. The package pins a default graph release and records its data versions. Output TSV columns, JSON input and CLI flags stay stable and change only additively (for example a data release field). Supported Python is 3.10 or newer, distributed on PyPI only.

## Consequences

Needs a download, checksum and offline-cache path, and a test that works without network. A breaking change requires a new ADR and a major version.
