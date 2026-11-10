# 0004. Distribute the graph as a versioned release asset; keep output schemas additive

- Status: accepted (2026-11-10), amended to match what shipped in 0.3.0
- Date: 2026-11-07
- Deciders: @antunderwood
- Kind: engineering

## Context

The 9 MB `graph.json` (13.7 MB for the 2026-11-09 release) is shipped inside the Python package, so every data refresh needs a package release. Existing users rely on the TSV output columns, the hAMRonization JSON input and the CLI flags.

## Options considered

1. Keep the graph in the package: simplest for users, couples data and code releases
2. Versioned GitHub release asset, downloaded and cached on first use, built by a scheduled GitHub Actions workflow

## Decision

Option 2 in part. Each graph is built by a GitHub Actions workflow and published as a versioned release asset (`graph-YYYY.MM.DD`) with its data manifest, diff report and database metadata. The package still bundles the graph of one release and pins it in `chamredb/data/graph_release.json` (release tag, sha256 of the graph and of each metadata file); `chamredb query` prints the data release, takes `--graph PATH` for any other release, and gives no names from a custom graph because its metadata is unknown. There is no download on first use: the bundled graph keeps the package usable offline, and a data refresh still needs a package release. Output TSV columns, JSON input and CLI flags stay stable and change only additively (the `data_release` column and the `data release:` line). Supported Python is 3.10 or newer, distributed on PyPI only.

## Consequences

Needs the pin to be checked on load (it is, by sha256) and tests that work without network (they do). A download and cache path was not built and can be added later without changing the output schemas. A breaking change requires a new ADR and a major version.
