#!/usr/bin/env bash
# Build only the three database metadata JSON files (no MMseqs2) from a pinned inputs file, in the active
# builder environment. Usage: build_metadata.sh <workdir> [inputs.json]   (default inputs: builder/refresh/inputs_resfinder.json)
# Writes <workdir>/data/db_metadata/{card,ncbi,resfinder}.metadata.json and <workdir>/data_manifest.json.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
builder=$(cd "$here/.." && pwd)
work=${1:?usage: build_metadata.sh <workdir> [inputs.json]}
case "$work" in ""|/|.|..|~) echo "refusing to use $work as workdir" >&2; exit 2;; esac
inputs=${2:-$builder/refresh/inputs_resfinder.json}
rm -rf "$work"; mkdir -p "$work"
cp "$builder/Snakefile" "$work/"; cp -r "$builder/scripts" "$work/scripts"
cd "$work"
python "$builder/reproduce_2021/fetch_inputs.py" data "$inputs"
python scripts/write_manifest.py --inputs "$inputs" --data data --out data_manifest.json
snakemake data/db_metadata/card.metadata.json data/db_metadata/ncbi.metadata.json data/db_metadata/resfinder.metadata.json \
  --cores 1 --rerun-triggers mtime --quiet rules
ls -la data/db_metadata/*.metadata.json
