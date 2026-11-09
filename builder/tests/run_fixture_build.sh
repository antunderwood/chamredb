#!/usr/bin/env bash
# Build the fixture graph in a scratch directory (nothing is written inside the repository).
# Usage: run_fixture_build.sh <workdir>   -> <workdir>/graph.json and <workdir>/versions.txt
set -euo pipefail
builder=$(cd "$(dirname "$0")/.." && pwd)
work=${1:?usage: run_fixture_build.sh <workdir>}
case "$work" in ""|/|.|..|~) echo "refusing to use $work as workdir" >&2; exit 2;; esac
rm -rf "$work"
mkdir -p "$work"
cp "$builder/Snakefile" "$work/"
cp -r "$builder/scripts" "$work/scripts"
cp -r "$builder/tests/fixtures/data" "$work/data"
cd "$work"

export MMSEQS_NUM_THREADS=1   # one thread keeps MMseqs2 output order deterministic

targets=""
for db in card ncbi resfinder; do targets="$targets data/db_metadata/$db.metadata.json"; done
for a in card ncbi resfinder; do
  for b in card ncbi resfinder; do
    [ "$a" = "$b" ] && continue
    targets="$targets data/mmseqs2_results/${a}_vs_${b}.rbh.tsv data/mmseqs2_results/${a}_vs_${b}.search.tsv"
  done
done
# the download rules' outputs are the fixture files, so only parsing, MMseqs2 and graph rules run
snakemake $targets --cores 1 --rerun-triggers mtime --quiet rules
python scripts/build_graph.py > build_graph.log

{
  echo "mmseqs2 $(mmseqs version)"
  echo "snakemake $(snakemake --version)"
  python - <<'PY'
import sys, pandas, networkx, Bio, pronto
from importlib.metadata import version
print("python", sys.version.split()[0])
for m in (pandas, networkx, Bio, pronto):
    print(m.__name__, m.__version__)
print("orfipy", version("orfipy"))
PY
} > versions.txt
