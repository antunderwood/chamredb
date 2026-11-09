#!/usr/bin/env bash
# Rebuild the graph from the archived 2021 inputs and compare it with the shipped graph (ADR 0003 step 1).
# Usage: run_reproduction.sh <workdir>  -> <workdir>/graph.json, data_manifest.json, report.md, timings.txt, versions.txt
# Needs the builder environment (builder/environment.yml, MMseqs2 13.45111 on Linux amd64) and network.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
builder=$(cd "$here/.." && pwd)
shipped=$(cd "$builder/../chamredb/data" && pwd)/graph.json
work=${1:?usage: run_reproduction.sh <workdir>}
case "$work" in ""|/|.|..|~) echo "refusing to use $work as workdir" >&2; exit 2;; esac
rm -rf "$work"; mkdir -p "$work"
cp "$builder/Snakefile" "$work/"; cp -r "$builder/scripts" "$work/scripts"
cd "$work"
export MMSEQS_NUM_THREADS=1   # one thread keeps MMseqs2 output order deterministic
t0=$(date +%s)
python "$here/fetch_inputs.py" data
python scripts/write_manifest.py --inputs "$here/inputs.json" --data data --out data_manifest.json
t1=$(date +%s)
targets=""
for db in card ncbi resfinder; do targets="$targets data/db_metadata/$db.metadata.json"; done
for a in card ncbi resfinder; do
  for b in card ncbi resfinder; do
    [ "$a" = "$b" ] && continue
    targets="$targets data/mmseqs2_results/${a}_vs_${b}.rbh.tsv data/mmseqs2_results/${a}_vs_${b}.search.tsv"
  done
done
snakemake $targets --cores 1 --rerun-triggers mtime --config bounds="$here/bounds.json" --quiet rules
t2=$(date +%s)
python scripts/build_graph.py > build_graph.log
t3=$(date +%s)
{ echo "fetch_s $((t1-t0))"; echo "parse_and_mmseqs2_s $((t2-t1))"; echo "build_graph_s $((t3-t2))"; echo "total_s $((t3-t0))"; } > timings.txt
{
  echo "mmseqs2 $(mmseqs version)"; echo "snakemake $(snakemake --version)"
  python - <<'PY'
import sys, pandas, networkx, Bio, pronto
from importlib.metadata import version
print("python", sys.version.split()[0])
for m in (pandas, networkx, Bio, pronto):
    print(m.__name__, m.__version__)
print("orfipy", version("orfipy"))
PY
} > versions.txt
python "$here/compare_graphs.py" "$shipped" graph.json report.md || status=$?
cat timings.txt
exit "${status:-0}"
