#!/usr/bin/env bash
# Rebuild the graph from the archived 2021 inputs and compare it with the shipped graph (ADR 0003 step 1).
# Usage: run_reproduction.sh <workdir>  -> <workdir>/graph.json, data_manifest.json, report.md, invariants.txt, timings.txt, versions.txt
# REPORT_ONLY=1 writes an informational diff report (for tool upgrades) instead of the pass/fail comparison.
# Run in builder/environment-mmseqs13.yml (MMseqs2 13.45111, Linux amd64) to reproduce the 2021 graph; any other
# environment (REPORT_ONLY=1) gives a diff report and network.
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
inputs=${INPUTS:-$here/inputs.json}   # INPUTS=<file> selects another set of pinned inputs (e.g. builder/refresh)
python "$here/fetch_inputs.py" data "$inputs"
python scripts/write_manifest.py --inputs "$inputs" --data data --out data_manifest.json
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
python - <<'PY' > invariants.txt || status=$?
import json, sys
sys.path.insert(0, "scripts/functions")
import graph_invariants
v = graph_invariants.check_graph(json.load(open("graph.json")))
print("graph invariants (predicates 1 to 3):", "hold" if not v else {k: len(x) for k, x in v.items()})
sys.exit(1 if v else 0)
PY
cat invariants.txt
if [ -n "${REPORT_ONLY:-}" ]; then
  python scripts/diff_graphs.py "$shipped" graph.json --out report.md
else
  python scripts/functions/graph_compare.py "$shipped" graph.json report.md || status=$?
fi
cat timings.txt
exit "${status:-0}"
