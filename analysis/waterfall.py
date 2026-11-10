"""Step summary for the release waterfall (spec #34, ticket #51).

Each transition compares two graphs of consecutive data steps with builder/scripts/functions/graph_diff (nodes and edges gained and
lost, edges whose type changed, ids absent and renamed). The summary is committed as docs/waterfall_steps.json, so the figure does not
need the graph files; each step records the sha256 of its graph and where it came from.

    python -m analysis.waterfall --out docs/waterfall_steps.json \\
        --step "2021 graph" PATH SOURCE --step "MMseqs2 18" PATH SOURCE ...
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "builder" / "scripts" / "functions"))
import graph_diff  # noqa: E402


def transition(old, new):
    d = graph_diff.diff(old, new)
    c = d["continuity"]
    return {
        "nodes_before": d["shipped_nodes"], "nodes_after": d["rebuilt_nodes"],
        "nodes_gained": len(d["nodes_gained"]), "nodes_lost": len(d["nodes_lost"]),
        "edges_before": d["shipped_edges"], "edges_after": d["rebuilt_edges"],
        "edges_gained": len(d["edges_gained"]), "edges_lost": len(d["edges_lost"]),
        "edges_type_changed": len(d["type_changed"]),
        "ids_absent": len(c["absent"]),
        "ids_absent_by_database": {db: n["absent"] for db, n in sorted(c["by_database"].items()) if n["absent"]},
        "ids_renamed": len(c["reassigned"]),
    }


def build(steps):
    """steps: list of (label, graph path, source description); the first is the starting graph."""
    graphs = [json.loads(Path(path).read_text()) for _, path, _ in steps]
    return {
        "steps": [{"label": label, "source": source, "graph_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                   "nodes": len(g["nodes"]), "edges": len(g["links"])} for (label, path, source), g in zip(steps, graphs)],
        "transitions": [{"from": steps[i][0], "to": steps[i + 1][0], **transition(graphs[i], graphs[i + 1])} for i in range(len(steps) - 1)],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True)
    parser.add_argument("--step", nargs=3, action="append", metavar=("LABEL", "GRAPH", "SOURCE"), required=True)
    args = parser.parse_args(argv)
    Path(args.out).write_text(json.dumps(build([tuple(s) for s in args.step]), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
