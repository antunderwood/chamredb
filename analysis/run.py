"""Write the README statistics and the two low-name-similarity tables from a graph file (spec #34, ticket #38).

    python -m analysis.run --graph PATH [--out docs]

Writes <out>/statistics.json and <out>/analysis/low_name_similarity{,_without_alleles}.csv. Deterministic: no timestamps,
sorted keys, rows in graph order, so two runs on the same graph give byte-identical files.
"""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from analysis import matches


def statistics(graph, result, graph_name, graph_sha256):
    database = lambda node_id: node_id.split(":", 1)[0]
    rbh_pairs = {frozenset((u, v)) for u, v, d in graph.edges(data=True) if d["type"] == "RBH"}
    by_databases = Counter("+".join(sorted({database(n) for n in pair})) for pair in rbh_pairs)

    def block(fit, predicted, high, low):
        return {"slope": fit[0], "intercept": fit[1], "predicted_similarity_at_threshold": predicted,
                "matches_at_or_above_threshold": len(high), "low_name_similarity_rows": len(low)}

    return {
        "graph": {
            "file": graph_name, "sha256": graph_sha256, "nodes": graph.number_of_nodes(), "edges": graph.number_of_edges(),
            "nodes_by_database": dict(sorted(Counter(database(n) for n in graph.nodes).items())),
            "edges_by_type": dict(sorted(Counter(d["type"] for _, _, d in graph.edges(data=True)).items())),
            "rbh_pairs": len(rbh_pairs), "rbh_pairs_by_databases": dict(sorted(by_databases.items())),
        },
        "analysis": {
            "identity_threshold": matches.IDENTITY_THRESHOLD, "matches_analysed": len(result["distance"]),
            "matches_skipped_blaPAO": result["skipped_blaPAO"],
            "full_names": block(result["linear_fit"], result["predicted_similarity"],
                                result["high_identity"], result["low_name_similarity"]),
            "without_alleles": block(result["linear_fit_without_alleles"],
                                     result["predicted_similarity_without_alleles"],
                                     result["high_identity_without_alleles"], result["low_name_similarity_without_alleles"]),
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--graph", required=True, help="graph JSON file (node-link, key 'links')")
    parser.add_argument("--out", default="docs", help="output directory (default docs)")
    args = parser.parse_args(argv)
    graph_path = Path(args.graph)
    graph = matches.read_graph(graph_path)
    result = matches.assess(graph)
    out = Path(args.out)
    (out / "analysis").mkdir(parents=True, exist_ok=True)
    sha = hashlib.sha256(graph_path.read_bytes()).hexdigest()
    (out / "statistics.json").write_text(
        json.dumps(statistics(graph, result, graph_path.name, sha), sort_keys=True, indent=2) + "\n")
    result["low_name_similarity"].to_csv(out / "analysis" / "low_name_similarity.csv", index=False)
    result["low_name_similarity_without_alleles"].to_csv(out / "analysis" / "low_name_similarity_without_alleles.csv", index=False)
    print(f"wrote {out}/statistics.json and two tables for {graph.number_of_nodes()} nodes")


if __name__ == "__main__":
    main()
