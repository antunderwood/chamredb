"""Write the README statistics and the two low-name-similarity tables from a graph file (spec #34, ticket #38).

    python -m analysis.run --graph PATH [--out docs]

Writes <out>/statistics.json, <out>/analysis/low_name_similarity{,_without_alleles}.csv and the figures under <out>/images. Deterministic: no timestamps,
sorted keys, rows in graph order, so two runs on the same graph give byte-identical files.
"""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from analysis import groups, matches


def group_statistics(graph):
    """The `groups` block of statistics.json: equivalence groups (proteins joined by RBH) by database profile."""
    found = groups.equivalence_groups(graph)
    p = groups.profiles(found)
    return {"total": len(found), "by_databases": p["groups"], "proteins_by_databases": p["proteins"],
            "size_distribution": {str(k): v for k, v in p["size_distribution"].items()},
            "largest_group_size": p["largest_group_size"], "largest_group_by_database": p["largest_group_by_database"]}


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
        "groups": group_statistics(graph),
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
    from analysis import figures  # needs matplotlib (optional extra `analysis`)
    import matplotlib.pyplot as plt
    figure = figures.name_similarity_vs_identity(result, statistics(graph, result, graph_path.name, sha), out / "images")
    plt.close(figure)
    plt.close(figures.overlap_upset(group_statistics(graph), out / "images"))
    print(f"wrote {out}/statistics.json, two tables and the figures for {graph.number_of_nodes()} nodes")


if __name__ == "__main__":
    main()
