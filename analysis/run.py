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


def pair_statistics(graph):
    """For each ordered pair of databases (source, target): the source proteins in the graph, those with a reciprocal best hit to the
    target, those with only a one-way hit, and those covered by either (the `pairs` block of statistics.json)."""
    database = lambda node_id: node_id.split(":", 1)[0]
    names = sorted({database(n) for n in graph.nodes})
    counts = {f"{a}>{b}": {"source_nodes": 0, "with_rbh": 0, "with_search_only": 0, "covered": 0} for a in names for b in names if a != b}
    for node in graph.nodes:
        a = database(node)
        kinds = {}
        for _, target, data in graph.out_edges(node, data=True):
            kinds.setdefault(database(target), set()).add(data["type"])
        for b in names:
            if b == a:
                continue
            row = counts[f"{a}>{b}"]
            row["source_nodes"] += 1
            has = kinds.get(b, set())
            row["with_rbh"] += "RBH" in has
            row["with_search_only"] += bool(has) and "RBH" not in has
            row["covered"] += bool(has)
    return counts


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
        "pairs": pair_statistics(graph),
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
    parser.add_argument("--manifest", help="data_manifest.json of the release; with it the pipeline diagram is drawn")
    args = parser.parse_args(argv)
    graph_path = Path(args.graph)
    graph = matches.read_graph(graph_path)
    result = matches.assess(graph)
    out = Path(args.out)
    (out / "analysis").mkdir(parents=True, exist_ok=True)
    sha = hashlib.sha256(graph_path.read_bytes()).hexdigest()
    stats = statistics(graph, result, graph_path.name, sha)
    (out / "statistics.json").write_text(json.dumps(stats, sort_keys=True, indent=2) + "\n")
    import platform
    from importlib.metadata import version
    (out / "analysis" / "environment.txt").write_text(
        f"python {platform.python_version()}\n" + "".join(f"{m} {version(m)}\n" for m in ("matplotlib", "networkx", "numpy", "pandas", "rich", "strsimpy")) +
        f"graph_sha256 {sha}\n")
    result["low_name_similarity"].to_csv(out / "analysis" / "low_name_similarity.csv", index=False)
    result["low_name_similarity_without_alleles"].to_csv(out / "analysis" / "low_name_similarity_without_alleles.csv", index=False)
    from analysis import families, figures  # needs matplotlib (optional extra `analysis`)
    import matplotlib.pyplot as plt
    images = out / "images"
    plt.close(figures.name_similarity_vs_identity(result, stats, images))
    plt.close(figures.name_distributions(result, stats, images))
    plt.close(figures.overlap_upset(group_statistics(graph), images))
    if all(n in graph for n in figures.GRAPH_EXAMPLE):  # the small test graphs lack the example nodes
        plt.close(figures.graph_example(graph, images))
    if len(families.select_families(graph)) == 6:
        plt.close(figures.real_families(graph, images))
    else:
        print("graph too small for the six real families: hero figure skipped")
    if len(families.select_disagreements(result["distance_ignoring_alleles"])) == 12:
        plt.close(figures.name_disagreement_gallery(result["distance_ignoring_alleles"], images))
    else:
        print("graph too small for the 12 name disagreements: gallery skipped")
    if args.manifest:
        root = Path(__file__).resolve().parents[1]
        pin = json.loads((root / "chamredb" / "data" / "graph_release.json").read_text())
        if pin["sha256"] != sha:
            print("pipeline diagram skipped: the graph is not the one pinned in graph_release.json")
        else:
            plt.close(figures.pipeline_diagram(json.loads(Path(args.manifest).read_text()), stats, pin,
                                               (root / "builder" / "Snakefile").read_text(),
                                               (root / "builder" / "environment.yml").read_text(), images))
    if args.manifest and json.loads((Path(__file__).resolve().parents[1] / "chamredb" / "data" / "graph_release.json").read_text())["sha256"] == sha:
        from analysis import screenshots
        for name, database, identifier in (("figure_query_qnrB2", "ncbi", "qnrB2"), ("figure_query_aac3IIIb", "resfinder", "aac(3)-IIIb")):
            command, markups = screenshots.capture_query(identifier, database)   # the default graph, as a user's first query
            plt.close(screenshots.render_screenshot(command, markups, images, name))
    print(f"wrote {out}/statistics.json, two tables and the figures for {graph.number_of_nodes()} nodes")


if __name__ == "__main__":
    main()
