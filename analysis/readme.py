"""Fill the numbers and tables of README.md and docs/appendices.md from docs/statistics.json (spec #34, ticket #43).

Numbers sit between `<!-- stat:KEY -->` and `<!-- /stat -->`; tables between `<!-- table:NAME -->` and `<!-- /table -->`.
    python -m analysis.readme [--stats docs/statistics.json]
"""
import argparse
import json
import re
from pathlib import Path

LABELS = {"card": "CARD", "ncbi": "NCBI", "resfinder": "ResFinder"}
ORDER = [("card", "ncbi"), ("card", "resfinder"), ("ncbi", "card"), ("ncbi", "resfinder"), ("resfinder", "card"), ("resfinder", "ncbi")]


def format_stat(stats, key):
    a, g, groups = stats["analysis"], stats["graph"], stats["groups"]
    high = a["full_names"]["matches_at_or_above_threshold"]
    three = groups["by_databases"]["card+ncbi+resfinder"]
    values = {
        "graph.nodes": f"{g['nodes']:,}", "graph.edges": f"{g['edges']:,}", "graph.rbh_pairs": f"{g['rbh_pairs']:,}",
        "graph.edges_rbh": f"{g['edges_by_type']['RBH']:,}", "graph.edges_owh": f"{g['edges_by_type']['OWH']:,}",
        "groups.total": f"{groups['total']:,}", "groups.all_three": f"{three:,}", "groups.all_three_share": f"{three / groups['total']:.0%}",
        "analysis.matches": f"{a['matches_analysed']:,}", "analysis.skipped": f"{a['matches_skipped_blaPAO']:,}",
        "analysis.high_identity": f"{high:,}",
        "analysis.full.predicted": f"{a['full_names']['predicted_similarity_at_threshold']:.2f}",
        "analysis.free.predicted": f"{a['without_alleles']['predicted_similarity_at_threshold']:.2f}",
        "analysis.full.low_share": f"{a['full_names']['low_name_similarity_rows'] / high:.1%}",
        "analysis.free.low_share": f"{a['without_alleles']['low_name_similarity_rows'] / high:.1%}",
        "analysis.full.low_rows": f"{a['full_names']['low_name_similarity_rows']:,}",
        "analysis.free.low_rows": f"{a['without_alleles']['low_name_similarity_rows']:,}",
    }
    return values[key]


def table(stats, name):
    rows = []
    if name == "rbh":
        rows += ["| Comparison | Proteins of the first database in the graph | With a reciprocal best hit |", "| --- | --- | --- |"]
        for a, b in ORDER:
            c = stats["pairs"][f"{a}>{b}"]
            rows.append(f"| {LABELS[a]} vs {LABELS[b]} | {c['source_nodes']:,} | {c['with_rbh']:,} ({c['with_rbh'] / c['source_nodes']:.1%}) |")
    elif name == "search":
        rows += ["| Comparison | With a one-way hit only | Covered by a reciprocal best hit or a one-way hit |", "| --- | --- | --- |"]
        for a, b in ORDER:
            c = stats["pairs"][f"{a}>{b}"]
            rows.append(f"| {LABELS[a]} vs {LABELS[b]} | {c['with_search_only']:,} | {c['covered']:,} ({c['covered'] / c['source_nodes']:.1%}) |")
    else:
        raise KeyError(name)
    return "\n".join(rows)


def render(text, stats):
    text = re.sub(r"<!-- stat:([\w.]+) -->.*?<!-- /stat -->", lambda m: f"<!-- stat:{m.group(1)} -->{format_stat(stats, m.group(1))}<!-- /stat -->", text)
    return re.sub(r"<!-- table:(\w+) -->\n.*?\n<!-- /table -->", lambda m: f"<!-- table:{m.group(1)} -->\n{table(stats, m.group(1))}\n<!-- /table -->", text, flags=re.S)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--stats", default="docs/statistics.json")
    args = parser.parse_args(argv)
    stats = json.loads(Path(args.stats).read_text())
    for name in ("README.md", "docs/appendices.md"):
        path = Path(name)
        path.write_text(render(path.read_text(), stats))


if __name__ == "__main__":
    main()
