"""The README, the appendices and the example outputs are generated from the release graph (ticket #43, spec #34 invariant 5).

Numbers in the README sit between markers `<!-- stat:KEY -->...<!-- /stat -->`; the test formats each KEY from docs/statistics.json with its
own table (not the generator's) and compares. The appendix tables sit between `<!-- table:NAME -->` markers and are compared with an independent
count from the graph file. The example TSVs in docs/data are compared with fresh output of the package on the bundled graph.
"""
import json
import re
from pathlib import Path

import pytest

from analysis import matches, readme, run

ROOT = Path(__file__).resolve().parents[2]
README = (ROOT / "README.md").read_text()
APPENDICES = (ROOT / "docs" / "appendices.md").read_text()
STATS = json.loads((ROOT / "docs" / "statistics.json").read_text())
BUNDLED = ROOT / "chamredb" / "data" / "graph.json"


def stat(key):
    a, g, grp = STATS["analysis"], STATS["graph"], STATS["groups"]
    n = a["full_names"]["matches_at_or_above_threshold"]
    table = {
        "graph.nodes": f"{g['nodes']:,}", "graph.edges": f"{g['edges']:,}", "graph.rbh_pairs": f"{g['rbh_pairs']:,}",
        "graph.edges_rbh": f"{g['edges_by_type']['RBH']:,}", "graph.edges_owh": f"{g['edges_by_type']['OWH']:,}",
        "groups.total": f"{grp['total']:,}", "groups.all_three": f"{grp['by_databases']['card+ncbi+resfinder']:,}",
        "groups.all_three_share": f"{grp['by_databases']['card+ncbi+resfinder'] / grp['total']:.0%}",
        "analysis.matches": f"{a['matches_analysed']:,}", "analysis.skipped": f"{a['matches_skipped_blaPAO']:,}",
        "analysis.high_identity": f"{n:,}",
        "analysis.full.predicted": f"{a['full_names']['predicted_similarity_at_threshold']:.2f}",
        "analysis.free.predicted": f"{a['without_alleles']['predicted_similarity_at_threshold']:.2f}",
        "analysis.full.low_share": f"{a['full_names']['low_name_similarity_rows'] / n:.1%}",
        "analysis.free.low_share": f"{a['without_alleles']['low_name_similarity_rows'] / n:.1%}",
        "analysis.full.low_rows": f"{a['full_names']['low_name_similarity_rows']:,}",
        "analysis.free.low_rows": f"{a['without_alleles']['low_name_similarity_rows']:,}",
    }
    return table[key]


MARKER = re.compile(r"<!-- stat:([\w.]+) -->(.*?)<!-- /stat -->")


def test_every_number_in_the_readme_equals_the_statistics_file():
    found = MARKER.findall(README)
    assert len(found) >= 12
    for key, value in found:
        assert value == stat(key), key


def test_the_readme_renders_to_itself():
    assert readme.render(README, STATS) == README       # nothing stale: regenerating changes nothing


@pytest.mark.parametrize("stale", ["19132", "19,132", "0.69", "0.86", "scripts/", "chamredb_network", "qnrB2.png)"])
def test_the_readme_has_no_stale_numbers_or_old_files(stale):
    assert stale not in README.replace("figure_query_qnrB2.png)", "")


def test_every_local_image_and_link_in_the_readme_exists():
    for target in re.findall(r"\]\((?!http|#|mailto)([^)\s]+)\)", README):
        assert (ROOT / target.split("#")[0].lstrip("/")).exists(), target   # an anchor (file.md#section) is checked as its file


def independent_pair_counts():
    raw = json.loads(BUNDLED.read_text())
    database = lambda i: i.split(":")[0]
    nodes = {}
    for n in raw["nodes"]:
        nodes.setdefault(database(n["id"]), set()).add(n["id"])
    rbh, owh = {}, {}
    for e in raw["links"]:
        key = (database(e["source"]), database(e["target"]))
        (rbh if e["type"] == "RBH" else owh).setdefault(key, set()).add(e["source"])
    out = {}
    for a in nodes:
        for b in nodes:
            if a != b:
                r, o = rbh.get((a, b), set()), owh.get((a, b), set())
                out[(a, b)] = {"nodes": len(nodes[a]), "rbh": len(r), "owh_only": len(o - r), "covered": len(r | o)}
    return out


def table_block(name):
    return re.search(rf"<!-- table:{name} -->\n(.*?)\n<!-- /table -->", APPENDICES, re.S).group(1)


def test_appendix_tables_equal_an_independent_count_on_the_graph():
    counts = independent_pair_counts()
    labels = {"card": "CARD", "ncbi": "NCBI", "resfinder": "ResFinder"}
    for name, cell in (("rbh", lambda c: f"{c['rbh']:,} ({c['rbh'] / c['nodes']:.1%})"),
                       ("search", lambda c: f"{c['owh_only']:,} | {c['covered']:,} ({c['covered'] / c['nodes']:.1%})")):
        block = table_block(name)
        for (a, b), c in counts.items():
            row = next(line for line in block.splitlines() if line.startswith(f"| {labels[a]} vs {labels[b]} |"))
            assert cell(c) in row, (name, a, b)
        assert len(block.splitlines()) == 2 + 6


def test_pair_statistics_equal_the_independent_count():
    stats = run.statistics(matches.read_graph(BUNDLED), matches.assess(matches.read_graph(BUNDLED)), "graph.json", "x")
    for (a, b), c in independent_pair_counts().items():
        assert stats["pairs"][f"{a}>{b}"] == {"source_nodes": c["nodes"], "with_rbh": c["rbh"], "with_search_only": c["owh_only"], "covered": c["covered"]}


def test_docs_examples_equal_fresh_output_on_the_bundled_graph(tmp_path):
    import subprocess
    import sys
    for name, args in (("card_ids.tsv", ["-d", "card", "-f", "docs/data/card_ids.txt"]),
                       ("hamronize_summary.tsv", ["-j", "docs/data/hamronize_summary.json"])):
        out = tmp_path / name
        subprocess.run([sys.executable, "-m", "chamredb.run_chamredb", "query", *args, "-o", str(out)], check=True, cwd=ROOT, capture_output=True)
        assert out.read_bytes() == (ROOT / "docs" / "data" / name).read_bytes(), name
