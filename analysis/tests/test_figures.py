"""Shared figure style and the name-similarity vs identity figure (ticket #39, spec #34).

The picture is reviewed by the maintainer; these tests check what can be checked: sizes in mm, font sizes, the palette,
that the source-data CSVs equal the statistics of #38, that numbers written on the figure equal the statistics, and that
two renders give byte-identical files.
"""
import json
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from analysis import figures, matches, run, style  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SUBSET = ROOT / "chamredb" / "tests" / "test_data" / "graph_2021_subset.json"
BUNDLED = ROOT / "chamredb" / "data" / "graph.json"
NAME = "figure_name_similarity_vs_identity"


def test_palette_is_okabe_ito():
    assert style.OKABE_ITO == ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
    assert style.DATABASE_COLOURS.keys() == {"card", "ncbi", "resfinder"}
    assert len(set(style.DATABASE_COLOURS.values())) == 3


def test_style_sets_fonts_of_at_least_7_pt_and_editable_text():
    style.apply_style()
    for key in ("font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize", "ytick.labelsize", "legend.fontsize"):
        assert matplotlib.rcParams[key] >= 7, key
    assert matplotlib.rcParams["pdf.fonttype"] == 42 and matplotlib.rcParams["ps.fonttype"] == 42


def test_figure_size_is_given_in_millimetres():
    fig = style.new_figure(180, 150)
    assert tuple(round(v * 25.4, 6) for v in fig.get_size_inches()) == (180.0, 150.0)


def build(graph_path, out):
    graph = matches.read_graph(graph_path)
    result = matches.assess(graph)
    stats = run.statistics(graph, result, graph_name=Path(graph_path).name, graph_sha256="x")
    fig = figures.name_similarity_vs_identity(result, stats, out)
    return result, stats, fig


def test_source_data_equals_the_statistics(tmp_path):
    result, stats, _ = build(SUBSET, tmp_path)
    a = stats["analysis"]
    density = pd.read_csv(tmp_path / f"source_data_{NAME}_C_density.csv")
    assert density["matches"].sum() == a["matches_analysed"] == len(result["distance"])
    for column, key in (("A_name_similarity", "name_similarity"), ("B_sequence_identity", "sequence_identity")):
        hist = pd.read_csv(tmp_path / f"source_data_{NAME}_{column[0]}_histogram.csv")
        assert hist["matches"].sum() == a["matches_analysed"]
    fit = pd.read_csv(tmp_path / f"source_data_{NAME}_C_fit.csv")
    full = fit[fit["analysis"] == "full_names"]
    slope = (full["fit"].iloc[-1] - full["fit"].iloc[0]) / (full["sequence_identity"].iloc[-1] - full["sequence_identity"].iloc[0])
    assert slope == pytest.approx(a["full_names"]["slope"], abs=1e-9)
    assert (full["lower_95"] <= full["fit"]).all() and (full["fit"] <= full["upper_95"]).all()


def test_numbers_written_on_the_figure_equal_the_statistics(tmp_path):
    _, stats, fig = build(SUBSET, tmp_path)
    texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
    a = stats["analysis"]
    assert f"{a['full_names']['predicted_similarity_at_threshold']:.2f}" in texts
    assert f"{a['without_alleles']['predicted_similarity_at_threshold']:.2f}" in texts
    assert f"n = {a['matches_analysed']:,}" in texts
    density_ax = fig.axes[2]
    assert density_ax.get_xlim() == (0, 1) and density_ax.get_ylim() == (0, 1)


def test_two_renders_are_byte_identical(tmp_path):
    build(SUBSET, tmp_path / "a")
    build(SUBSET, tmp_path / "b")
    files = sorted(p.name for p in (tmp_path / "a").iterdir())
    assert f"{NAME}.png" in files and f"{NAME}.svg" in files and len(files) == 6
    for name in files:
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


def assert_same_statistics(got, expected, path="analysis"):
    """Integers and strings exactly; floats within 1e-9. Another platform's numpy gives regression coefficients that agree to
    about 1e-13 but not bit for bit, so exact equality is only expected between runs in one environment."""
    if isinstance(expected, dict):
        assert got.keys() == expected.keys(), path
        for key in expected:
            assert_same_statistics(got[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, float):
        assert got == pytest.approx(expected, abs=1e-9), path
    else:
        assert got == expected, path


def test_figure_on_the_bundled_release_graph_matches_the_committed_statistics(tmp_path):
    result, stats, _ = build(BUNDLED, tmp_path)
    committed = json.loads((ROOT / "docs" / "statistics.json").read_text())
    assert_same_statistics(stats["analysis"], committed["analysis"])
    assert pd.read_csv(tmp_path / f"source_data_{NAME}_C_density.csv")["matches"].sum() == committed["analysis"]["matches_analysed"]


# ---- overlap (UpSet) figure, ticket #49 ----

OVERLAP = "figure_overlap_upset"


def hand_group_stats():
    from analysis.tests.test_groups import HAND
    return run.group_statistics(HAND)


def test_overlap_source_data_equal_the_group_statistics_on_the_hand_graph(tmp_path):
    stats = hand_group_stats()
    figures.overlap_upset(stats, tmp_path)
    inter = pd.read_csv(tmp_path / f"source_data_{OVERLAP}_intersections.csv")
    assert inter["groups"].sum() == stats["total"] == 8 and inter["proteins"].sum() == 24 and len(inter) == 7
    assert inter.iloc[0]["groups"] == 2 and inter.iloc[0]["databases"] == "card+ncbi"   # largest first
    row = inter[inter["databases"] == "card+resfinder"].iloc[0]
    assert (row["card"], row["ncbi"], row["resfinder"], row["groups"], row["proteins"]) == (1, 0, 1, 1, 2)
    sets = pd.read_csv(tmp_path / f"source_data_{OVERLAP}_set_sizes.csv").set_index("database")["groups"].to_dict()
    # groups containing each database, by hand: card 1+2+1+1 = 5 (G1, G2, G3, G5, G6 -> 5), ncbi G1 G2 G3 G4 G8 = 5, resfinder G1 G6 G7 G8 = 4
    assert sets == {"card": 5, "ncbi": 5, "resfinder": 4}


def test_overlap_figure_numbers_size_and_determinism(tmp_path):
    stats = hand_group_stats()
    fig = figures.overlap_upset(stats, tmp_path / "a")
    figures.overlap_upset(stats, tmp_path / "b")
    assert tuple(round(v * 25.4, 6) for v in fig.get_size_inches()) == (180.0, 105.0)
    texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
    assert "8 groups" in texts and "24 proteins" in texts
    for count in stats["by_databases"].values():
        assert str(count) in texts
    for name in sorted(p.name for p in (tmp_path / "a").iterdir()):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


def test_overlap_figure_on_the_release_graph_matches_the_committed_statistics(tmp_path):
    graph = matches.read_graph(BUNDLED)
    stats = run.group_statistics(graph)
    committed = json.loads((ROOT / "docs" / "statistics.json").read_text())["groups"]
    assert stats == committed
    figures.overlap_upset(stats, tmp_path)
    inter = pd.read_csv(tmp_path / f"source_data_{OVERLAP}_intersections.csv")
    assert inter["groups"].sum() == committed["total"] and inter["proteins"].sum() == graph.number_of_nodes()


# ---- real-gene hero figure, ticket #41 ----

HERO = "figure_real_families"


def test_hero_source_data_is_read_back_from_the_graph_file(tmp_path):
    graph = matches.read_graph(BUNDLED)
    figures.real_families(graph, tmp_path)
    raw = json.loads(BUNDLED.read_text())                          # independent of networkx
    names = {n["id"]: n["name"] for n in raw["nodes"]}
    edges = {(e["source"], e["target"]): e for e in raw["links"]}
    nodes = pd.read_csv(tmp_path / f"source_data_{HERO}_nodes.csv")
    drawn = pd.read_csv(tmp_path / f"source_data_{HERO}_edges.csv")
    assert nodes["panel"].nunique() == drawn["panel"].nunique() == 6
    for row in nodes.itertuples():
        assert names[row.node_id] == row.name and row.node_id.split(":")[0] == row.database
    for row in drawn.itertuples():
        e = edges[(row.source, row.target)]
        assert (e["type"], e["identity"], e["coverage"]) == (row.type, row.identity, row.coverage)


def test_hero_figure_names_identities_size_and_determinism(tmp_path):
    graph = matches.read_graph(BUNDLED)
    fig = figures.real_families(graph, tmp_path / "a")
    figures.real_families(graph, tmp_path / "b")
    assert tuple(round(v * 25.4, 6) for v in fig.get_size_inches()) == (180.0, 100.0)
    texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
    nodes = pd.read_csv(tmp_path / "a" / f"source_data_{HERO}_nodes.csv")
    for name in nodes["name"]:
        assert name in texts, name
    assert "0.48" in texts                                         # the one-way hit identity written on its edge
    for name in sorted(p.name for p in (tmp_path / "a").iterdir()):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


# ---- distributions and the name disagreement gallery, ticket #40 ----

DIST = "figure_name_distributions"
GALLERY = "figure_name_disagreements"


def test_distribution_summary_equals_the_statistics_and_an_independent_median(tmp_path):
    import statistics as py_statistics
    graph = matches.read_graph(BUNDLED)
    result = matches.assess(graph)
    stats = run.statistics(graph, result, graph_name="graph.json", graph_sha256="x")
    figures.name_distributions(result, stats, tmp_path)
    summary = pd.read_csv(tmp_path / f"source_data_{DIST}_summary.csv").set_index("panel")
    a = stats["analysis"]
    assert summary.loc["A", "n"] == summary.loc["C", "n"] == a["full_names"]["matches_at_or_above_threshold"]
    assert summary.loc["B", "n"] == summary.loc["D", "n"] == a["without_alleles"]["matches_at_or_above_threshold"]
    high = result["high_identity"]
    high_free = result["high_identity_without_alleles"]
    for panel, values in (("A", high["name_similarity"]), ("B", high_free["name_similarity_without_alleles"]),
                          ("C", high["similarity_difference"]), ("D", high_free["similarity_difference"])):
        assert summary.loc[panel, "median"] == pytest.approx(py_statistics.median(values.tolist()), abs=1e-12)
        assert summary.loc[panel, "min"] <= summary.loc[panel, "p2_5"] <= summary.loc[panel, "median"] <= summary.loc[panel, "p97_5"] <= summary.loc[panel, "max"]
    hist = pd.read_csv(tmp_path / f"source_data_{DIST}_histogram.csv")
    for panel in "ABCD":
        assert hist[hist["panel"] == panel]["matches"].sum() == summary.loc[panel, "n"]


def test_distribution_figure_numbers_size_and_determinism(tmp_path):
    graph = matches.read_graph(BUNDLED)
    result = matches.assess(graph)
    stats = run.statistics(graph, result, graph_name="graph.json", graph_sha256="x")
    fig = figures.name_distributions(result, stats, tmp_path / "a")
    figures.name_distributions(result, stats, tmp_path / "b")
    assert tuple(round(v * 25.4, 6) for v in fig.get_size_inches()) == (180.0, 115.0)
    texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
    a = stats["analysis"]
    for block in ("full_names", "without_alleles"):
        n = a[block]["matches_at_or_above_threshold"]
        share = a[block]["low_name_similarity_rows"] / n
        assert f"n = {n:,}" in texts
        assert f"{share:.1%} below {a[block]['predicted_similarity_at_threshold']:.2f}" in texts
    for name in sorted(p.name for p in (tmp_path / "a").iterdir()):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


def test_gallery_rows_numbers_size_and_determinism(tmp_path):
    from analysis import families
    graph = matches.read_graph(BUNDLED)
    table = matches.assess(graph)["distance_ignoring_alleles"]
    chosen = families.select_disagreements(table, per_pair=4)
    fig = figures.name_disagreement_gallery(table, tmp_path / "a")
    figures.name_disagreement_gallery(table, tmp_path / "b")
    assert tuple(round(v * 25.4, 6) for v in fig.get_size_inches()) == (180.0, 95.0)
    rows = pd.read_csv(tmp_path / "a" / f"source_data_{GALLERY}_rows.csv")
    assert rows[["name_a", "name_b"]].values.tolist() == chosen[["name_a", "name_b"]].values.tolist()
    eligible = pd.read_csv(tmp_path / "a" / f"source_data_{GALLERY}_eligible.csv").iloc[0]
    everything = families.select_disagreements(table, per_pair=10 ** 9)
    assert eligible["eligible_matches"] == len(everything)
    assert eligible["with_similarity_0"] == int((everything["name_similarity_without_alleles"] == 0).sum())
    texts = " ".join(t.get_text() for t in fig.texts) + " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
    for name in list(rows["name_a"]) + list(rows["name_b"]):
        assert name in texts, name
    assert f"{int(eligible['eligible_matches']):,}" in texts
    for name in sorted(p.name for p in (tmp_path / "a").iterdir()):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


def test_distribution_quantiles_on_a_small_untied_example_worked_out_by_hand(tmp_path):
    # the real data are so tied that median, p25 and p75 coincide in several panels, so they cannot catch a wrong quantile;
    # five untied values per panel, with linear-interpolation percentiles worked out by hand:
    #   [0.1, 0.2, 0.3, 0.4, 0.9]: min 0.1, p2.5 = 0.1 + 0.1 * 0.1 = 0.11, p25 0.2, median 0.3, p75 0.4, p97.5 = 0.4 + 0.9 * 0.5 = 0.85, max 0.9
    #   [0.05, 0.1, 0.2, 0.4, 0.8]: p2.5 = 0.05 + 0.1 * 0.05 = 0.055, p25 0.1, median 0.2, p75 0.4, p97.5 = 0.4 + 0.9 * 0.4 = 0.76
    sim = [0.1, 0.2, 0.3, 0.4, 0.9]
    diff = [0.05, 0.1, 0.2, 0.4, 0.8]
    result = {"high_identity": pd.DataFrame({"name_similarity": sim, "similarity_difference": diff}),
              "high_identity_without_alleles": pd.DataFrame({"name_similarity_without_alleles": sim, "similarity_difference": diff})}
    block = {"matches_at_or_above_threshold": 5, "predicted_similarity_at_threshold": 0.8, "low_name_similarity_rows": 2}
    stats = {"analysis": {"identity_threshold": 0.95, "full_names": block, "without_alleles": block}}
    figures.name_distributions(result, stats, tmp_path)
    summary = pd.read_csv(tmp_path / f"source_data_{DIST}_summary.csv").set_index("panel")
    for panel, expected in (("A", [0.1, 0.11, 0.2, 0.3, 0.4, 0.85, 0.9]), ("C", [0.05, 0.055, 0.1, 0.2, 0.4, 0.76, 0.8])):
        got = summary.loc[panel, ["min", "p2_5", "p25", "median", "p75", "p97_5", "max"]].tolist()
        assert got == pytest.approx(expected, abs=1e-12), panel


# ---- pipeline diagram, ticket #50 ----

PIPE = "figure_pipeline"


def pipeline_inputs():
    manifest = {"files": [
        {"database": "card", "role": "sequences", "release": "9.9.9", "retrieved": "2030-01-02", "card_release": "9.9.9"},
        {"database": "card", "role": "ontology", "release": "9.9.9", "retrieved": "2030-01-02", "card_release": "9.9.9"},
        {"database": "ncbi", "role": "catalogue", "release": "AMRFinderPlus database 1.2, 2030-01-01.1", "retrieved": "2030-01-02"},
        {"database": "resfinder", "role": "phenotypes", "release": "abcdef1234567890", "retrieved": "2030-01-02"}]}
    stats = {"graph": {"nodes": 1234, "edges": 5678, "edges_by_type": {"OWH": 3000, "RBH": 2678}, "rbh_pairs": 1339,
                       "nodes_by_database": {"card": 500, "ncbi": 600, "resfinder": 134}},
             "groups": {"total": 900}}
    pin = {"data_release": "2030-01-02", "release_tag": "graph-2030.01.02", "sha256": "0123456789abcdef" * 4}
    snakefile = ("rule mmseqs_easy_rbh:\n\tparams:\n\t\tsensitivity = 6.5\nrule mmseqs_search:\n\tparams:\n\t\tsensitivity = 5.5\n\t\tmax_accept = 4\n")
    environment = "dependencies:\n    - mmseqs2=99.abcde\n    - python=3.9\n"
    return manifest, stats, pin, snakefile, environment


def test_pipeline_facts_on_hand_inputs_worked_out_by_hand():
    facts = figures.pipeline_facts(*pipeline_inputs())
    got = {(f["box"], f["item"]): f["value"] for f in facts}
    assert got[("inputs", "card")] == "CARD 9.9.9"
    assert got[("inputs", "ncbi")] == "NCBI AMRFinderPlus database 1.2, 2030-01-01.1"
    assert got[("inputs", "resfinder")] == "ResFinder db commit abcdef1"
    assert got[("match", "tool")] == "MMseqs2 99.abcde"
    assert got[("match", "rbh")] == "reciprocal best hits: easy-rbh, sensitivity 6.5"
    assert got[("match", "search")] == "one-way hits: search, sensitivity 5.5, max-accept 4"
    assert got[("graph", "nodes")] == "1,234 nodes"
    assert got[("graph", "edges")] == "5,678 edges (2,678 RBH, 3,000 one-way)"
    assert got[("graph", "groups")] == "1,339 RBH pairs, 900 equivalence groups"
    assert got[("release", "tag")] == "graph-2030.01.02"
    assert got[("release", "sha256")] == "sha256 01234567"


def test_pipeline_diagram_numbers_come_from_the_committed_files(tmp_path):
    graph_stats = json.loads((ROOT / "docs" / "statistics.json").read_text())
    manifest = json.loads((ROOT / "docs" / "data_manifest.json").read_text())
    pin = json.loads((ROOT / "chamredb" / "data" / "graph_release.json").read_text())
    snakefile = (ROOT / "builder" / "Snakefile").read_text()
    environment = (ROOT / "builder" / "environment.yml").read_text()
    fig = figures.pipeline_diagram(manifest, graph_stats, pin, snakefile, environment, tmp_path)
    assert tuple(round(v * 25.4, 6) for v in fig.get_size_inches()) == (180.0, 50.0)
    # an independent read of the same files
    g = graph_stats["graph"]
    texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
    for expected in (f"{g['nodes']:,} nodes", f"{g['rbh_pairs']:,} RBH pairs", f"{graph_stats['groups']['total']:,} equivalence groups",
                     f"{g['edges_by_type']['RBH']:,} RBH", pin["release_tag"], "sensitivity 7.5", "sensitivity 7.0", "max-accept 3", "MMseqs2 18.8cc5c",
                     "CARD 4.0.2", "AMRFinderPlus database 4.2, 2026-10-29.1", "commit eecf0aa"):
        assert expected in texts, expected
    rows = pd.read_csv(tmp_path / f"source_data_{PIPE}.csv")
    assert set(rows.columns) == {"box", "item", "value", "source"} and rows["value"].str.len().min() > 0
    assert pin["nodes"] == g["nodes"] and pin["edges"] == g["edges"]    # the pin and the statistics describe the same graph


def test_pipeline_diagram_is_byte_identical_across_renders(tmp_path):
    args = pipeline_inputs()
    figures.pipeline_diagram(*args, tmp_path / "a")
    figures.pipeline_diagram(*args, tmp_path / "b")
    for name in sorted(p.name for p in (tmp_path / "a").iterdir()):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


# ---- the graph example (the original README hero, redrawn from the graph) ----

EXAMPLE = "figure_graph_example"


def raw_node_and_edges():
    raw = json.loads(BUNDLED.read_text())
    nodes = {n["id"]: n for n in raw["nodes"]}
    edges = {(e["source"], e["target"]): e for e in raw["links"]}
    return nodes, edges


def test_graph_example_objects_equal_the_graph_file(tmp_path):
    graph = matches.read_graph(BUNDLED)
    figures.graph_example(graph, tmp_path)
    shown = json.loads((tmp_path / f"source_data_{EXAMPLE}.json").read_text())
    nodes, edges = raw_node_and_edges()
    assert [n["id"] for n in shown["nodes"]] == list(figures.GRAPH_EXAMPLE)
    for n in shown["nodes"]:
        original = {k: v for k, v in nodes[n["id"]].items() if k != "id"}
        for key, value in n.items():
            if key != "id":
                assert original[key] == value, (n["id"], key)       # every shown field is the graph's value
        assert n["database"] == n["id"].split(":")[0] and n["name"] == original["name"]
    types = {(e["source"], e["target"]): e for e in shown["edges"]}
    assert {k: (v["type"], v["identity"], v["coverage"]) for k, v in types.items()} == {
        k: (edges[k]["type"], edges[k]["identity"], edges[k]["coverage"]) for k in types}
    assert sorted(v["type"] for v in types.values()) == ["OWH", "RBH", "RBH"]


def test_graph_example_figure_names_size_and_determinism(tmp_path):
    graph = matches.read_graph(BUNDLED)
    fig = figures.graph_example(graph, tmp_path / "a")
    figures.graph_example(graph, tmp_path / "b")
    assert tuple(round(v * 25.4, 6) for v in fig.get_size_inches()) == (180.0, 100.0)
    texts = " ".join(" ".join(t.get_text().split()) for ax in fig.axes for t in ax.texts)
    for needle in ("blaB-21", "Chryseobacterium meningosepticum BlaB", "blaB-11", '"RBH"', '"OWH"', "0.976", "WP_029728367.1"):
        assert needle in texts, needle
    for name in sorted(p.name for p in (tmp_path / "a").iterdir()):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


def test_graph_example_refuses_a_graph_without_the_example(tmp_path):
    from networkx import DiGraph
    with pytest.raises(ValueError, match="example"):
        figures.graph_example(DiGraph(), tmp_path)
