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


def test_figure_on_the_bundled_release_graph_matches_the_committed_statistics(tmp_path):
    result, stats, _ = build(BUNDLED, tmp_path)
    committed = json.loads((ROOT / "docs" / "statistics.json").read_text())
    assert stats["analysis"] == committed["analysis"]
    assert pd.read_csv(tmp_path / f"source_data_{NAME}_C_density.csv")["matches"].sum() == committed["analysis"]["matches_analysed"]
