"""Terminal screenshots rendered from the real `chamredb query` output (ticket #42, spec #34 invariant 7).

The text is captured from the package's own query function (nothing is typed in), so the picture and the tool cannot
disagree. The PNG is drawn by matplotlib in DejaVu Sans Mono with the CLI's colours; the emoji icons, which no bundled font
has, are replaced by plain symbols of the same cell width (the SVG beside it keeps the real emoji).
"""
import json
from pathlib import Path

import pytest
from matplotlib import font_manager
from matplotlib.ft2font import FT2Font
from networkx.readwrite import json_graph

from analysis import screenshots

ROOT = Path(__file__).resolve().parents[2]
BUNDLED = ROOT / "chamredb" / "data" / "graph.json"
PIN = json.loads((ROOT / "chamredb" / "data" / "graph_release.json").read_text())


def raw_matches(database, name, threshold=0.9):
    """Independent of networkx and of the package: the targets of the node with this name that pass both thresholds."""
    raw = json.loads(BUNDLED.read_text())
    names = {n["id"]: n["name"] for n in raw["nodes"]}
    source = next(n["id"] for n in raw["nodes"] if n["database"] == database and n["name"] == name)
    return [(names[e["target"]], e["type"]) for e in raw["links"]
            if e["source"] == source and e["identity"] >= threshold and e["coverage"] >= threshold]


@pytest.mark.parametrize("database,name", [("ncbi", "qnrB2"), ("resfinder", "aac(3)-IIIb")])
def test_text_has_the_data_release_the_node_and_every_match_of_the_graph(database, name):
    command, markups = screenshots.capture_query(name, database)
    text = screenshots.plain_text(markups)
    assert command == f"chamredb query -d {database} -i {name}" or command == f'chamredb query -d {database} -i "{name}"'
    assert f"data release: {PIN['data_release']}" in text                       # invariant 7
    assert name in text
    expected = raw_matches(database, name)
    assert expected
    for target_name, edge_type in expected:
        assert target_name in text and edge_type in text


def test_screenshot_files_text_and_determinism(tmp_path):
    command, markups = screenshots.capture_query("qnrB2", "ncbi")
    fig = screenshots.render_screenshot(command, markups, tmp_path / "a", "figure_query_qnrB2")
    screenshots.render_screenshot(command, markups, tmp_path / "b", "figure_query_qnrB2")
    names = sorted(p.name for p in (tmp_path / "a").iterdir())
    assert names == ["figure_query_qnrB2.png", "figure_query_qnrB2.svg", "source_data_figure_query_qnrB2.txt"]
    for name in names:
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name
    saved = (tmp_path / "a" / "source_data_figure_query_qnrB2.txt").read_text()
    assert saved.startswith("$ chamredb query") and f"data release: {PIN['data_release']}" in saved
    import html
    import re
    svg = (tmp_path / "a" / "figure_query_qnrB2.svg").read_text()
    svg_text = "".join(html.unescape(s) for s in re.findall(r">([^<]+)</text>", svg)).replace("\xa0", " ")   # rich splits words and uses non-breaking spaces
    assert f"data release: {PIN['data_release']}" in svg_text and "🧬" in svg_text   # the SVG keeps the real emoji
    width_mm = fig.get_size_inches()[0] * 25.4
    assert round(width_mm, 6) == 165.0


def test_every_character_drawn_in_the_png_has_a_glyph(tmp_path):
    font = FT2Font(font_manager.findfont("DejaVu Sans Mono", fallback_to_default=False))
    command, markups = screenshots.capture_query("aac(3)-IIIb", "resfinder")
    drawn = set(screenshots.drawn_text(command, markups))
    missing = [c for c in drawn if c not in " \n" and font.get_char_index(ord(c)) == 0]
    assert not missing, missing


def test_a_custom_graph_is_captured_with_its_own_release_label_and_the_cli_thresholds(tmp_path):
    def node(i, name, **extra):
        return {"id": i, "database": i.split(":")[0], "name": name, **extra}
    # the query shows the best match per database, so the weak match must be the only candidate in its database
    graph = {"directed": True, "multigraph": False, "graph": {},
             "nodes": [node("card:ARO:1", "a", phenotype="p"), node("resfinder:r1", "strongmatch", phenotype="q"),
                       node("ncbi:WP_2", "weakmatch", product="w")],
             "links": [{"source": "card:ARO:1", "target": "resfinder:r1", "type": "RBH", "identity": 1.0, "coverage": 1.0},
                       {"source": "resfinder:r1", "target": "card:ARO:1", "type": "RBH", "identity": 1.0, "coverage": 1.0},
                       # identity 0.8 is below the CLI default of 0.9: it must not be shown (a lower threshold would show it)
                       {"source": "card:ARO:1", "target": "ncbi:WP_2", "type": "OWH", "identity": 0.8, "coverage": 1.0}]}
    path = tmp_path / "g.json"
    path.write_text(json.dumps(graph))
    _, markups = screenshots.capture_query("a", "card", graph_path=path)
    text = screenshots.plain_text(markups)
    assert "data release: unknown (custom graph)" in text and "strongmatch" in text and "RBH" in text
    assert "weakmatch" not in text
