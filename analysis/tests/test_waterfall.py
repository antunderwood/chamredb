"""Release waterfall: the four accepted data steps from the 2021 graph to graph-2026.11.09 (ticket #51, spec #34).

Level 1: a hand-built pair of graphs whose differences are worked out by hand. Level 5: the committed step summary
(docs/waterfall_steps.json) adds up exactly: each step starts where the previous one ended, start + gained - lost = end, and the end
is the bundled graph. The counts also equal the numbers of the reviewed reports (pull requests 28 to 31), typed in here from those reports.
"""
import json
from pathlib import Path

import pytest

from analysis import figures, waterfall

ROOT = Path(__file__).resolve().parents[2]
STEPS = ROOT / "docs" / "waterfall_steps.json"


def node(i, name):
    return {"id": i, "name": name}


def edge(s, t, ty, ident, cov):
    return {"source": s, "target": t, "type": ty, "identity": ident, "coverage": cov}


OLD = {"nodes": [node("card:ARO:1", "a"), node("card:ARO:2", "b"), node("card:ARO:3", "c"), node("ncbi:N1", "x"), node("ncbi:N2", "y"), node("resfinder:r1", "r")],
       "links": [edge("card:ARO:1", "ncbi:N1", "RBH", 0.95, 1.0), edge("ncbi:N1", "card:ARO:1", "RBH", 0.95, 1.0),
                 edge("card:ARO:2", "resfinder:r1", "OWH", 0.80, 0.9), edge("ncbi:N2", "resfinder:r1", "OWH", 0.70, 0.9),
                 edge("card:ARO:3", "ncbi:N1", "OWH", 0.85, 0.9)]}
NEW = {"nodes": [node("card:ARO:1", "a"), node("card:ARO:2", "B2"), node("card:ARO:30", "c"), node("ncbi:N1", "x"), node("ncbi:N9", "z"), node("resfinder:r1", "r")],
       "links": [edge("card:ARO:1", "ncbi:N1", "RBH", 0.95, 1.0), edge("ncbi:N1", "card:ARO:1", "RBH", 0.97, 1.0),
                 edge("card:ARO:2", "resfinder:r1", "RBH", 0.80, 0.9), edge("resfinder:r1", "card:ARO:2", "RBH", 0.80, 0.9),
                 edge("card:ARO:30", "ncbi:N1", "OWH", 0.85, 0.9), edge("ncbi:N9", "resfinder:r1", "OWH", 0.60, 0.9)]}


def test_transition_counts_on_the_hand_graphs():
    t = waterfall.transition(OLD, NEW)
    # worked out by hand: nodes 6 -> 6 (gained card:ARO:30 and ncbi:N9; lost card:ARO:3 and ncbi:N2); edges 5 -> 6 (gained 3, lost 2);
    # one type change (card:ARO:2 -> resfinder:r1 OWH to RBH); ids absent: card:ARO:3 and ncbi:N2; renamed: card:ARO:2 (b to B2)
    assert (t["nodes_before"], t["nodes_after"], t["nodes_gained"], t["nodes_lost"]) == (6, 6, 2, 2)
    assert (t["edges_before"], t["edges_after"], t["edges_gained"], t["edges_lost"]) == (5, 6, 3, 2)
    assert t["edges_type_changed"] == 1 and t["ids_absent"] == 2 and t["ids_renamed"] == 1
    assert t["ids_absent_by_database"] == {"card": 1, "ncbi": 1}


def test_a_transition_adds_up():
    t = waterfall.transition(OLD, NEW)
    assert t["nodes_before"] + t["nodes_gained"] - t["nodes_lost"] == t["nodes_after"]
    assert t["edges_before"] + t["edges_gained"] - t["edges_lost"] == t["edges_after"]


@pytest.fixture(scope="module")
def steps():
    return json.loads(STEPS.read_text())


def test_the_committed_steps_chain_and_add_up(steps):
    transitions = steps["transitions"]
    assert [t["to"] for t in transitions] == ["MMseqs2 18", "CARD 4.0.2", "NCBI 2026-10-29.1", "ResFinder 2.6.0"]
    assert (transitions[0]["nodes_before"], transitions[0]["edges_before"]) == (11718, 33632)         # the shipped 2021 graph
    for previous, current in zip(transitions, transitions[1:]):
        assert current["nodes_before"] == previous["nodes_after"] and current["edges_before"] == previous["edges_after"]
    for t in transitions:
        assert t["nodes_before"] + t["nodes_gained"] - t["nodes_lost"] == t["nodes_after"], t["to"]
        assert t["edges_before"] + t["edges_gained"] - t["edges_lost"] == t["edges_after"], t["to"]
    pin = json.loads((ROOT / "chamredb" / "data" / "graph_release.json").read_text())
    assert (transitions[-1]["nodes_after"], transitions[-1]["edges_after"]) == (pin["nodes"], pin["edges"])
    assert steps["steps"][-1]["graph_sha256"] == pin["sha256"]


def test_the_counts_equal_the_reviewed_reports(steps):
    # numbers typed from the reports reviewed in pull requests 28 (MMseqs2), 29 (CARD), 30 (NCBI) and 31 (ResFinder)
    reviewed = {
        "MMseqs2 18":        dict(nodes_gained=0, nodes_lost=0, edges_gained=802, edges_lost=802, edges_type_changed=0, ids_absent=0, ids_renamed=0),
        "CARD 4.0.2":        dict(nodes_gained=1492, nodes_lost=17, edges_gained=8762, edges_lost=1112, edges_type_changed=29, ids_absent=17, ids_renamed=108),
        "NCBI 2026-10-29.1": dict(nodes_gained=2639, nodes_lost=31, edges_gained=14738, edges_lost=4069, edges_type_changed=15, ids_absent=31, ids_renamed=23),
        "ResFinder 2.6.0":   dict(nodes_gained=100, nodes_lost=42, edges_gained=2876, edges_lost=2856, edges_type_changed=51, ids_absent=42, ids_renamed=0),
    }
    for t in steps["transitions"]:
        for key, value in reviewed[t["to"]].items():
            assert t[key] == value, (t["to"], key)


def test_the_waterfall_figure_source_data_numbers_size_and_determinism(steps, tmp_path):
    import pandas as pd
    fig = figures.release_waterfall(steps, tmp_path / "a")
    figures.release_waterfall(steps, tmp_path / "b")
    assert tuple(round(v * 25.4, 6) for v in fig.get_size_inches()) == (180.0, 110.0)
    rows = pd.read_csv(tmp_path / "a" / "source_data_figure_release_waterfall.csv")
    assert rows["quantity"].nunique() == 2 and len(rows) == 2 * 6          # start, four steps, end for nodes and for edges
    texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
    for t in steps["transitions"]:
        assert f"+{t['nodes_gained']:,}" in texts and f"+{t['edges_gained']:,}" in texts
        assert f"{t['ids_absent']:,}" in texts
    assert f"{steps['transitions'][-1]['nodes_after']:,}" in texts
    for name in sorted(p.name for p in (tmp_path / "a").iterdir()):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name
