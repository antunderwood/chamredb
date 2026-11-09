"""Characterisation harness for the analysis library (ticket #35, spec #34).

Two oracles:
1. Level 1: a hand-built 12-node graph whose best matches, name similarities and regression are worked out by hand
   in this file (exact fractions, regression by closed-form least squares, not numpy).
2. Level 6: the unmodified scripts/assess_matches.py run on the committed 2021 subset graph
   (chamredb/tests/test_data/graph_2021_subset.json); its two CSVs and statistics are frozen in analysis/tests/expected.
   Characterisation proves "unchanged", not "correct"; the defects it preserves are listed in README.md next to it.
"""
from fractions import Fraction
from pathlib import Path

import pandas as pd
import pytest
from networkx.readwrite import json_graph

from analysis import matches

HERE = Path(__file__).resolve().parent
SUBSET = HERE.parents[1] / "chamredb" / "tests" / "test_data" / "graph_2021_subset.json"


def node(i, name):
    return {"id": i, "name": name, "database": i.split(":")[0]}


def edge(s, t, ty, identity, coverage):
    return {"source": s, "target": t, "type": ty, "identity": identity, "coverage": coverage}


def pair(a, b, identity, coverage=1.0):
    return [edge(a, b, "RBH", identity, coverage), edge(b, a, "RBH", identity, coverage)]


# 12 nodes: card c1..c4, ncbi n1..n4, resfinder r1..r4
NODES = [
    node("card:c1", "qnrB2"), node("card:c2", "Staphylococcus aureus mupB conferring resistance to mupirocin"),
    node("card:c3", "blaTEM-1"), node("card:c4", "vanA"),
    node("ncbi:n1", "qnrB20"), node("ncbi:n2", "mupB"), node("ncbi:n3", "TEM-1"), node("ncbi:n4", "vanA"),
    node("resfinder:r1", "qnrB2"), node("resfinder:r2", "TEM-1"), node("resfinder:r3", "vanHAX"), node("resfinder:r4", "blaPAO"),
]
LINKS = (
    pair("card:c1", "ncbi:n1", 0.95) + pair("card:c1", "resfinder:r1", 1.0) + pair("card:c2", "ncbi:n2", 0.99)
    + pair("card:c3", "ncbi:n3", 1.0) + pair("card:c3", "resfinder:r2", 1.0) + pair("card:c4", "ncbi:n4", 1.0)
    # one-way hits from c4 to resfinder, in this order: r3 first (identity 0.9, coverage 0.95), then r2 (identity 0.95,
    # coverage 0.9). The existing rule keeps r3: a later target must beat BOTH the best identity and the best coverage.
    + [edge("card:c4", "resfinder:r3", "OWH", 0.9, 0.95), edge("card:c4", "resfinder:r2", "OWH", 0.95, 0.9)]
    # n3 reaches the node named blaPAO, which the analysis skips and counts
    + [edge("ncbi:n3", "resfinder:r4", "OWH", 0.7, 0.9)]
)
GRAPH = json_graph.node_link_graph(
    {"directed": True, "multigraph": False, "graph": {}, "nodes": NODES, "links": LINKS}, edges="links")

# (source id, target database, identity, name similarity as an exact fraction), worked out by hand:
# cleaned names: qnrB2 -> qnrb2; qnrB20 -> qnrb20; the CARD species name -> mupb (first species pattern); blaTEM-1 and TEM-1 -> tem1
# (so the bla prefix removal matters: without it these pairs would score 4/7);
# vanA -> vana; vanHAX -> vanhax. Similarity = 1 - Levenshtein distance / longer length:
# qnrb2 vs qnrb20: distance 1, length 6 -> 5/6; vana vs vanhax: distance 2 (insert h, insert x), length 6 -> 2/3; the rest identical -> 1.
HAND = [
    ("card:c1", "ncbi", "0.95", Fraction(5, 6)), ("card:c1", "resfinder", "1.0", Fraction(1)),
    ("card:c2", "ncbi", "0.99", Fraction(1)),
    ("card:c3", "ncbi", "1.0", Fraction(1)), ("card:c3", "resfinder", "1.0", Fraction(1)),
    ("card:c4", "ncbi", "1.0", Fraction(1)), ("card:c4", "resfinder", "0.9", Fraction(2, 3)),
    ("ncbi:n1", "card", "0.95", Fraction(5, 6)), ("ncbi:n2", "card", "0.99", Fraction(1)),
    ("ncbi:n3", "card", "1.0", Fraction(1)), ("ncbi:n4", "card", "1.0", Fraction(1)),
    ("resfinder:r1", "card", "1.0", Fraction(1)), ("resfinder:r2", "card", "1.0", Fraction(1)),
]


def test_best_matches_and_name_similarities_on_the_hand_graph():
    table, skipped = matches.best_matches(GRAPH)
    assert skipped == 1  # n3 -> blaPAO
    assert len(table) == len(HAND) == 13
    got = {(r.source_id, r.target_database): (r.sequence_identity, r.name_similarity) for r in table.itertuples()}
    assert set(got) == {(s, d) for s, d, _, _ in HAND}
    for source, database, identity, similarity in HAND:
        assert got[(source, database)][0] == pytest.approx(float(identity), abs=1e-12)
        assert got[(source, database)][1] == pytest.approx(float(similarity), abs=1e-12), (source, database)


def test_best_match_must_beat_both_best_identity_and_best_coverage_in_edge_order():
    table, _ = matches.best_matches(GRAPH)
    row = table[(table.source_id == "card:c4") & (table.target_database == "resfinder")].iloc[0]
    assert row.target_name == "vanHAX" and row.sequence_identity == pytest.approx(0.9)  # not blaTEM-1 at 0.95


def test_regression_equals_closed_form_least_squares_on_the_hand_table():
    table, _ = matches.best_matches(GRAPH)
    fit = matches.assess(GRAPH)["linear_fit"]
    xs = [Fraction(i) for _, _, i, _ in HAND]
    ys = [s for _, _, _, s in HAND]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    intercept = my - slope * mx
    assert fit[0] == pytest.approx(float(slope), abs=1e-12) and fit[1] == pytest.approx(float(intercept), abs=1e-12)
    assert matches.assess(GRAPH)["predicted_similarity"] == pytest.approx(float(slope * Fraction(95, 100) + intercept), abs=1e-12)


# ---- characterisation of the unmodified script on the 2021 subset graph (level 6) ----

@pytest.fixture(scope="module")
def result():
    # legacy=True: the behaviour of the unmodified script; ADR 0006 changed two rules (bla prefix, allele-free regression)
    return matches.assess(matches.read_graph(SUBSET), legacy=True)


def test_statistics_equal_the_unmodified_script(result):
    # frozen 2026-11-09 from scripts/assess_matches.py run unmodified (networkx 3.4.2, pandas 3.1.1, numpy 2.4.6, strsimpy 0.4.1)
    assert len(result["distance"]) == 264 and result["skipped_blaPAO"] == 0
    assert result["linear_fit"] == pytest.approx([1.255012027241102, -0.3092489903689918], abs=1e-12)
    assert result["linear_fit_without_alleles"] == pytest.approx([1.2922737244509965, -0.3152604168172698], abs=1e-12)
    assert result["predicted_similarity"] == pytest.approx(0.8830124355100551, abs=1e-12)
    # preserved defect: the allele-free threshold is predicted from the full-name regression, so it is the same number
    assert result["predicted_similarity_without_alleles"] == result["predicted_similarity"]
    assert (len(result["high_identity"]), len(result["low_name_similarity"]), len(result["low_name_similarity_without_alleles"])) == (205, 26, 15)


@pytest.mark.parametrize("key,name", [("low_name_similarity", "low_name_similarity.csv"),
                                       ("low_name_similarity_without_alleles", "low_name_similarity_without_alleles.csv")])
def test_csv_tables_equal_the_unmodified_script(result, key, name):
    expected = pd.read_csv(HERE / "expected" / name)
    got = result[key].reset_index(drop=True)
    pd.testing.assert_frame_equal(got[expected.columns.tolist()], expected, check_exact=False, rtol=1e-12, atol=1e-12)
