"""Statistics, decisions of ADR 0006 and the run command (ticket #38, spec #34 invariants 1, 2 and 6).

Hand graph from test_matches (12 nodes): edges are 6 RBH pairs (12 directed RBH edges) plus 3 one-way hits = 15 edges;
nodes 4 per database; RBH pairs card-ncbi 4 and card-resfinder 2; 13 matches analysed, 1 skipped (blaPAO).
"""
import hashlib
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

from analysis import matches, run
from analysis.tests.test_matches import GRAPH as HAND, HAND as HAND_TABLE

ROOT = Path(__file__).resolve().parents[2]
SUBSET = ROOT / "chamredb" / "tests" / "test_data" / "graph_2021_subset.json"
BUNDLED = ROOT / "chamredb" / "data" / "graph.json"


def test_statistics_on_the_hand_graph():
    stats = run.statistics(HAND, matches.assess(HAND), graph_name="hand.json", graph_sha256="x")
    g = stats["graph"]
    assert g["nodes"] == 12 and g["nodes_by_database"] == {"card": 4, "ncbi": 4, "resfinder": 4}
    assert g["edges"] == 15 and g["edges_by_type"] == {"OWH": 3, "RBH": 12}
    assert g["rbh_pairs"] == 6 and g["rbh_pairs_by_databases"] == {"card+ncbi": 4, "card+resfinder": 2}
    a = stats["analysis"]
    assert a["matches_analysed"] == 13 and a["matches_skipped_blaPAO"] == 1 and a["identity_threshold"] == 0.95


# ---- decisions of ADR 0006 ----

def test_bla_is_removed_only_as_a_leading_prefix():
    assert matches.clean_name("cblA", "ncbi") == "cbla"          # was "c"
    assert matches.clean_name("elBla2-1", "ncbi") == "elbla21"   # was "el21"
    assert matches.clean_name("blaTEM-1", "ncbi") == "tem1"
    assert matches.clean_name("blaOXA/blaLRA13", "ncbi") == "oxa/blalra13"  # only the leading one
    assert matches.clean_name("cblA", "ncbi", legacy=True) == "c"           # the old behaviour stays reachable


def test_allele_free_threshold_comes_from_its_own_regression():
    # allele-free names on the hand graph (worked out by hand): qnrb2 and qnrb20 both become qnrb, blaTEM-1 and TEM-1 both
    # tem, mupb and vana unchanged, vanhax unchanged; so every similarity is 1 except vana vs vanhax = 2/3
    result = matches.assess(HAND)
    xs = [Fraction(i) for _, _, i, _ in HAND_TABLE]
    ys = [Fraction(2, 3) if (s, d) == ("card:c4", "resfinder") else Fraction(1) for s, d, _, _ in HAND_TABLE]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    intercept = my - slope * mx
    assert result["linear_fit_without_alleles"] == pytest.approx([float(slope), float(intercept)], abs=1e-12)
    assert result["predicted_similarity_without_alleles"] == pytest.approx(float(slope * Fraction(95, 100) + intercept), abs=1e-12)
    assert result["predicted_similarity_without_alleles"] != pytest.approx(result["predicted_similarity"], abs=1e-6)
    legacy = matches.assess(HAND, legacy=True)
    assert legacy["predicted_similarity_without_alleles"] == legacy["predicted_similarity"]


# ---- the run command ----

def run_cli(graph, out):
    subprocess.run([sys.executable, "-m", "analysis.run", "--graph", str(graph), "--out", str(out)],
                   check=True, cwd=ROOT, capture_output=True)


NAME = "figure_name_similarity_vs_identity"
FILES = ["statistics.json", "analysis/low_name_similarity.csv", "analysis/low_name_similarity_without_alleles.csv",
         f"images/{NAME}.png", f"images/{NAME}.svg"] + [f"images/source_data_{NAME}_{s}.csv" for s in ("A_histogram", "B_histogram", "C_density", "C_fit")]


def test_run_writes_the_three_files_and_is_byte_identical_across_runs(tmp_path):
    run_cli(SUBSET, tmp_path / "a")
    run_cli(SUBSET, tmp_path / "b")
    for name in FILES:
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name
    stats = json.loads((tmp_path / "a" / "statistics.json").read_text())
    assert stats["graph"]["sha256"] == hashlib.sha256(SUBSET.read_bytes()).hexdigest()
    assert stats["graph"]["file"] == "graph_2021_subset.json"


def test_run_is_byte_identical_across_runs_on_the_bundled_release_graph(tmp_path):
    run_cli(BUNDLED, tmp_path / "a")
    run_cli(BUNDLED, tmp_path / "b")
    for name in FILES:
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name
    stats = json.loads((tmp_path / "a" / "statistics.json").read_text())
    pin = json.loads((ROOT / "chamredb" / "data" / "graph_release.json").read_text())
    assert (stats["graph"]["nodes"], stats["graph"]["edges"]) == (pin["nodes"], pin["edges"])
