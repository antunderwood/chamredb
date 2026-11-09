"""Name cleaning and similarity invariants (ticket #37, spec #34 invariants 1 to 4).

Level 1: hand-worked cleaning cases. Level 4: the normalised Levenshtein similarity checked against an independent
pure-Python dynamic-programming implementation. Level 5: invariants on the hand graph, the 2021 subset graph and the
bundled release graph.
"""
from pathlib import Path

import pytest

from analysis import matches
from analysis.tests.test_matches import GRAPH as HAND_GRAPH

ROOT = Path(__file__).resolve().parents[2]
SUBSET = ROOT / "chamredb" / "tests" / "test_data" / "graph_2021_subset.json"
BUNDLED = ROOT / "chamredb" / "data" / "graph.json"

# (name, database, remove_alleles, expected), each worked out by hand from the rules in clean_name
CASES = [
    # CARD species pattern one: "<Genus> <species> <gene> conferring resistance to <drug>" keeps the gene
    ("Staphylococcus aureus mupB conferring resistance to mupirocin", "card", False, "mupb"),
    # CARD species pattern two: "<gene> from <Genus> <species>" keeps the gene (real CARD name; red before #37)
    ("dfrA6 from Proteus mirabilis", "card", False, "dfra6"),
    # CARD species pattern three: "<Genus> <species> <gene>" keeps the gene
    ("Clostridium butyricum catB", "card", False, "catb"),
    # species names are removed for CARD only
    ("Staphylococcus aureus mupB conferring resistance to mupirocin", "ncbi", False,
     "staphylococcus aureus mupb conferring resistance to mupirocin"),
    # lower case, bla prefix, brackets, hyphens
    ("blaTEM-1", "ncbi", False, "tem1"),
    ("tet(37)", "card", False, "tet37"),
    ("aac(6')-Ib", "ncbi", False, "aac6'ib"),   # the apostrophe is kept
    # alleles removed: trailing digits and what follows the first digit run
    ("qnrB20", "ncbi", True, "qnrb"),
    ("blaTEM-1", "resfinder", True, "tem"),
    ("OKP-B-6", "card", True, "okpb"),           # okpb6 -> okpb
    ("aac(6')-Ib", "ncbi", True, "aac6'ib"),     # the apostrophe stops the allele pattern: unchanged
    ("vanA", "card", True, "vana"),              # no digits: unchanged
]


@pytest.mark.parametrize("name,database,remove_alleles,expected", CASES)
def test_clean_name_hand_cases(name, database, remove_alleles, expected):
    assert matches.clean_name(name, database, remove_alleles=remove_alleles) == expected


@pytest.mark.xfail(strict=True, reason=(
    "FINDING (ticket #37): clean_name removes 'bla' wherever it occurs, not only as a prefix (the README says prefix). "
    "6 of 15,859 names in the 2026-11-09 graph are affected, for example cblA (lower case cbla) becomes c. Not fixed here: "
    "it changes the numbers, so it needs your decision (comment on #36)."))
def test_bla_is_removed_only_as_a_prefix():
    assert matches.clean_name("cblA", "ncbi") == "cbla"


def levenshtein(a, b):
    """Independent implementation: two-row dynamic programming."""
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        for j, cb in enumerate(b, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (ca != cb)))
        previous = current
    return previous[-1]


def reference_similarity(a, b):
    longest = max(len(a), len(b))
    return 1.0 if longest == 0 else 1 - levenshtein(a, b) / longest


PAIRS = [("qnrb2", "qnrb20"), ("vana", "vanhax"), ("tem1", "tem1"), ("", "a"), ("", ""), ("kitten", "sitting"),
         ("flaw", "lawn"), ("a", "b"), ("abc", "")]


@pytest.mark.parametrize("a,b", PAIRS)
def test_similarity_equals_one_minus_distance_over_longest_length(a, b):
    assert matches._levenshtein.similarity(a, b) == pytest.approx(reference_similarity(a, b), abs=1e-12)


def test_similarity_matches_the_reference_on_every_cleaned_name_pair_of_the_subset_graph():
    table, _ = matches.best_matches(matches.read_graph(SUBSET))
    for row in table.itertuples():
        assert row.name_similarity == pytest.approx(
            reference_similarity(row.cleaned_source_name, row.cleaned_target_name), abs=1e-12)


def pairs_with_an_edge(graph):
    """Independent set logic: (node, other database) pairs that have at least one out-edge."""
    return {(u, v.split(":")[0]) for u, v in graph.edges()}


def check_invariants(graph):
    table, skipped = matches.best_matches(graph)
    assert len(table) + skipped == len(pairs_with_an_edge(graph))                      # invariant 1
    assert ((table.sequence_identity > 0) & (table.sequence_identity <= 1)).all()    # invariant 2
    assert ((table.name_similarity >= 0) & (table.name_similarity <= 1)).all()


def test_invariants_on_the_hand_graph():
    check_invariants(HAND_GRAPH)


def test_invariants_on_the_2021_subset_graph():
    check_invariants(matches.read_graph(SUBSET))


def test_invariants_on_the_bundled_release_graph():
    check_invariants(matches.read_graph(BUNDLED))
