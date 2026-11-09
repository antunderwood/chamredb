"""Characterisation harness for the builder (ticket #9).

Builds the graph from the tiny fixture end to end and compares it with the frozen expected graph.
Expected values prove "unchanged", not "correct"; the designed cases below add hand-checkable
expectations derived from the fixture sequences. make_fixtures.py reads the shipped 2021 graph only
to choose decoys; nothing in these tests does.

Not covered here (see README): the coverage invariant 0 < coverage <= 1 as a general property
(ticket #3), per-query sort order, reverse-strand ResFinder genes, the download rules.
"""
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "fixtures" / "data"
EXPECTED = HERE / "expected" / "fixture_graph.json"
MMSEQS2_VERSION = "18.8cc5c"
TOLERANCE = 1e-3


@pytest.fixture(scope="session")
def built(tmp_path_factory):
    if os.environ.get("FIXTURE_BUILD_DIR"):  # reuse a finished build; may be stale, see README
        return Path(os.environ["FIXTURE_BUILD_DIR"])
    if shutil.which("mmseqs") is None:
        pytest.fail(f"needs MMseqs2 {MMSEQS2_VERSION} on Linux amd64 (use `make test-docker` on macOS)")
    work = tmp_path_factory.mktemp("fixture_build") / "work"
    subprocess.run(["bash", str(HERE / "run_fixture_build.sh"), str(work)], check=True)
    return work


@pytest.fixture(scope="session")
def graph(built):
    return json.loads((built / "graph.json").read_text())


def read_versions(text):
    out = {}
    for line in text.splitlines():
        key, _, value = line.partition(" ")
        if key in {"mmseqs2", "snakemake", "python", "pandas", "networkx", "Bio", "pronto", "orfipy"}:
            out[key] = ".".join(value.split(".")[:2]) if key == "python" else value
    return out


def test_tool_versions_match_frozen_set(built):
    got = read_versions((built / "versions.txt").read_text())
    want = read_versions((HERE / "expected" / "VERSIONS.txt").read_text())
    assert want["mmseqs2"] == MMSEQS2_VERSION
    assert got == want, "tool versions differ from the set the expected graph was frozen with"


def test_fixture_files_unchanged():
    lines = (HERE / "expected" / "SHA256SUMS").read_text().splitlines()
    want = {name: digest for digest, name in (line.split("  ", 1) for line in lines)}
    got = {str(p.relative_to(FIXTURE)): hashlib.sha256(p.read_bytes()).hexdigest()
           for p in sorted(FIXTURE.rglob("*")) if p.is_file()}
    assert got == want


def test_matches_frozen_graph(graph):
    expected = json.loads(EXPECTED.read_text())
    assert len(graph["nodes"]) == len(expected["nodes"])
    assert len(graph["links"]) == len(expected["links"])
    assert ({n["id"]: n for n in graph["nodes"]}
            == {n["id"]: n for n in expected["nodes"]})
    got = {(e["source"], e["target"]): e for e in graph["links"]}
    want = {(e["source"], e["target"]): e for e in expected["links"]}
    assert got.keys() == want.keys()
    for key, edge in want.items():
        assert got[key]["type"] == edge["type"], key
        assert got[key]["identity"] == pytest.approx(edge["identity"], abs=TOLERANCE), key
        assert got[key]["coverage"] == pytest.approx(edge["coverage"], abs=TOLERANCE), key


# ---- designed cases: hand-checkable expectations from how the fixture was constructed ----

CARD = lambda i: f"card:ARO:{i}"
NCBI = lambda i: f"ncbi:{i}"
RES = lambda i: f"resfinder:{i}"
TRIPLES = [  # (card, ncbi, resfinder) mutual best hits in all three databases
    (CARD("3002864"), NCBI("WP_063844476.1"), RES("dfrB1")),
    (CARD("3003021"), NCBI("WP_063844477.1"), RES("dfrB2")),
    (CARD("3003022"), NCBI("WP_010890144.1"), RES("dfrB3")),
    (CARD("3003023"), NCBI("WP_063844479.1"), RES("dfrB6")),
    (CARD("3002871"), NCBI("WP_063855713.1"), RES("tet(37)")),
    (CARD("3004556"), NCBI("WP_063844480.1"), RES("dfrB8")),  # named dfrB7 in CARD, dfrB8 elsewhere
]
CARD_NCBI_ONLY = [
    (CARD("3000124"), NCBI("WP_000369216.1")),  # mecI
    (CARD("3004253"), NCBI("WP_021418929.1")),  # vanU-G
]
PARTIAL = (CARD("3002718"), NCBI("WP_017900990.1"))  # NCBI protein is the CARD protein minus residue 1


def edges(graph):
    return {(e["source"], e["target"]): e for e in graph["links"]}


def types(graph):
    return {k: e["type"] for k, e in edges(graph).items()}


def test_mutual_best_hit_triples(graph):
    got = types(graph)
    for trio in TRIPLES:
        for a in trio:
            for b in trio:
                if a != b:
                    assert got.get((a, b)) == "RBH", (a, b)


def test_card_ncbi_only_pairs(graph):
    got = types(graph)
    for card, ncbi in CARD_NCBI_ONLY:
        assert got[(card, ncbi)] == got[(ncbi, card)] == "RBH"
        assert not [1 for (s, t) in got
                    if (s in (card, ncbi) or t in (card, ncbi)) and "resfinder" in (s.split(":")[0], t.split(":")[0])]


def test_dfrB9_pair_has_only_one_way_hits_to_resfinder(graph):
    got = types(graph)
    card, ncbi = CARD("3005350"), NCBI("WP_071846200.1")
    assert got[(card, ncbi)] == got[(ncbi, card)] == "RBH"
    for src in (card, ncbi):
        to_res = {t: ty for (s, t), ty in got.items() if s == src and t.startswith("resfinder")}
        assert to_res and set(to_res.values()) == {"OWH"}


def test_resfinder_dfrB7_has_no_reciprocal_best_hit(graph):
    got = types(graph)
    out = {t: ty for (s, t), ty in got.items() if s == RES("dfrB7")}
    assert out and set(out.values()) == {"OWH"}
    # others may hit dfrB7 one way (dfrB9 does), but no RBH touches it in either direction
    assert not any(ty == "RBH" and RES("dfrB7") in (s, t) for (s, t), ty in got.items())


def fasta_by(path, field):
    out, key = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:].split("|")[field]
            out[key] = ""
        else:
            out[key] += line.strip()
    return out


def test_coverage_is_alignment_length_over_query_end(graph):
    # CARD protein = residue 1 + the NCBI protein, so the ungapped alignment covers query positions 2..215:
    # alnlen = 214, qend = 215, coverage = 214/215 (not 1.0). RBH edges carry the first query's values.
    card = fasta_by(FIXTURE / "protein_fasta_protein_homolog_model.fasta", 2)["ARO:3002718"]
    ncbi = fasta_by(FIXTURE / "db_fastas" / "ncbi.protein.raw.fasta", 1)["WP_017900990.1"]
    assert card.index(ncbi) == 1 and len(card) == len(ncbi) + 1  # precondition, from the sequences alone
    expected = len(ncbi) / (card.index(ncbi) + len(ncbi))
    e = edges(graph)
    for key in ((PARTIAL[0], PARTIAL[1]), (PARTIAL[1], PARTIAL[0])):
        assert e[key]["type"] == "RBH"
        assert e[key]["identity"] == pytest.approx(1.0, abs=TOLERANCE)
        assert e[key]["coverage"] == pytest.approx(expected, abs=TOLERANCE)
    assert expected < 0.999


def test_rbh_edges_are_reciprocal(graph):
    got = types(graph)
    for (s, t), ty in got.items():
        if ty == "RBH":
            assert got.get((t, s)) == "RBH", (s, t)


def test_no_edge_within_one_database(graph):
    assert all(e["source"].split(":")[0] != e["target"].split(":")[0] for e in graph["links"])


def test_only_designed_nodes_plus_documented_extras(graph):
    designed = {n for trio in TRIPLES for n in trio} | {n for pair in CARD_NCBI_ONLY for n in pair} | set(PARTIAL) | {
        CARD("3005350"), NCBI("WP_071846200.1"), RES("dfrB7")}
    # extras, not designed: card ARO:3001205 and ncbi WP_004201167.1 (ble, a survivor of a four-accession
    # name group) have identical sequences and form an RBH pair; card ARO:3005036 reaches that NCBI protein
    # at 0.616 identity as a one-way hit. Which accession survives de-duplication decides this pair.
    extra = {CARD("3001205"), NCBI("WP_004201167.1"), CARD("3005036")}
    assert {n["id"] for n in graph["nodes"]} == designed | extra


# ---- parser behaviour and invariants (spec predicates 4 to 6, checked on the built intermediate files) ----

def built_ids(built, db):
    path = built / "data" / "db_fastas" / f"{db}.protein.fasta"
    return [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]


def ncbi_rows():
    with open(FIXTURE / "db_metadata" / "ncbi.metadata.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_protein_ids_unique_per_database(built):
    for db in ("card", "ncbi", "resfinder"):
        ids = built_ids(built, db)
        assert ids and len(ids) == len(set(ids)), db


def test_graph_nodes_are_parsed_proteins(built, graph):
    ids = {db: set(built_ids(built, db)) for db in ("card", "ncbi", "resfinder")}
    for n in graph["nodes"]:
        db, ident = n["id"].split(":", 1)
        assert ident in ids[db], n["id"]


def test_ncbi_proteins_equal_unique_names_among_core_amr(built):
    # computed from the raw fixture table, not with the parser's code path: core AMR rows, one protein per name
    raw = fasta_by(FIXTURE / "db_fastas" / "ncbi.protein.raw.fasta", 1)
    names = {r["allele"] or r["gene_family"] for r in ncbi_rows()
             if r["scope"] == "core" and r["type"] == "AMR" and r["subtype"] == "AMR"
             and r["refseq_protein_accession"] in raw}
    assert len(built_ids(built, "ncbi")) == len(names)


def test_non_core_proteins_dropped_and_one_protein_kept_per_shared_name(built):
    prot = set(built_ids(built, "ncbi"))
    rows = ncbi_rows()
    plus = {r["refseq_protein_accession"] for r in rows if r["scope"] != "core"}
    assert plus and not (plus & prot)
    groups = {}
    for r in rows:
        if r["scope"] == "core" and not r["allele"]:
            groups.setdefault(r["gene_family"], []).append(r["refseq_protein_accession"])
    shared = [accs for accs in groups.values() if len(accs) > 1]
    assert shared
    for accs in shared:
        assert len(set(accs) & prot) == 1
