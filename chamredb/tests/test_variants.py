"""The bundled variants.json, the truth set, mutation parsing and the `chamredb variants` CLI."""
import csv
import json
import os
import subprocess
import sys

import pytest

from chamredb import graph_loader
from chamredb.query_variants import parse_hamronization_mutations
from chamredb.variants import Variants, normalise_mutation, resolve_variants, split_query

HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = os.path.join(HERE, "test_data", "variant_truth_set.tsv")


@pytest.fixture(scope="module")
def variants():
    return Variants.load()


def _truth_rows():
    with open(TRUTH) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def test_bundled_variants_match_their_pin():
    path, release = resolve_variants()
    pin = json.loads((graph_loader.PACKAGE_DATA / "graph_release.json").read_text())["variants"]
    assert release == pin["data_release"]
    data = json.loads(path.read_text())
    assert len(data["records"]) == pin["records"]


def test_truth_set_covers_every_class():
    classes = {r["expected_class"] for r in _truth_rows()}
    single = {"same_mutation", "same_mutation_different_reference", "possible_same_mutation_offset",
              "same_mutation_only_in_combination", "same_position_different_change", "different_mutation",
              "same_mutation_related_organism", "locus_not_covered", "organism_not_covered"}
    combination = {"same_combination", "partial_combination", "components_listed_separately", "combination_not_listed"}
    assert single | combination <= classes, (single | combination) - classes


@pytest.mark.parametrize("row", _truth_rows(), ids=lambda r: f"{r['database']}:{r['query']}->{r['other_database']}")
def test_truth_set(variants, row):
    found = variants.find(row["query"], database=row["database"], taxon=row["taxon"])
    assert found, f"{row['query']} not found in {row['database']} for {row['taxon']}"
    results = [variants.comparison(r)[row["other_database"]] for r in found]
    classes = [m["class"] for m in results]
    assert row["expected_class"] in classes, f"{classes} ({row['rationale']})"
    if row["expected_match"]:
        m = results[classes.index(row["expected_class"])]
        notations = [h["native_notation"] for h in m["matched"]]
        assert any(row["expected_match"] in n for n in notations), notations


@pytest.mark.parametrize("text,expected", [
    ("S83L", "S83L"), ("p.Ser83Leu", "S83L"), ("p.S83L", "S83L"), ("Q10Ter", "Q10*"), ("Q10*", "Q10*"),
    ("g2032t", "G2032T"), ("n.2059A>G", "A2059G"), ("c.-15C>T", "C-15T"), ("-15c>t", "C-15T"), ("A2075G", "A2075G"),
    ("SER83LEU", "S83L"), ("p.Leu197=", "L197L"), ("R447Var", "R447?"), ("G129fs", "G129fs"),
    # indels: the residues beside the positions are optional, three-letter residues are read anywhere
    ("p.Asn437_Asn438del", "437_438del"), ("N437_N438del", "437_438del"), ("437_438del", "437_438del"),
    ("p.Arg177_Thr178insIleGly", "177_178insIG"), ("D126del", "126del"), ("c.983del", "983del"),
    ("p.Asn474delinsLysGly", "474delinsKG"), ("474_474delinsKG", "474delinsKG"), ("-14_-13insGT", "-14_-13insGT"),
    ("131_132insTERY", "131_132insTERY"),   # one-letter residues that spell Ter are not a stop
])
def test_normalise_mutation(text, expected):
    assert normalise_mutation(text) == expected


@pytest.mark.parametrize("query,expected", [
    ("gyrA_S83L", ("gyrA", "S83L")), ("gyrA S83L", ("gyrA", "S83L")), ("gyrA:p.Ser83Leu", ("gyrA", "S83L")),
    ("23S A2059G", ("23S", "A2059G")), ("S83L", (None, "S83L")), ("ampC_C-42T", ("ampC", "C-42T")),
    ("16S_rrsB_C1192T", ("16S_rrsB", "C1192T")), ("Rv0678_A110fs", ("Rv0678", "A110fs")), ("pbp5 467_468insD", ("pbp5", "467_468insD")),
    # a range keeps both ends
    ("D63_S67del", (None, "63_67del")), ("gyrA_D63_S67del", ("gyrA", "63_67del")), ("rpoB_427_429del", ("rpoB", "427_429del")),
    ("ampC_-14_-13insGT", ("ampC", "-14_-13insGT")), ("liaF:p.Arg177_Thr178insIle", ("liaF", "177_178insI")),
])
def test_split_query(query, expected):
    assert split_query(query) == expected


def test_find_accepts_hgvs_and_native_forms(variants):
    a = {r["record_id"] for r in variants.find("gyrA p.Ser83Leu", database="ncbi", taxon="Escherichia")}
    b = {r["record_id"] for r in variants.find("gyrA_S83L", database="ncbi", taxon="Escherichia")}
    assert a and a == b


def test_every_single_mutation_is_found_by_its_own_notations(variants):
    """Native notation, gene with compact form (space or underscore) and gene with HGVS all find the record."""
    missed = []
    for r in variants.records:
        if len(r["components"]) != 1:
            continue
        c = r["components"][0]
        gene = c["gene_label"]
        for query in (r["native_notation"], f"{gene} {c['native_compact']}", f"{gene}_{c['native_compact']}",
                      f"{gene} {c['native_hgvs']}", f"{gene}:{c['native_hgvs']}"):
            if r not in variants.find(query, database=r["database"]):
                missed.append((r["database"], query))
    assert not missed, missed[:10]


def test_locus_lookup(variants):
    assert variants.locus_ids("Escherichia:gyrA") == ["Escherichia:gyrA:protein"]
    assert variants.locus_ids("Neisseria gonorrhoeae:23S") == ["Neisseria_gonorrhoeae:23S rRNA:rna"]


def _cli(*args):
    env = {**os.environ, "COLUMNS": "200", "PYTHONIOENCODING": "utf-8"}
    return subprocess.run([sys.executable, "-m", "chamredb.run_chamredb", "variants", *args],
                          capture_output=True, text=True, env=env)


def test_cli_single_mutation():
    out = _cli("-d", "ncbi", "-i", "gyrA_S83L", "-t", "Escherichia")
    assert out.returncode == 0, out.stderr
    for expected in ("data release: 2026-11-19", "p.Ser83Leu", "same mutation, different reference sequence", "CARD", "PointFinder"):
        assert expected in out.stdout, expected


def test_cli_unknown_mutation_exits_1():
    out = _cli("-d", "ncbi", "-i", "gyrA_W999Y", "-t", "Escherichia")
    assert out.returncode == 1 and "no mutation matching" in out.stdout


def test_cli_locus_table():
    out = _cli("-l", "Escherichia:gyrA")
    assert out.returncode == 0, out.stderr
    assert "p.Ser83Leu" in out.stdout and "gyrA_S83L" in out.stdout


def test_cli_id_file(tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("gyrA_S83L\nparC_S80I\nnot_a_mutation\n")
    out_tsv = tmp_path / "out.tsv"
    out = _cli("-d", "ncbi", "-t", "Escherichia", "-f", str(ids), "-o", str(out_tsv))
    assert out.returncode == 0, out.stderr
    rows = list(csv.DictReader(open(out_tsv), delimiter="\t"))
    assert [r["found"] for r in rows] == ["yes", "yes", "no"]
    assert rows[0]["card_status"] == "same_mutation_different_reference"
    assert rows[1]["pointfinder_status"] == "same_mutation_only_in_combination"
    assert all(r["data_release"] == "2026-11-19" for r in rows)


HAMRONIZE = [
    {"input_file_name": "s1", "gene_symbol": "gyrA_S83L", "analysis_software_name": "amrfinderplus",
     "reference_accession": "WP_001281243.1", "genetic_variation_type": "protein_variant_detected"},
    {"input_file_name": "s2", "gene_symbol": "gyrA", "analysis_software_name": "rgi", "reference_accession": "3003294",
     "amino_acid_mutation": "S83L,D87N", "genetic_variation_type": "protein_variant_detected"},
    {"input_file_name": "s3", "gene_symbol": "blaTEM-1", "analysis_software_name": "amrfinderplus",
     "reference_accession": "WP_000027057.1", "genetic_variation_type": "gene_presence_detected"},
]


def test_parse_hamronization_mutations(tmp_path):
    p = tmp_path / "h.json"
    p.write_text(json.dumps(HAMRONIZE))
    items = parse_hamronization_mutations(p)
    assert [(i["database"], i["id"]) for i in items] == [("ncbi", "gyrA S83L"), ("card", "gyrA S83L,gyrA D87N")]


def test_cli_hamronization(tmp_path):
    p = tmp_path / "h.json"
    p.write_text(json.dumps(HAMRONIZE))
    out_tsv = tmp_path / "out.tsv"
    out = _cli("-j", str(p), "-t", "Escherichia", "-o", str(out_tsv))
    assert out.returncode == 0, out.stderr
    rows = list(csv.DictReader(open(out_tsv), delimiter="\t"))
    assert rows[0]["input_file"] == "s1" and rows[0]["ncbi_notation"] == "gyrA_S83L"
    combo = [r for r in rows if r["input_file"] == "s2"]
    assert combo and combo[0]["pointfinder_status"] == "same_combination"
