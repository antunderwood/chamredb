"""Parsers and checks fail loudly on broken input (ticket #4).

Level 1: tiny broken files made by hand; each case names the exact problem a user would see. The first
cases run the real parser scripts as subprocesses; they need the builder environment.
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

BUILDER = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BUILDER / "scripts" / "functions"))
import validation as v  # noqa: E402

NCBI_HEADER = ["allele", "gene_family", "product_name", "subclass", "refseq_protein_accession",
               "refseq_nucleotide_accession", "scope", "type", "subtype"]


def write_ncbi(path, header, rows):
    path.write_text("\n".join("\t".join(r) for r in [header] + rows) + "\n")


def run(script, *args):
    return subprocess.run([sys.executable, str(BUILDER / "scripts" / script), *map(str, args)],
                          capture_output=True, text=True)


# ---- pure checks ----

def test_require_columns_names_the_missing_ones():
    with pytest.raises(ValueError, match=r"x\.tsv.*\['b'\]"):
        v.require_columns(["a", "c"], ["a", "b"], "x.tsv")
    v.require_columns(["a", "b"], ["a", "b"], "x.tsv")


def test_check_count_bounds():
    v.check_count(5, {"min": 5, "max": 5}, "x")
    v.check_count(5, None, "x")
    for n in (4, 6):
        with pytest.raises(ValueError, match="outside the stated bounds"):
            v.check_count(n, {"min": 5, "max": 5}, "x")


def test_check_unique():
    v.check_unique(["a", "b"], "x")
    with pytest.raises(ValueError, match=r"1 duplicate ids, e\.g\. \['a'\]"):
        v.check_unique(["a", "b", "a"], "x")


def manifest(card_release="3.1.4", obo_release="3.1.4"):
    entry = lambda db, role, rel, path: {"database": db, "role": role, "url": "u", "release": rel,
                                         "retrieved": "2026-11-09", "sha256": "0" * 64, "path": path}
    card = lambda role, rel, path: {**entry("card", role, rel, path), "card_release": rel}
    return {"files": [card("sequences", card_release, "a"), card("ontology", obo_release, "b"),
                      entry("ncbi", "catalogue", "2021-09-30.1", "c")]}


def test_manifest_requires_one_card_release():
    v.check_manifest(manifest())
    with pytest.raises(ValueError, match="different releases"):
        v.check_manifest(manifest(obo_release="4.0.2"))


def test_manifest_requires_card_release_declared():
    m = manifest()
    del m["files"][1]["card_release"]
    with pytest.raises(ValueError, match="declare card_release"):
        v.check_manifest(m)


def test_manifest_requires_fields_and_matching_hashes(tmp_path):
    m = manifest()
    m["files"][2]["sha256"] = ""
    with pytest.raises(ValueError, match="missing sha256"):
        v.check_manifest(m)
    (tmp_path / "f").write_text("hello")
    good = {"files": [{"database": "card", "role": "r", "url": "u", "release": "1", "retrieved": "d",
                       "card_release": "1", "path": "f", "sha256": v.sha256_file(tmp_path / "f")}]}
    v.check_manifest(good, tmp_path)
    (tmp_path / "f").write_text("changed")
    with pytest.raises(ValueError, match="sha256 differs"):
        v.check_manifest(good, tmp_path)


# ---- the parser scripts accept broken fixtures today (red), and must refuse them ----

def test_ncbi_protein_parser_refuses_renamed_column(tmp_path):
    header = [("refseq_protein_acc" if c == "refseq_protein_accession" else c) for c in NCBI_HEADER]
    meta = tmp_path / "m.tsv"
    write_ncbi(meta, header, [["a1", "g", "p", "S", "WP_1", "NG_1", "core", "AMR", "AMR"]])
    fa = tmp_path / "p.fa"
    fa.write_text(">x|WP_1|x|x|a1\nMKV\n")
    r = run("parse_ncbi_proteins.py", "--meta", meta, "--input", fa, "--output", tmp_path / "o.fa")
    assert r.returncode != 0 and "missing required columns" in r.stderr and "refseq_protein_accession" in r.stderr


def test_ncbi_metadata_parser_refuses_renamed_column(tmp_path):
    header = [("product_name" if c == "product_name" else c) for c in NCBI_HEADER if c != "subclass"]
    meta = tmp_path / "m.tsv"
    write_ncbi(meta, header, [["a1", "g", "p", "WP_1", "NG_1", "core", "AMR", "AMR"]])
    r = run("parse_ncbi_metadata.py", "--meta", meta, "--output", tmp_path / "o.json")
    assert r.returncode != 0 and "missing required columns" in r.stderr and "subclass" in r.stderr


def test_ncbi_protein_parser_refuses_count_outside_bounds(tmp_path):
    meta = tmp_path / "m.tsv"
    write_ncbi(meta, NCBI_HEADER, [["a1", "g", "p", "S", "WP_1", "NG_1", "core", "AMR", "AMR"]])
    fa = tmp_path / "p.fa"
    fa.write_text(">x|WP_1|x|x|a1\nMKV\n")
    bounds = tmp_path / "b.json"
    bounds.write_text(json.dumps({"ncbi": {"min": 2, "max": 10}}))
    r = run("parse_ncbi_proteins.py", "--meta", meta, "--input", fa, "--output", tmp_path / "o.fa", "--bounds", bounds)
    assert r.returncode != 0 and "outside the stated bounds" in r.stderr
    bounds.write_text(json.dumps({"ncbi": {"min": 1, "max": 10}}))
    assert run("parse_ncbi_proteins.py", "--meta", meta, "--input", fa, "--output", tmp_path / "o.fa", "--bounds", bounds).returncode == 0


def test_card_protein_parser_refuses_duplicate_ids_and_bad_headers(tmp_path):
    fa = tmp_path / "c.fa"
    fa.write_text(">gb|X1|ARO:1|a\nMK\n>gb|X2|ARO:1|b\nMK\n")
    r = run("parse_card_proteins.py", "--input", fa, "--output", tmp_path / "o.fa")
    assert r.returncode != 0 and "duplicate ids" in r.stderr
    fa.write_text(">no_pipes_here\nMK\n")
    r = run("parse_card_proteins.py", "--input", fa, "--output", tmp_path / "o.fa")
    assert r.returncode != 0 and "header" in r.stderr


def test_resfinder_metadata_parser_refuses_missing_column(tmp_path):
    meta = tmp_path / "r.tsv"
    meta.write_text("Gene_accession no.\tClass\nblaX_1_A\tbeta\n")
    r = run("parse_resfinder_metadata.py", "--meta", meta, "--output", tmp_path / "o.json")
    assert r.returncode != 0 and "missing required columns" in r.stderr and "Phenotype" in r.stderr


def test_write_manifest_records_every_input_with_its_hash(tmp_path):
    data = tmp_path / "data"
    names = ["protein_fasta_protein_homolog_model.fasta", "db_metadata/card.metadata.obo",
             "db_metadata/ncbi.metadata.tsv", "db_fastas/ncbi.protein.raw.fasta", "db_fastas/ncbi.nucl.raw.fasta",
             "db_metadata/resfinder.metadata.tsv", "db_fastas/resfinder.nucl.fasta"]
    for n in names:
        (data / n).parent.mkdir(parents=True, exist_ok=True)
        (data / n).write_text(n)
    inputs = tmp_path / "inputs.json"
    inputs.write_text(json.dumps({
        "card": {"release": "3.1.4", "url": "u1"}, "ontology": {"release": "3.1.4", "url": "u2"},
        "ncbi": {"release": "2021-09-30.1", "base_url": "u3/"}, "resfinder": {"git_url": "u4", "commit": "abc"}}))
    out = tmp_path / "m.json"
    r = run("write_manifest.py", "--inputs", inputs, "--data", data, "--out", out, "--retrieved", "2026-11-09")
    assert r.returncode == 0, r.stderr
    files = json.loads(out.read_text())["files"]
    assert sorted(f["path"] for f in files) == sorted(names)
    assert all(f["retrieved"] == "2026-11-09" and f["sha256"] == v.sha256_file(data / f["path"]) for f in files)
    (data / "db_metadata/card.metadata.obo").unlink()
    assert run("write_manifest.py", "--inputs", inputs, "--data", data, "--out", out).returncode != 0


# ---- ticket #12: CARD release agreement is checked against the URLs, not only declared ----

def test_manifest_rejects_card_files_whose_urls_name_different_releases():
    m = manifest(card_release="4.0.2", obo_release="4.0.2")
    m["files"][0]["url"] = "https://card.mcmaster.ca/download/0/broadstreet-v4.0.2.tar.bz2"
    m["files"][1]["url"] = "https://card.mcmaster.ca/download/5/ontology-v3.1.4.tar.bz2"
    with pytest.raises(ValueError, match="URL.*3.1.4.*4.0.2"):
        v.check_manifest(m)
    m["files"][1]["url"] = "https://card.mcmaster.ca/download/5/ontology-v4.0.2.tar.bz2"
    v.check_manifest(m)
