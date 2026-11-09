"""ResFinder parsing (ticket #14): duplicate records, reverse-strand genes, allele and phenotype alignment.

Level 1: tiny hand-made inputs whose expected values are worked out by hand in this file. Runs the real scripts.
"""
import json
import subprocess
import sys
from pathlib import Path

BUILDER = Path(__file__).resolve().parents[1]

# hand codon table (one codon per amino acid) so the expected protein is written out by hand, not computed
CODON = {"M": "ATG", "K": "AAA", "V": "GTT", "L": "CTG", "A": "GCG", "G": "GGC", "I": "ATT", "Q": "CAG",
         "P": "CCG", "S": "AGC", "W": "TGG", "T": "ACC", "E": "GAA", "D": "GAT", "R": "CGT"}
PROTEIN = "MKVLAAGIVALLLAQPSWTEDRMKVLAAGIVALLLAQPSWTEDR"  # 44 aa
FORWARD = "".join(CODON[a] for a in PROTEIN) + "TAA"
REVERSE = FORWARD[::-1].translate(str.maketrans("ACGT", "TGCA"))  # reverse complement, written without Bio
OTHER = "".join(CODON[a] for a in PROTEIN[::-1]) + "TAA"


def run(script, *args):
    return subprocess.run([sys.executable, str(BUILDER / "scripts" / script), *map(str, args)], capture_output=True, text=True)


def fasta(path, records):
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


def proteins(path):
    out, key = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:]
            out[key] = ""
        else:
            out[key] += line
    return out


def test_reverse_strand_gene_gives_the_same_protein_as_the_forward_gene(tmp_path):
    nucl = tmp_path / "n.fa"
    fasta(nucl, [("fwd_1_A", FORWARD), ("rev_1_B", REVERSE)])
    r = run("resfinder_orfs_to_proteins.py", "--nucl", nucl, "--protein", tmp_path / "p.fa")
    assert r.returncode == 0, r.stderr
    got = proteins(tmp_path / "p.fa")
    assert got["fwd"].rstrip("*") == PROTEIN
    assert got["rev"].rstrip("*") == PROTEIN


def test_identical_duplicate_records_are_collapsed_and_reported_within_bound(tmp_path):
    nucl = tmp_path / "n.fa"
    fasta(nucl, [("g1_1_A", FORWARD), ("g2_1_B", OTHER), ("g1_1_A", FORWARD)])
    r = run("resfinder_orfs_to_proteins.py", "--nucl", nucl, "--protein", tmp_path / "p.fa")
    assert r.returncode == 0, r.stderr
    assert set(proteins(tmp_path / "p.fa")) == {"g1", "g2"} and "1 duplicate" in r.stderr


def test_doubled_fasta_is_refused(tmp_path):
    # the all.fsa case: every record present twice (cat *.fsa including all.fsa)
    recs = [(f"g{i}_1_A", FORWARD) for i in range(60)]
    nucl = tmp_path / "n.fa"
    fasta(nucl, recs + recs)
    r = run("resfinder_orfs_to_proteins.py", "--nucl", nucl, "--protein", tmp_path / "p.fa")
    assert r.returncode != 0 and "duplicate" in r.stderr


def test_same_id_with_different_sequences_is_refused(tmp_path):
    nucl = tmp_path / "n.fa"
    fasta(nucl, [("g1_1_A", FORWARD), ("g1_1_A", OTHER)])
    r = run("resfinder_orfs_to_proteins.py", "--nucl", nucl, "--protein", tmp_path / "p.fa")
    assert r.returncode != 0 and "different sequences" in r.stderr


PHENO_HEADER = "Gene_accession no.\tClass\tPhenotype\tPMID\tMechanism of resistance\tNotes\tRequired_gene\n"


def pheno(path, rows):
    path.write_text(PHENO_HEADER + "".join(f"{a}\tC\t{p}\t1\tm\t\t\n" for a, p in rows))


def test_metadata_phenotype_comes_from_the_allele_whose_sequence_is_used(tmp_path):
    nucl = tmp_path / "n.fa"
    fasta(nucl, [("gA_1_X", FORWARD), ("gA_2_Y", OTHER), ("gB_1_Z", OTHER)])
    meta = tmp_path / "m.tsv"
    pheno(meta, [("gA_1_X", "Amikacin, Kanamycin"), ("gA_2_Y", "Amikacin"), ("gB_1_Z", "Tobramycin")])
    r = run("parse_resfinder_metadata.py", "--meta", meta, "--nucl", nucl, "--output", tmp_path / "o.json")
    assert r.returncode == 0, r.stderr
    got = json.loads((tmp_path / "o.json").read_text())
    assert got["gA"]["phenotype"] == "confers resistance to Amikacin, Kanamycin"  # allele _1 is the sequence used
    assert got["gB"]["phenotype"] == "confers resistance to Tobramycin"


def test_metadata_missing_for_the_used_allele_is_refused(tmp_path):
    nucl = tmp_path / "n.fa"
    fasta(nucl, [("gA_1_X", FORWARD)])
    meta = tmp_path / "m.tsv"
    pheno(meta, [("gA_2_Y", "Amikacin")])
    r = run("parse_resfinder_metadata.py", "--meta", meta, "--nucl", nucl, "--output", tmp_path / "o.json")
    assert r.returncode != 0 and "gA_1_X" in r.stderr
