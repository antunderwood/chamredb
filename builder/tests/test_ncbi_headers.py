"""AMRProt FASTA header layouts (ticket #13).

2021 (11 fields):  gi|accession|1|1|allele|gene_family|...   e.g. 1028100430|WP_063856565.1|1|1|vanB|vanB|...
2026 (10 fields):  accession|1|1|allele|gene_family|...      e.g. WP_032675805.1|1|1|blaMIR-9|blaMIR|...
The unmodified parser read field 1 as the accession, so on the 2026 file it matched nothing and wrote an empty FASTA.
Level 1: headers copied from the real files by hand.
"""
import subprocess
import sys
from pathlib import Path

import pytest

BUILDER = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BUILDER / "scripts" / "functions"))
import ncbi_headers as nh  # noqa: E402

OLD = "1028100430|WP_063856565.1|1|1|vanB|vanB|re-structuring_of_the_cell_wall|2|VANCOMYCIN|GLYCOPEPTIDE|D-alanine--(R)-lactate_ligase_VanB"
NEW = "WP_032675805.1|1|1|blaMIR-9|blaMIR|hydrolase|2|CEPHALOSPORIN|BETA-LACTAM|cephalosporin-hydrolyzing_class_C_beta-lactamase_MIR-9"


def test_both_layouts_give_accession_and_allele():
    assert nh.parse_header(OLD) == ("WP_063856565.1", "vanB")
    assert nh.parse_header(NEW) == ("WP_032675805.1", "blaMIR-9")


def test_unknown_layout_is_an_error():
    with pytest.raises(ValueError, match="unexpected AMRProt header"):
        nh.parse_header("a|b|c")


def test_protein_parser_reads_the_2026_layout(tmp_path):
    meta = tmp_path / "m.tsv"
    cols = ["allele", "gene_family", "product_name", "subclass", "refseq_protein_accession",
            "refseq_nucleotide_accession", "scope", "type", "subtype"]
    meta.write_text("\t".join(cols) + "\n" + "\t".join(["blaMIR-9", "blaMIR", "p", "CEPHALOSPORIN", "WP_032675805.1", "NG_1", "core", "AMR", "AMR"]) + "\n")
    fa = tmp_path / "p.fa"
    fa.write_text(f">{NEW}\nMKV\n")
    out = tmp_path / "o.fa"
    r = subprocess.run([sys.executable, str(BUILDER / "scripts" / "parse_ncbi_proteins.py"), "--meta", str(meta),
                        "--input", str(fa), "--output", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert out.read_text() == ">WP_032675805.1\nMKV\n"
