"""Mutation notation, with examples taken from CARD 4.0.2, NCBI 2026-10-29.1 and pointfinder_db 3bd7d5e,
and the alignment position map.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from functions import liftover  # noqa: E402
from functions import mutation_notation as mn  # noqa: E402


def c(molecule, kind, start, end, ref="", alt=""):
    return mn.change(molecule, kind, start, end, ref, alt)


@pytest.mark.parametrize("token,molecule,expected", [
    ("S83L", "protein", c("protein", "substitution", 83, 83, "S", "L")),
    ("C52Ter", "protein", c("protein", "nonsense", 52, 52, "C", "*")),
    ("R447Var", "protein", c("protein", "any_substitution", 447, 447, "R", "?")),
    ("F50F", "protein", c("protein", "synonymous", 50, 50, "F", "F")),
    ("G129fs", "protein", c("protein", "frameshift", 129, 129, "G")),
    ("D63_S67del", "protein", c("protein", "deletion", 63, 67)),
    ("D126del", "protein", c("protein", "deletion", 126, 126, "D")),
    ("R97_G98insA", "protein", c("protein", "insertion", 97, 98, "RG", "A")),
    ("P627dup", "protein", c("protein", "duplication", 627, 627, "P")),
    ("V624_P626dup", "protein", c("protein", "duplication", 624, 626)),
    ("g1065c", "protein", c("dna", "substitution", 1065, 1065, "G", "C")),
    ("g2032t", "rna", c("rna", "substitution", 2032, 2032, "G", "T")),
    ("-11a>g", "protein", c("promoter", "substitution", -11, -11, "A", "G")),
])
def test_card(token, molecule, expected):
    assert mn.parse_card(token, molecule) == expected


@pytest.mark.parametrize("token", ["Q10Z", "Q142X", "N330_R335insTIPY"])
def test_card_rejects(token):
    with pytest.raises(mn.NotationError):
        mn.parse_card(token, "protein")


def test_card_codependent():
    assert mn.parse_card_codependent("45612,N87K,V199A+45613,D481E") == [("45612", ["N87K", "V199A"]), ("45613", ["D481E"])]


@pytest.mark.parametrize("symbol,position,molecule,expected", [
    ("S83L", 83, "protein", c("protein", "substitution", 83, 83, "S", "L")),
    # conventional numbering in the symbol, sequence position in mutation_position (E. faecalis parC)
    ("S80I", 82, "protein", c("protein", "substitution", 82, 82, "S", "I")),
    ("Q163Ter", 163, "protein", c("protein", "nonsense", 163, 163, "Q", "*")),
    ("ASLD153del", 153, "protein", c("protein", "deletion", 153, 156, "ASLD")),
    ("G585GTG", 585, "protein", c("protein", "insertion", 585, 586, "G", "TG")),
    ("T314AfsTer16ins2", 314, "protein", c("protein", "frameshift", 314, 314, "T")),
    ("A2074G", 2074, "rna", c("rna", "substitution", 2074, 2074, "A", "G")),
    ("T-35G", 204, "promoter", c("promoter", "substitution", -35, -35, "T", "G")),
])
def test_ncbi(symbol, position, molecule, expected):
    got, index = mn.parse_ncbi(symbol, position, molecule)
    assert got == expected and index == position


def test_pointfinder_alternatives_are_split():
    got = mn.parse_pointfinder("AA", "83", "S", "L,W,A,V", "TCG", False)
    assert [x["alt"] for x in got] == ["L", "W", "A", "V"] and all(x["start"] == 83 for x in got)


@pytest.mark.parametrize("args,expected", [
    (("AA", "103", "K", "*", "AAA", False), c("protein", "nonsense", 103, 103, "K", "*")),
    (("AA", "131", "ins", "SG", "-", False), c("protein", "insertion", 131, 132, "", "SG")),
    (("NUC", "-42", "C", "T", "C", True), c("promoter", "substitution", -42, -42, "C", "T")),
    (("NUC", "-13", "ins", "GT", "-", True), c("promoter", "insertion", -13, -12, "", "GT")),
    (("NUC", "382", "del", "CAGGCC", "-", False), c("dna", "deletion", 382, 387, "CAGGCC")),
])
def test_pointfinder(args, expected):
    assert mn.parse_pointfinder(*args) == [expected]


def test_pointfinder_required():
    assert mn.parse_pointfinder_required("gyrA_1_CP073768.1_S83L.W.A.V") == [[("gyrA_1_CP073768.1", "S", 83, ["L", "W", "A", "V"])]]
    assert mn.parse_pointfinder_required("gyrA_S83Y.F;gyrA_D87N") == [[("gyrA", "S", 83, ["Y", "F"])], [("gyrA", "D", 87, ["N"])]]


def test_promoter_coordinates_have_no_zero():
    assert mn._advance(-1, 1) == 1 and mn._advance(-3, 2) == -1 and mn._advance(5, 2) == 7


@pytest.mark.parametrize("change,text", [
    (c("protein", "substitution", 83, 83, "S", "L"), "p.Ser83Leu"),
    (c("protein", "nonsense", 52, 52, "C", "*"), "p.Cys52Ter"),
    (c("protein", "insertion", 97, 98, "RG", "A"), "p.Arg97_Gly98insAla"),
    (c("rna", "substitution", 2059, 2059, "A", "G"), "n.2059A>G"),
    (c("promoter", "substitution", -15, -15, "C", "T"), "c.-15C>T"),
])
def test_hgvs(change, text):
    assert mn.hgvs(change) == text


def test_hgvs_on_anchor_coordinates():
    assert mn.hgvs(c("protein", "substitution", 85, 85, "S", "L"), 83, 83) == "p.Ser83Leu"


GYRA = "MSDLAREITPVNIEEELKSSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMNVLGNDWNKAYKKSARVVGDVIGKYHPHGDSAVYDTIVRMAQPF"


def test_identical_sequences_map_one_to_one():
    aln = liftover.align(GYRA, GYRA, "protein")
    assert aln["identity"] == 1.0 and aln["map"][83] == 83


def test_extension_shifts_the_map():
    aln = liftover.align(GYRA, "MA" + GYRA, "protein")
    assert aln["map"][85] == 83 and aln["identity"] == 1.0
    assert liftover.lift(aln["map"], 85, 85) == (83, 83)


def test_deletion_in_other_leaves_a_gap():
    other = GYRA[:50] + GYRA[55:]
    aln = liftover.align(GYRA, other, "protein")
    assert aln["map"][51] == 56 and 52 not in aln["map"].values()


def test_nucleotides_and_kmer_filter():
    a = "ACGTACGGTTAGCATCGATCGATTAGCTAGCTAGGCTAGCTAGCATCGAT"
    assert liftover.align(a, a[5:], "nucleotide")["map"][1] == 6
    assert liftover.kmer_similarity(a, a) == 1.0 and liftover.kmer_similarity(a, "TTTTTTTTTT") == 0.0


def test_pointfinder_negative_position_is_a_promoter_change_whatever_the_type():
    # mtrC-promoter-size-174bp CYS--120 is typed AA in pointfinder_db 3bd7d5e: C>T at -120, not p.Cys-120Thr
    assert mn.parse_pointfinder("AA", "-120", "C", "T", "C", True) == [c("promoter", "substitution", -120, -120, "C", "T")]
    assert mn.parse_pointfinder("AA", "-290", "ins", "ATT", "-", True) == [c("promoter", "insertion", -290, -289, "", "ATT")]


def test_pointfinder_any_stop_codon_row_is_not_a_mutation_at_position_1():
    # katG 'ASN-1, 1, N, N, *' under '#Stop codons' lists any premature stop in the gene
    with pytest.raises(mn.NotationError, match="any premature stop codon"):
        mn.parse_pointfinder("AA", "1", "N", "*", "N", False)
    assert mn.parse_pointfinder("AA", "1", "M", "R", "ATG", False) == [c("protein", "substitution", 1, 1, "M", "R")]


def test_hgvs_of_a_one_residue_delins():
    assert mn.hgvs(c("protein", "delins", 474, 474, "N", "KG")) == "p.Asn474delinsLysGly"


def test_ncbi_insertion_index_is_the_residue_before_the_insertion():
    change, index = mn.parse_ncbi("AG585AGTG", 584, "protein")
    assert (change["start"], change["end"], change["ref"], change["alt"], index) == (585, 586, "G", "TG", 585)


def test_reverse_strand_reference_is_aligned_as_its_reverse_complement():
    import random
    rng = random.Random(1)
    anchor = "".join(rng.choice("ACGT") for _ in range(300))
    other = liftover.revcomp(anchor[10:290])
    aln = liftover.align(anchor, other, "nucleotide")
    assert aln["strand"] == -1 and aln["identity"] == 1.0 and aln["coverage"] == 1.0
    assert aln["map"][1] == 290 and aln["map"][280] == 11          # first base of other is the last aligned of anchor
    assert liftover.lift(aln["map"], 5, 7) == (284, 286)             # a range is turned around
    assert liftover.revcomp("AC?G") == "C?GT"
    assert liftover.align(anchor, anchor[10:290], "nucleotide")["strand"] == 1
