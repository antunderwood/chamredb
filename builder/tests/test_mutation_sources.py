"""The three readers on small files laid out like the real ones."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from functions import mutation_sources as ms  # noqa: E402

PROT = "MSDLAREITPVNIEEELKSSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMNVLGNDWNKAYKKSARVVGDVIGKYHPHGDSAVYDTIVRMAQPF"
#       position 83 is S (as in E. coli gyrA S83), 87 is D


def _codons(protein):
    table = {"M": "ATG", "S": "TCG", "D": "GAT", "L": "CTG", "A": "GCG", "R": "CGT", "E": "GAA", "I": "ATT", "T": "ACC",
             "P": "CCG", "V": "GTG", "N": "AAC", "K": "AAA", "Y": "TAT", "G": "GGC", "H": "CAT", "W": "TGG", "Q": "CAG",
             "F": "TTT", "C": "TGC"}
    return "".join(table[a] for a in protein)


def write_card(path):
    model = lambda mid, aro, aid, short, mtype, params, seqs: {
        "model_id": mid, "ARO_accession": aro, "ARO_id": aid, "model_name": short, "CARD_short_name": short,
        "model_type": mtype, "model_param": params, "ARO_category": {"1": {"category_aro_class_name": "Drug Class",
                                                                           "category_aro_name": "fluoroquinolone antibiotic"}},
        "model_sequences": {"sequence": {"1": seqs}}}
    tax = {"NCBI_taxonomy": {"NCBI_taxonomy_name": "Escherichia coli str. K-12"}}
    protein = {"protein_sequence": {"accession": "P1", "sequence": PROT},
               "dna_sequence": {"accession": "D1", "sequence": _codons(PROT)}, **tax}
    card = {
        "_version": "test",
        "1": model("1", "3000001", "100", "Ecol_gyrA_FLO", "protein variant model", {
            "blastp_bit_score": {"param_type": "BLASTP bit-score", "param_value": "500"},
            # no param_type on 'snp' (as in some 4.0.2 models); S83L listed by two sources (one record)
            "snp": {"param_value": {"1": "S83L", "2": "Q10Z"}, "Curated-R": {"3": "S83L", "4": "D87N"}},
            "40330": {"param_type": "multiple resistance variants", "param_value": {"5": "S83L,D87N"}},
            "41339": {"param_type": "nucleotide substitution in promoter region", "param_value": {"6": "-15c>t"}},
            "40438": {"param_type": "co-dependent single resistance variant", "param_value": {"7": "100,S83L+200,S80I"}},
        }, protein),
        "2": model("2", "3000002", "200", "Ecol_parC_FLO", "protein variant model",
                   {"snp": {"param_type": "single resistance variant", "param_value": {"8": "S80I"}}},
                   {"protein_sequence": {"accession": "P2", "sequence": "M" * 79 + "S" + "A" * 20}, **tax}),
        "3": model("3", "3000003", "300", "Ecol_23S_MLS", "rRNA gene variant model",
                   {"snp": {"param_type": "single resistance variant", "param_value": {"9": "a5g"}}},
                   {"dna_sequence": {"accession": "R1", "sequence": "CCCCACCCC"}, **tax}),
        "4": model("4", "3000004", "400", "x", "protein homolog model", {}, protein),
    }
    path.write_text(json.dumps(card))


def test_card(tmp_path):
    write_card(tmp_path / "card.json")
    out = ms.read_card(tmp_path / "card.json")
    recs = {r["native_notation"]: r for r in out["records"]}
    assert set(recs) == {"S83L", "D87N", "S83L,D87N", "-15c>t", "100,S83L+200,S80I", "S80I", "a5g"}
    assert recs["S83L"]["sources"] == ["Curated-R", "param_value"]
    assert out["rows_read"] == out["rows_parsed"] + len(out["excluded"])
    assert [e["reason"] for e in out["excluded"]] == ["CARD mutation with a residue that is not an amino acid: 'Q10Z'"]
    codep = recs["100,S83L+200,S80I"]["components"]
    assert [x["gene_label"] for x in codep] == ["gyrA", "parC"]
    assert recs["-15c>t"]["components"][0]["index"] is None          # no upstream sequence in CARD
    assert out["references"][recs["a5g"]["components"][0]["reference"]]["sequence"][4] == "A"


def write_ncbi(d):
    (d / "AMRProt.fa").write_text(f">WP_1.1|1|1|gyrA|gyrA|x\n{PROT}\n>WP_2.1|1|1|parC|parC|x\n{'M' + 'A' * 80 + 'S' + 'A' * 10}\n")
    (d / "AMRProt-mutation.tsv").write_text(
        "#taxgroup\taccession_version\tmutation_position\tstandard_mutation_symbol\treported_mutation_symbol\tclass\tsubclass\tmutated_protein_name\n"
        "Escherichia\tWP_1.1\t83\tgyrA_S83L\tgyrA_S83L\tQUINOLONE\tQUINOLONE\tx\n"
        "Escherichia\tWP_2.1\t82\tparC_S80I\tparC_S80I\tQUINOLONE\tQUINOLONE\tx\n"
        "Escherichia\tWP_9.1\t10\tfoo_A10V\tfoo_A10V\tQUINOLONE\tQUINOLONE\tx\n")
    region = "NC_1.1@ampC/blaEC_promoter_region@ampC:100-140"
    (d / "AMR_DNA-Escherichia.fa").write_text(f">{region}\n{'A' * 20 + 'C' + 'G' * 20}\n")
    (d / "AMR_DNA-Escherichia.tsv").write_text(
        "#accesion_version@gene_name@gene_symbol:start-stop\tmutation_position\tstandard_mutation_symbol\treported_mutation_symbol\tclass\tsubclass\tmutated_gene_name\n"
        f"{region}\t21\tampC_C-42T\tampC_C-42T\tBETA-LACTAM\tCEPHALOSPORIN\tx\n")


def test_ncbi(tmp_path):
    write_ncbi(tmp_path)
    out = ms.read_ncbi(tmp_path / "AMRProt-mutation.tsv", tmp_path / "AMRProt.fa", tmp_path)
    recs = {r["native_notation"]: r for r in out["records"]}
    assert set(recs) == {"gyrA_S83L", "parC_S80I", "ampC_C-42T"}
    assert out["rows_read"] == 4 and out["rows_parsed"] == 3 and len(out["excluded"]) == 1
    parc = recs["parC_S80I"]
    assert parc["components"][0]["start"] == 82 and "symbol numbering 80, sequence position 82" in parc["sources"]
    amp = recs["ampC_C-42T"]["components"][0]
    assert (amp["molecule"], amp["start"], amp["index"]) == ("promoter", -42, 21)
    assert out["references"][amp["reference"]]["sequence"][amp["index"] - 1] == "C"


def write_pointfinder(d):
    sp = d / "escherichia_coli"
    sp.mkdir()
    gene = _codons(PROT)
    (sp / "gyrA.fsa").write_text(f">gyrA_1_CP1.1\n{gene}\n")
    (sp / "parC.fsa").write_text(f">parC_1_CP2.1\n{_codons('M' * 79 + 'S' + 'A' * 20)}\n")
    promoter = "T" * 11 + "C" + "T" * 41          # size 53: position -42 is index 53 - 42 + 1 = 12
    (sp / "ampC-promoter-size-53bp.fsa").write_text(f">ampC-promoter_1_CP3.1\n{promoter}ATG\n")
    (sp / "RNA_genes.txt").write_text("23S\n")
    header = "#Gene_accession\tTypeGene\tMutation ID\tCodon_pos\tRef_nuc\tRef_codon\tRes_codon\tClass\tPhenotype\tPMID\tMechanism of resistance\tNotes\tRequired mutation\t% abundance mutation required\n"
    rows = [
        "gyrA_1_CP1.1\tAA\tSER-83\t83\tTCG\tS\tL,W\tQuinolone\tCiprofloxacin\t1\tTarget\t\t\t",
        "gyrA_1_CP1.1\tAA\tSER-83\t83\tTCG\tS\tL,W\tQuinolone\tNalidixic acid\t1\tTarget\t\t\t",   # same change, other phenotype
        "parC_1_CP2.1\tAA\tSER-80\t80\tAGC\tS\tI\tQuinolone\tCiprofloxacin\t1\tTarget\t\tgyrA_1_CP1.1_S83L.W\t",
        "ampC-promoter-size-53bp_1_CP3.1\tNUC\tCYS--42\t-42\tC\tC\tT\tBeta-lactam\tAmpicillin\t1\tx\t\t\t",
        "gyrA_2-CP9.1\tAA\tASP-87\t87\tGAT\tD\tN\tQuinolone\tCiprofloxacin\t1\tx\t\t\t",          # no sequence
    ]
    (sp / "phenotypes.txt").write_text(header + "\n".join(rows) + "\n")
    pl = d / "plasmodium_falciparum"
    pl.mkdir()
    (pl / "phenotypes.txt").write_text(header + "DHPS_1_X\tAA\tILE-431\t431\tATA\tI\tV\tx\tx\t1\tx\t\t\t\n")


def test_pointfinder(tmp_path):
    write_pointfinder(tmp_path)
    out = ms.read_pointfinder(tmp_path)
    recs = {r["native_notation"]: r for r in out["records"]}
    assert set(recs) == {"gyrA_S83L", "gyrA_S83W", "parC_S80I+gyrA_S83L", "parC_S80I+gyrA_S83W", "ampC_C-42T"}
    assert recs["gyrA_S83L"]["phenotype_native"] == "Ciprofloxacin; Nalidixic acid"
    assert out["rows_read"] == 6 and out["rows_parsed"] == 4
    assert sorted(e["reason"] for e in out["excluded"]) == ["no sequence for gyrA_2-CP9.1", "not a bacterium"]
    amp = recs["ampC_C-42T"]["components"][0]
    assert amp["index"] == 12 and out["references"][amp["reference"]]["sequence"][11] == "C"
    gyra = recs["gyrA_S83L"]["components"][0]
    assert out["references"][gyra["reference"]]["sequence"][82] == "S"      # translated reference
