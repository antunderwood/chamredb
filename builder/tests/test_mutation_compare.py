"""Match classes on a small designed set of records, one case per class. The sequences are part of E. coli gyrA;
PointFinder's reference has two extra N-terminal residues, and one CARD record gives a residue one position off.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from functions import mutation_compare as mc  # noqa: E402
from functions import mutation_invariants  # noqa: E402
from functions import mutation_notation as mn  # noqa: E402

GYRA = "MSDLAREITPVNIEEELKSSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMNVLGNDWNKAYKKSARVVGDVIGKYHPHGDSAVYDTIVRMAQPFSLRYMLVDGQGNFGSIDGDSAAAMRYTEIRLAKIAHELMADLEKETVDFVDNYDGTEKIPDVMPTKIPNLLVNG"
SAL = GYRA[:9] + "K" + GYRA[10:]           # Salmonella: one substitution away
PARC = "MSDMAERLALHEFTENAYLNYSMYVIMDRALPFIGDGLKPVQRRIVYAMSELGLNASAKFKKSARTVGDVLGKYHPHGDSACYEAMVLMAQPFSYRYPLVDGQGNWG"

TAXA = [dict(database=db, label=label, group=group, group_source="", group_family=fam, note="")
        for db, label, group, fam in [
            ("ncbi", "Escherichia", "Escherichia", "Enterobacteriaceae"),
            ("ncbi", "Salmonella", "Salmonella", "Enterobacteriaceae"),
            ("card", "Escherichia coli", "Escherichia", "Enterobacteriaceae"),
            ("card", "Mycobacterium tuberculosis", "Mycobacterium tuberculosis", "Mycobacteriaceae"),
            ("pointfinder", "escherichia_coli", "Escherichia", "Enterobacteriaceae")]]

REFS = {
    "ncbi|N1|protein": ("ncbi", "gyrA", "Escherichia", GYRA),
    "ncbi|N2|protein": ("ncbi", "parC", "Escherichia", PARC),
    "ncbi|S1|protein": ("ncbi", "gyrA", "Salmonella", SAL),
    "card|C1|protein": ("card", "gyrA", "Escherichia coli K-12", GYRA),
    "card|M1|protein": ("card", "rpoB", "Mycobacterium tuberculosis H37Rv", "MSAAA"),
    "pointfinder|F1|protein": ("pointfinder", "gyrA", "escherichia_coli", "MA" + GYRA),
}


def sub(ref, pos, alt, residue_from=None):
    """Substitution at pos of ref; residue_from states the original residue read at another position (an error)."""
    seq = REFS[ref][3]
    r = seq[(residue_from or pos) - 1]
    c = mn.substitution("protein", pos, r, alt)
    c.update(gene_label=REFS[ref][1], reference=ref, index=pos)
    return c


def record(ref, *components):
    db, _, organism, _ = REFS[ref]
    comps = list(components)
    return {"database": db, "record_id": ms_id(db, comps), "native_notation": "+".join(mn.compact(c) for c in comps),
            "native_label": ref, "organism_label": organism, "components": comps, "phenotype_native": "",
            "drug_classes_native": [], "sources": []}


def ms_id(db, comps):
    return f"{db}|" + "+".join(sorted(f"{c['reference']}:{mn.hgvs(c)}" for c in comps))


@pytest.fixture(scope="module")
def result():
    records = [
        record("ncbi|N1|protein", sub("ncbi|N1|protein", 83, "L")),     # A: in all three
        record("card|C1|protein", sub("card|C1|protein", 83, "L")),     # identical reference -> same_mutation
        record("pointfinder|F1|protein", sub("pointfinder|F1|protein", 85, "L")),   # +2 numbering -> different reference
        record("ncbi|N1|protein", sub("ncbi|N1|protein", 87, "N")),     # B
        record("card|C1|protein", sub("card|C1|protein", 88, "N", residue_from=87)),  # residue of 87 written at 88
        record("ncbi|N1|protein", sub("ncbi|N1|protein", 100, "W")),    # C: PointFinder has 100 -> F
        record("pointfinder|F1|protein", sub("pointfinder|F1|protein", 102, "F")),
        record("ncbi|N1|protein", sub("ncbi|N1|protein", 120, "W")),    # D: only NCBI
        record("ncbi|N1|protein", sub("ncbi|N1|protein", 40, "W")),     # E: PointFinder only in a combination
        record("pointfinder|F1|protein", sub("pointfinder|F1|protein", 42, "W"), sub("pointfinder|F1|protein", 52, "F")),
        record("card|C1|protein", sub("card|C1|protein", 40, "W"), sub("card|C1|protein", 50, "F")),
        record("ncbi|N2|protein", sub("ncbi|N2|protein", 80, "I")),     # F: parC only in NCBI
        record("card|M1|protein", sub("card|M1|protein", 3, "V")),       # G: M. tuberculosis only in CARD
        record("ncbi|S1|protein", sub("ncbi|S1|protein", 83, "L")),     # H: Salmonella, related to E. coli
    ]
    references = {k: {"database": v[0], "accession": k.split("|")[1], "molecule": "protein", "gene_label": v[1],
                      "organism_label": v[2], "sequence": v[3]} for k, v in REFS.items()}
    mc.assign_groups(records, TAXA)
    mc.normalise_classes(records, [])
    loci, joined = mc.match_loci(records, references, {})
    mc.place(records, loci, references, mc.lift_maps(loci, references))
    matches = mc.classify(records, loci, references)
    clusters = mc.cluster(records, matches)
    index = {(m["record_id"], m["other_database"]): m for m in matches}
    by_note = {r["native_notation"] + "@" + r["native_label"]: r for r in records}
    data = {"sources": {}, "excluded": [], "records": records, "loci": loci, "matches": matches, "clusters": clusters}
    return index, by_note, data


def cls(result, note, other):
    index, by_note, _ = result
    return index[(by_note[note]["record_id"], other)]["class"]


@pytest.mark.parametrize("note,other,expected", [
    ("S83L@ncbi|N1|protein", "card", "same_mutation"),
    ("S83L@ncbi|N1|protein", "pointfinder", "same_mutation_different_reference"),
    ("S85L@pointfinder|F1|protein", "ncbi", "same_mutation_different_reference"),
    ("D87N@ncbi|N1|protein", "card", "possible_same_mutation_offset"),
    ("W100W@ncbi|N1|protein".replace("W100W", f"{GYRA[99]}100W"), "pointfinder", "same_position_different_change"),
    (f"{GYRA[119]}120W@ncbi|N1|protein", "card", "different_mutation"),
    (f"{GYRA[39]}40W@ncbi|N1|protein", "pointfinder", "same_mutation_only_in_combination"),
    (f"{PARC[79]}80I@ncbi|N2|protein", "card", "locus_not_covered"),
    ("S3V@card|M1|protein".replace("S3V", "A3V"), "ncbi", "organism_not_covered"),
    ("S83L@ncbi|S1|protein", "card", "same_mutation_related_organism"),
])
def test_single_classes(result, note, other, expected):
    assert cls(result, note, other) == expected


def test_combinations(result):
    pf = f"{GYRA[39]}42W+{GYRA[49]}52F@pointfinder|F1|protein"
    card = f"{GYRA[39]}40W+{GYRA[49]}50F@card|C1|protein"
    assert cls(result, card, "pointfinder") == "same_combination"
    assert cls(result, card, "ncbi") == "combination_not_listed"      # NCBI has the first change only


def test_offset_needs_a_reference_check_failure(result):
    index, by_note, _ = result
    card = by_note[f"{GYRA[86]}88N@card|C1|protein"]
    assert card["components"][0]["ref_check"] == "mismatch"
    assert index[(card["record_id"], "ncbi")]["offset"] == -1


def test_positions_are_reported_on_the_card_anchor(result):
    _, by_note, _ = result
    pf = by_note["S85L@pointfinder|F1|protein"]["components"][0]
    assert (pf["anchor_start"], pf["hgvs"]) == (83, "p.Ser83Leu")


def test_clusters_join_same_mutations(result):
    _, by_note, data = result
    c = next(c for c in data["clusters"] if by_note["S83L@ncbi|N1|protein"]["record_id"] in c["members"])
    assert c["status"] == {"ncbi": "listed", "card": "listed", "pointfinder": "listed"} and c["hgvs"] == "p.Ser83Leu"


def test_invariants_hold_and_catch_a_broken_table(result):
    _, _, data = result
    # the designed CARD numbering error is 1 of 5 CARD changes, above the bound meant for real releases
    assert set(mutation_invariants.check_variants(data)) == {"self_check_mismatch_bound"}
    broken = dict(data, matches=data["matches"][1:])
    assert "one_match_per_record_and_database" in mutation_invariants.check_variants(broken)


def test_unknown_organism_label_is_an_error():
    r = record("ncbi|N1|protein", sub("ncbi|N1|protein", 83, "L"))
    r["organism_label"] = "Unknown_bug"
    with pytest.raises(ValueError, match="taxa.tsv"):
        mc.assign_groups([r], TAXA)


def test_class_orders_do_not_overwrite_each_other():
    # both lists hold locus_not_covered and organism_not_covered, at different ranks
    assert mc.PRIORITY["same_mutation_related_organism"] < mc.PRIORITY["locus_not_covered"] < mc.PRIORITY["organism_not_covered"]
    assert mc.PRIORITY["different_mutation"] < mc.PRIORITY["locus_not_covered"]
    assert mc.COMBINATION_PRIORITY["combination_not_listed"] < mc.COMBINATION_PRIORITY["locus_not_covered"]


def _rrna_result(pointfinder_sequence, pointfinder_change):
    """CARD and PointFinder list one 16S change; PointFinder's reference is given by the caller."""
    import random
    rng = random.Random(7)
    rrs = "".join(rng.choice("ACGT") for _ in range(400))
    refs = {"card|R1|rna": ("card", "16S", "Escherichia coli K-12", rrs),
            "pointfinder|R2|rna": ("pointfinder", "16S", "escherichia_coli", pointfinder_sequence(rrs))}
    references = {k: {"database": v[0], "accession": k.split("|")[1], "molecule": "rna", "gene_label": v[1],
                      "organism_label": v[2], "sequence": v[3]} for k, v in refs.items()}

    def rec(ref, c):
        c.update(gene_label="16S", reference=ref, index=c["start"])
        return {"database": refs[ref][0], "record_id": ms_id(refs[ref][0], [c]), "native_notation": mn.compact(c),
                "native_label": ref, "organism_label": refs[ref][2], "components": [c], "phenotype_native": "",
                "drug_classes_native": [], "sources": []}
    alt = "A" if rrs[99] != "A" else "C"
    records = [rec("card|R1|rna", mn.substitution("rna", 100, rrs[99], alt)), rec("pointfinder|R2|rna", pointfinder_change(rrs, alt))]
    mc.assign_groups(records, TAXA)
    mc.normalise_classes(records, [])
    loci, _ = mc.match_loci(records, references, {})
    mc.place(records, loci, references, mc.lift_maps(loci, references))
    return records, loci, mc.classify(records, loci, references)


def test_reference_on_the_other_strand_is_placed_and_complemented():
    from functions import liftover
    # PointFinder stores the reverse complement: anchor position 100 is its position 301, residues complemented
    records, loci, matches = _rrna_result(
        liftover.revcomp, lambda rrs, alt: mn.substitution("rna", 301, liftover.revcomp(rrs[99]), liftover.revcomp(alt)))
    assert len(loci) == 1
    pf = records[1]["components"][0]
    assert (pf["anchor_start"], pf["anchor_strand"], pf["ref_check"], pf["anchor_ref_check"]) == (100, -1, "ok", "ok")
    assert pf["hgvs"] == records[0]["components"][0]["hgvs"]
    assert {m["class"] for m in matches if m["other_database"] in ("card", "pointfinder")} == {"same_mutation_different_reference"}


def test_same_notation_on_opposite_strands_is_left_for_review():
    from functions import liftover
    # PointFinder writes the same number and residues as CARD but its reference is the other strand
    records, loci, matches = _rrna_result(liftover.revcomp, lambda rrs, alt: mn.substitution("rna", 100, rrs[99], alt))
    assert {m["class"] for m in matches if m["other_database"] in ("card", "pointfinder")} == {"possible_same_mutation_offset"}


def test_unrelated_reference_is_not_lifted_on_a_few_aligned_columns():
    import random
    other = "".join(random.Random(99).choice("ACGT") for _ in range(400))
    records, loci, matches = _rrna_result(lambda rrs: other, lambda rrs, alt: mn.substitution("rna", 100, other[99], "A" if other[99] != "A" else "C"))
    assert records[1]["components"][0]["anchor_start"] is None
    entry = [e for l in loci.values() for e in l["references"] if e["reference"] == "pointfinder|R2|rna"][0]
    assert entry["lifted"] is False


def test_any_residue_record_does_not_join_different_changes_into_one_cluster():
    def rec(db, rid, kind):
        return {"record_id": rid, "database": db, "taxon_group": "g", "locus_ids": ["g:katG:protein"], "kind": kind,
                "components": [{"kind": kind, "hgvs": rid, "native_hgvs": rid}]}
    records = [rec("card", "c|Var", "any_substitution"), rec("card", "c|T", "substitution"),
               rec("pointfinder", "p|T", "substitution"), rec("pointfinder", "p|N", "substitution")]

    def match(a, db, other, hits):
        return {"record_id": a, "database": db, "other_database": other, "class": "same_mutation_different_reference", "matched": hits}
    matches = [match("c|Var", "card", "pointfinder", ["p|T", "p|N"]), match("c|T", "card", "pointfinder", ["p|T"]),
               match("p|T", "pointfinder", "card", ["c|T", "c|Var"]), match("p|N", "pointfinder", "card", ["c|Var"])]
    members = sorted(sorted(c["members"]) for c in mc.cluster(records, matches))
    assert members == [["c|T", "p|T"], ["c|Var"], ["p|N"]]
    var = [c for c in mc.cluster(records, matches) if c["members"] == ["c|Var"]][0]
    assert var["status"]["pointfinder"] == "same_mutation_different_reference"   # still reported as covered


def _oriented_rrna(card_changes, pointfinder_changes, card_stored_antisense=True):
    """A 16S gene with conserved stretches (so its strand can be told); CARD stores it antisense if asked, PointFinder
    sense. The callers give (position, original residue, new residue) as each database writes them."""
    import random
    from functions import liftover
    rng = random.Random(11)
    part = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    sense = part(200) + "ACGGGAGGCAGCAG" + part(100) + "CACACCGCCCGTC" + part(73)
    seqs = {"card|R1|rna": liftover.revcomp(sense) if card_stored_antisense else sense, "pointfinder|R2|rna": sense}
    org = {"card": "Escherichia coli K-12", "pointfinder": "escherichia_coli"}
    references = {k: {"database": k.split("|")[0], "accession": k.split("|")[1], "molecule": "rna", "gene_label": "16S",
                      "organism_label": org[k.split("|")[0]], "sequence": v} for k, v in seqs.items()}
    records = []
    for ref, changes in (("card|R1|rna", card_changes(sense)), ("pointfinder|R2|rna", pointfinder_changes(sense))):
        for pos, old, new in changes:
            c = mn.substitution("rna", pos, old, new)
            c.update(gene_label="16S", reference=ref, index=pos)
            db = ref.split("|")[0]
            records.append({"database": db, "record_id": ms_id(db, [c]), "native_notation": mn.compact(c), "native_label": ref,
                            "organism_label": org[db], "components": [c], "phenotype_native": "", "drug_classes_native": [],
                            "sources": []})
    mc.assign_groups(records, TAXA)
    mc.normalise_classes(records, [])
    loci, _ = mc.match_loci(records, references, {})
    strands = mc.orient_rrna(records, loci, references)
    mc.place(records, loci, references, mc.lift_maps(loci, references))
    return sense, records, strands, mc.classify(records, loci, references)


def _other(base):
    return "A" if base != "A" else "C"


def test_rrna_strand_is_told_from_conserved_stretches():
    import random
    from functions import liftover
    rng = random.Random(3)
    sense = "".join(rng.choice("ACGT") for _ in range(300)) + "GCCAGCAGCCGCGGTAA" + "".join(rng.choice("ACGT") for _ in range(300))
    assert liftover.rrna_strand(sense) == 1
    assert liftover.rrna_strand(liftover.revcomp(sense).lower().replace("t", "u")) == -1
    assert liftover.rrna_strand(sense[:300]) == 0


def test_antisense_stored_rrna_is_read_and_anchored_on_the_sense_strand():
    # CARD stores the reverse complement but numbers the rRNA itself, as PointFinder does on its sense sequence
    site = lambda sense: [(100, sense[99], _other(sense[99]))]
    sense, records, strands, matches = _oriented_rrna(site, site)
    card = [s for s in strands if s["database"] == "card"][0]
    assert (card["stored_strand"], card["numbering_strand"]) == ("antisense", "sense")
    c = records[0]["components"][0]
    assert (c["anchor_start"], c["anchor_strand"], c["ref_check"], c["anchor_ref_check"]) == (100, 1, "ok", "ok")
    assert c["hgvs"] == f"n.100{sense[99]}>{_other(sense[99])}" == records[1]["components"][0]["hgvs"]
    assert {m["class"] for m in matches if m["other_database"] in ("card", "pointfinder")} == {"same_mutation_different_reference"}


def test_rrna_numbered_along_the_stored_antisense_strand_is_turned_onto_the_sense_strand():
    from functions import liftover

    def card_changes(sense):   # three sites written on the stored (antisense) strand that do not fit the sense strand
        anti = liftover.revcomp(sense)
        sites = [p for p in range(20, 180) if anti[p - 1] != sense[p - 1]][:3]
        return [(p, anti[p - 1], _other(anti[p - 1])) for p in sites]

    def pointfinder_changes(sense):   # the same changes, written on the sense strand
        return [(len(sense) - p + 1, liftover.revcomp(old), liftover.revcomp(new)) for p, old, new in card_changes(sense)]

    sense, records, strands, matches = _oriented_rrna(card_changes, pointfinder_changes)
    card = [s for s in strands if s["database"] == "card"][0]
    assert (card["numbering_strand"], card["sites_fit_sense"], card["sites_fit_antisense"]) == ("antisense", 0, 3)
    for r in records[:3]:
        c = r["components"][0]
        assert (c["anchor_start"], c["anchor_strand"], c["ref_check"]) == (len(sense) - c["start"] + 1, -1, "ok")
    assert {m["class"] for m in matches if m["other_database"] in ("card", "pointfinder")} == {"same_mutation_different_reference"}


def test_two_sites_fitting_the_antisense_strand_do_not_turn_the_numbering():
    from functions import liftover

    def card_changes(sense):
        anti = liftover.revcomp(sense)
        sites = [p for p in range(20, 180) if anti[p - 1] != sense[p - 1]][:2]
        return [(p, anti[p - 1], _other(anti[p - 1])) for p in sites]

    _, records, strands, _ = _oriented_rrna(card_changes, lambda sense: [(100, sense[99], _other(sense[99]))])
    assert [s for s in strands if s["database"] == "card"][0]["numbering_strand"] == "sense"
    assert {r["components"][0]["ref_check"] for r in records[:2]} == {"mismatch"}
