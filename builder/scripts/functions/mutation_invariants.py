"""Checks on a variants.json. check_variants returns {check: [examples]} for every check that fails."""
from collections import Counter

DATABASES = ("card", "ncbi", "pointfinder")
SINGLE_CLASSES = {
    "same_mutation", "same_mutation_different_reference", "possible_same_mutation_offset",
    "same_mutation_only_in_combination", "same_position_different_change", "different_mutation", "same_mutation_related_organism",
    "locus_not_covered", "organism_not_covered",
}
COMBINATION_CLASSES = {
    "same_combination", "partial_combination", "components_listed_separately", "combination_not_listed",
    "locus_not_covered", "organism_not_covered",
}
# Largest share of self-check mismatches per database (2026-11-19: NCBI 0%, PointFinder 0.6%, CARD 6.9%).
SELF_CHECK_MAX_MISMATCH = {"ncbi": 0.01, "card": 0.09, "pointfinder": 0.02}


def check_variants(data):
    v = {}

    def fail(name, example):
        v.setdefault(name, []).append(example)

    for db, s in data["sources"].items():
        whole = sum(1 for e in data["excluded"] if e["database"] == db and not e.get("partial"))
        if s["rows_read"] != s["rows_parsed"] + whole:
            fail("rows_accounted", f"{db}: read {s['rows_read']} != parsed {s['rows_parsed']} + excluded {whole}")
        if s["rows_parsed"] == 0:
            fail("rows_accounted", f"{db}: no rows parsed")

    ids = Counter(r["record_id"] for r in data["records"])
    for rid, n in ids.items():
        if n > 1:
            fail("unique_record_ids", rid)
    records = {r["record_id"]: r for r in data["records"]}

    for r in data["records"]:
        if not r.get("taxon_group"):
            fail("record_has_group", r["record_id"])
        for c in r["components"]:
            if c["locus_id"] not in data["loci"]:
                fail("component_has_locus", r["record_id"])
            if (c["anchor_start"] is None) != (c["anchor_end"] is None):
                fail("anchor_start_and_end", r["record_id"])

    seen = Counter((m["record_id"], m["other_database"]) for m in data["matches"])
    for r in data["records"]:
        for db in DATABASES:
            if db != r["database"] and seen[(r["record_id"], db)] != 1:
                fail("one_match_per_record_and_database", f"{r['record_id']} vs {db}: {seen[(r['record_id'], db)]}")
    for m in data["matches"]:
        r = records.get(m["record_id"])
        if r is None:
            fail("match_record_exists", m["record_id"])
            continue
        allowed = SINGLE_CLASSES if len(r["components"]) == 1 else COMBINATION_CLASSES
        if m["class"] not in allowed:
            fail("valid_class", f"{m['record_id']}: {m['class']}")
        for h in m["matched"]:
            if h not in records:
                fail("matched_record_exists", h)
            elif records[h]["database"] != m["other_database"]:
                fail("matched_record_in_other_database", h)

    # same-mutation classes are symmetric for single records
    index = {(m["record_id"], m["other_database"]): m for m in data["matches"]}
    same = {"same_mutation", "same_mutation_different_reference"}
    for m in data["matches"]:
        if m["class"] in same:
            for h in m["matched"]:
                back = index.get((h, m["database"]))
                if back is None or back["class"] not in same or m["record_id"] not in back["matched"]:
                    fail("same_mutation_symmetric", f"{m['record_id']} -> {h}")
        if m["class"] == "same_mutation":
            sha = records[m["record_id"]]["components"][0]["reference_sha256"]
            if not any(records[h]["components"][0]["reference_sha256"] == sha for h in m["matched"]):
                fail("same_mutation_same_reference", m["record_id"])

    in_clusters = Counter(rid for c in data["clusters"] for rid in c["members"])
    for rid in records:
        if in_clusters[rid] != 1:
            fail("record_in_one_cluster", f"{rid}: {in_clusters[rid]}")

    for db, bound in SELF_CHECK_MAX_MISMATCH.items():
        checks = Counter(c["ref_check"] for r in data["records"] if r["database"] == db for c in r["components"])
        checked = checks["ok"] + checks["mismatch"]
        if checked and checks["mismatch"] / checked > bound:
            fail("self_check_mismatch_bound", f"{db}: {checks['mismatch']}/{checked} above {bound}")
    return v
