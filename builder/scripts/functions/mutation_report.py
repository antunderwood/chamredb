"""Summary numbers and the markdown report of a mutation comparison."""
from collections import Counter, defaultdict

DATABASES = ("card", "ncbi", "pointfinder")
NAMES = {"ncbi": "NCBI", "card": "CARD", "pointfinder": "PointFinder"}
SAME = ("same_mutation", "same_mutation_different_reference")
CLASS_ORDER = (
    "same_mutation", "same_mutation_different_reference", "possible_same_mutation_offset",
    "same_mutation_only_in_combination", "same_position_different_change", "different_mutation",
    "same_mutation_related_organism",
    "locus_not_covered", "organism_not_covered",
)
COMBO_ORDER = ("same_combination", "partial_combination", "components_listed_separately", "combination_not_listed",
               "locus_not_covered", "organism_not_covered")


def summary(data):
    records = data["records"]
    by_id = {r["record_id"]: r for r in records}
    s = {"databases": {}, "groups": {}, "pairs": {}, "combination_pairs": {}, "clusters": {}, "phenotype": {}}
    for db in DATABASES:
        mine = [r for r in records if r["database"] == db]
        comps = [c for r in mine for c in r["components"]]
        s["databases"][db] = {
            "records": len(mine),
            "single": sum(1 for r in mine if len(r["components"]) == 1),
            "combination": sum(1 for r in mine if len(r["components"]) > 1),
            "groups": len({r["taxon_group"] for r in mine}),
            "molecule": dict(Counter(c["molecule"] for c in comps)),
            "self_check": dict(Counter(c["ref_check"] for c in comps)),
            "anchor_check": dict(Counter(c["anchor_ref_check"] for c in comps)),
            "not_placed": sum(1 for c in comps if c["anchor_start"] is None),
            "excluded": sum(1 for e in data["excluded"] if e["database"] == db and not e.get("partial")),
        }
    groups = defaultdict(Counter)
    for r in records:
        groups[r["taxon_group"]][r["database"]] += 1
    s["groups"] = {g: dict(c) for g, c in sorted(groups.items())}
    s["groups_in_two_or_more"] = sorted(g for g, c in groups.items() if len(c) >= 2)
    s["groups_in_all_three"] = sorted(g for g, c in groups.items() if len(c) == 3)

    for m in data["matches"]:
        r = by_id[m["record_id"]]
        key = f"{m['database']}->{m['other_database']}"
        bucket = s["pairs"] if len(r["components"]) == 1 else s["combination_pairs"]
        bucket.setdefault(key, Counter())[m["class"]] += 1
        if m["class"] in SAME and m["phenotype"]:
            s["phenotype"].setdefault(key, Counter())[m["phenotype"]] += 1
    for bucket in ("pairs", "combination_pairs", "phenotype"):
        s[bucket] = {k: dict(v) for k, v in sorted(s[bucket].items())}

    loci_dbs = Counter(len({ref["database"] for ref in l["references"]}) for l in data["loci"].values())
    s["loci"] = {"total": len(data["loci"]), "by_databases": dict(loci_dbs),
                 "joined_by_sequence": len(data["loci_joined_by_sequence"])}

    clusters = data["clusters"]
    s["clusters"]["total"] = len(clusters)
    s["clusters"]["by_databases"] = dict(Counter(c["databases"] for c in clusters))
    pattern = defaultdict(Counter)
    for c in clusters:
        if len(c["members"]) and c["kind"] != "combination":
            listed = "+".join(db for db in DATABASES if c["status"][db] == "listed")
            pattern[c["group"]][listed] += 1
    s["clusters"]["patterns_by_group"] = {g: dict(v) for g, v in sorted(pattern.items()) if g in s["groups_in_two_or_more"]}
    agree = {}
    for a in DATABASES:
        for b in DATABASES:
            if a >= b:
                continue
            both = [c for c in clusters if c["kind"] != "combination" and c["group"] in
                    {g for g, cnt in groups.items() if cnt.get(a) and cnt.get(b)}
                    and (c["status"][a] == "listed" or c["status"][b] == "listed")]
            n_both = sum(1 for c in both if c["status"][a] == "listed" and c["status"][b] == "listed")
            agree[f"{a}+{b}"] = {"clusters": len(both), "listed_by_both": n_both}
    s["pair_overlap_shared_groups"] = agree
    s["systematic_offsets"] = len(data["systematic_offsets"])
    s["self_check_offsets"] = len(data["self_check_offsets"])
    return s


def _pct(n, d):
    return f"{100 * n / d:.1f}%" if d else "n/a"


def _table(header, rows):
    out = ["| " + " | ".join(header) + " |", "| " + " | ".join("---" for _ in header) + " |"]
    out += ["| " + " | ".join(str(x) for x in row) + " |" for row in rows]
    return "\n".join(out)


def markdown(data, violations):
    s = data["summary"]
    L = [f"# Mutation comparison report{(' (data release ' + data['data_release'] + ')') if data['data_release'] else ''}", ""]
    L.append("Built by `builder/scripts/compare_mutations.py`. Positions are in the numbering of each locus's "
             "anchor database (CARD, else NCBI, else PointFinder); rRNA positions are on the sense strand of the anchor.")
    L.append("")
    if data["inputs"]:
        L += ["## Inputs", "", _table(["Database", "Role", "Release"],
                                       [[NAMES.get(db, db), role, v["release"]] for db, roles in sorted(data["inputs"].items())
                                        for role, v in sorted(roles.items()) if db.endswith("mutations") or db in ("pointfinder",)
                                        or role.startswith("mutation")]), ""]
    L += ["## Invariants", "", "All hold." if not violations else
          "\n".join(f"- **{k}**: {len(v)} (e.g. {v[:3]})" for k, v in violations.items()), ""]

    L += ["## Rows and records", ""]
    rows = []
    for db in DATABASES:
        d, src = s["databases"][db], data["sources"][db]
        rows.append([NAMES[db], src["rows_read"], src["rows_parsed"], d["excluded"], d["records"], d["single"], d["combination"], d["groups"]])
    L += [_table(["Database", "Upstream rows", "Parsed", "Excluded", "Records", "Single", "Combination", "Organism groups"], rows), ""]
    L.append("Records are distinct changes: rows listing the same change twice (CARD once per evidence source, "
             "PointFinder once per phenotype) are merged, and rows listing several alternative residues are split.")
    L.append("")
    reasons = Counter((e["database"], e["reason"] if len(e["reason"]) < 90 else e["reason"][:87] + "...")
                      for e in data["excluded"] if not e.get("partial"))
    if reasons:
        L += ["### Excluded rows", "", _table(["Database", "Reason", "Rows"],
                                             [[NAMES[db], r.replace("|", "/"), n] for (db, r), n in sorted(reasons.items())]), ""]

    L += ["## Reference checks", ""]
    L.append("Self check: the residue each database states as the original one, against its own reference sequence. "
             "Anchor check: the same residue against the anchor sequence at the converted position.")
    L.append("")
    rows = []
    for db in DATABASES:
        d = s["databases"][db]
        sc, ac = d["self_check"], d["anchor_check"]
        rows.append([NAMES[db], sc.get("ok", 0), sc.get("mismatch", 0), _pct(sc.get("mismatch", 0), sc.get("ok", 0) + sc.get("mismatch", 0)),
                     sc.get("unchecked", 0), ac.get("ok", 0), ac.get("mismatch", 0), d["not_placed"]])
    L += [_table(["Database", "Self ok", "Self mismatch", "Mismatch share", "Unchecked", "Anchor ok", "Anchor mismatch", "Not placed"], rows), ""]
    if data["self_check_offsets"]:
        L += ["### Self-check mismatches explained by one constant shift", "",
              "A shift that makes at least three mismatching residues agree suggests the database numbers this gene in "
              "another convention (often *E. coli* numbering). Not corrected: listed for review.", "",
              _table(["Organism", "Gene", "Changes", "Mismatches", "Shift", "Explained"],
                     [[o["organism"], o["gene"], o["changes"], o["mismatches"], f"{o['shift']:+d}", o["mismatches_explained"]]
                      for o in data["self_check_offsets"]]), ""]

    odd = [o for o in data.get("rrna_strands", []) if (o["stored_strand"], o["numbering_strand"]) != ("sense", "sense")]
    if odd:
        L += ["### rRNA references not stored or numbered along the sense strand", "",
              "The stored strand is found from conserved rRNA stretches. Anchor positions of rRNA are always on the "
              "sense strand (the rRNA itself). Positions are read along the sense strand unless the stated original "
              f"residues fit the antisense strand at {data['settings']['rrna_antisense_margin']} more sites; a "
              "reference with unknown strand is read as stored.", "",
              _table(["Locus", "Database", "Reference", "Stored", "Positions read along", "Sites", "Fit sense", "Fit antisense"],
                     [[o["locus_id"], NAMES[o["database"]], o["reference"].split("|")[1], o["stored_strand"], o["numbering_strand"],
                       o["sites"], "" if o["sites_fit_sense"] is None else o["sites_fit_sense"],
                       "" if o["sites_fit_antisense"] is None else o["sites_fit_antisense"]] for o in odd]), ""]

    L += ["## Loci", ""]
    lb = s["loci"]["by_databases"]
    L.append(f"{s['loci']['total']} loci: {lb.get(3, 0)} with references from all three databases, {lb.get(2, 0)} from two, "
             f"{lb.get(1, 0)} from one. {s['loci']['joined_by_sequence']} pairs of differently named references were joined "
             f"because their sequences are near-identical.")
    L.append("")
    if data["loci_joined_by_sequence"]:
        L += [_table(["Group", "Molecule", "Reference", "Joined with", "Identity", "Coverage"],
                     [[j["group"], j["molecule"], j["a"], j["b"], j["identity"], j["coverage"]] for j in data["loci_joined_by_sequence"]]), ""]

    L += ["## Organism coverage", "", _table(["Organism group"] + [NAMES[d] for d in DATABASES],
                                             [[g] + [c.get(d, 0) for d in DATABASES] for g, c in s["groups"].items() if len(c) >= 2]),
          "", f"Only organism groups with mutations in at least two databases are shown ({len(s['groups_in_two_or_more'])} of "
          f"{len(s['groups'])}).", ""]

    L += ["## Match classes (single mutations)", ""]
    pairs = list(s["pairs"].keys())
    L += [_table(["Class"] + pairs, [[c] + [s["pairs"][p].get(c, 0) for p in pairs] for c in CLASS_ORDER]), ""]
    if s["combination_pairs"]:
        cp = list(s["combination_pairs"].keys())
        L += ["## Match classes (combinations)", "",
              _table(["Class"] + cp, [[c] + [s["combination_pairs"][p].get(c, 0) for p in cp] for c in COMBO_ORDER]), ""]

    L += ["## Agreement where both databases cover the organism", ""]
    L += [_table(["Pair", "Mutations (clusters) in shared organism groups", "Listed by both", "Share"],
                 [[k.replace("+", " and ").replace("ncbi", "NCBI").replace("card", "CARD").replace("pointfinder", "PointFinder"),
                   v["clusters"], v["listed_by_both"], _pct(v["listed_by_both"], v["clusters"])]
                  for k, v in s["pair_overlap_shared_groups"].items()]), ""]

    L += ["## Mutations per organism group and databases listing them", ""]
    pats = sorted({p for v in s["clusters"]["patterns_by_group"].values() for p in v})
    L += [_table(["Group"] + pats, [[g] + [v.get(p, 0) for p in pats] for g, v in s["clusters"]["patterns_by_group"].items()]), ""]

    if s["phenotype"]:
        L += ["## Drug classes of matched mutations", "", _table(
            ["Pair", "same_classes", "shared_class", "no_shared_class", "unknown"],
            [[k] + [v.get(c, 0) for c in ("same_classes", "shared_class", "no_shared_class", "unknown")] for k, v in s["phenotype"].items()]), ""]

    if data["systematic_offsets"]:
        L += ["## Systematic offsets between databases", "",
              "The same change at a constant distance on the anchor, seen at least three times: probably a numbering "
              "convention, not different mutations. Review in `offset_review.tsv`.", "",
              _table(["Locus", "Database", "Other", "Offset", "Count"],
                     [[o["locus_id"], NAMES[o["database"]], NAMES[o["other_database"]], f"{o['offset']:+d}", o["count"]]
                      for o in data["systematic_offsets"]]), ""]
    return "\n".join(L) + "\n"
