"""Differences between two variants.json files: records added and removed, changed match classes, and loci whose
anchor changed (which moves the reported positions at that locus).
"""
from collections import Counter

NAMES = {"ncbi": "NCBI", "card": "CARD", "pointfinder": "PointFinder"}


def diff(old, new):
    o = {r["record_id"]: r for r in old["records"]}
    n = {r["record_id"]: r for r in new["records"]}
    om = {(m["record_id"], m["other_database"]): m["class"] for m in old["matches"]}
    nm = {(m["record_id"], m["other_database"]): m["class"] for m in new["matches"]}
    changed = sorted((k, om[k], nm[k]) for k in set(om) & set(nm) if om[k] != nm[k])
    anchors = []
    for lid in sorted(set(old["loci"]) & set(new["loci"])):
        a, b = old["loci"][lid], new["loci"][lid]
        if (a["anchor_database"], a["anchor_accession"], a["anchor_sha256"]) != (b["anchor_database"], b["anchor_accession"], b["anchor_sha256"]):
            anchors.append((lid, f"{a['anchor_database']} {a['anchor_accession']}", f"{b['anchor_database']} {b['anchor_accession']}"))
    return {
        "added": sorted(set(n) - set(o)), "removed": sorted(set(o) - set(n)), "class_changes": changed,
        "anchor_changes": anchors, "loci_added": sorted(set(new["loci"]) - set(old["loci"])),
        "loci_removed": sorted(set(old["loci"]) - set(new["loci"])),
        "records": (len(o), len(n)), "clusters": (len(old["clusters"]), len(new["clusters"])),
        "by_database": {db: (sum(1 for r in o.values() if r["database"] == db), sum(1 for r in n.values() if r["database"] == db))
                        for db in NAMES},
        "_new": n, "_old": o,
    }


def markdown(d, limit=50):
    L = ["## Mutation comparison changes", ""]
    L.append(f"Records {d['records'][0]} -> {d['records'][1]}; clusters {d['clusters'][0]} -> {d['clusters'][1]}.")
    L += ["", "| Database | Before | After |", "| --- | --- | --- |"]
    L += [f"| {NAMES[db]} | {a} | {b} |" for db, (a, b) in d["by_database"].items()]
    L.append("")

    def section(title, rows):
        L.append(f"### {title} ({len(rows)})")
        L.append("")
        L.extend(f"- {r}" for r in rows[:limit])
        if len(rows) > limit:
            L.append(f"- ... and {len(rows) - limit} more")
        L.append("")

    describe = lambda rid, recs: f"{NAMES[recs[rid]['database']]} {recs[rid]['native_notation']} ({recs[rid]['taxon_group']})"
    section("Anchor changes (reported positions move at these loci)", [f"{l}: {a} -> {b}" for l, a, b in d["anchor_changes"]])
    section("Records added", [describe(r, d["_new"]) for r in d["added"]])
    section("Records removed", [describe(r, d["_old"]) for r in d["removed"]])
    counts = Counter((a, b) for _, a, b in d["class_changes"])
    section("Match class changes", [f"{a} -> {b}: {n}" for (a, b), n in counts.most_common()])
    section("Loci added", d["loci_added"])
    section("Loci removed", d["loci_removed"])
    return "\n".join(L) + "\n"
