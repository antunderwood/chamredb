"""README numbers and figures for the mutation comparison, from a variants.json.

    python -m analysis.variants --variants chamredb/data/variants.json [--out docs]

Writes <out>/variants/statistics.json and <out>/images/figure_variant_{coverage,classes}.{svg,png} with their source data.
Deterministic: sorted keys and rows, no timestamps.
"""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

DATABASES = ("card", "ncbi", "pointfinder")
LABELS = {"ncbi": "NCBI", "card": "CARD", "pointfinder": "PointFinder"}
SAME = ("same_mutation", "same_mutation_different_reference")
CLASS_ORDER = ("same_mutation", "same_mutation_different_reference", "possible_same_mutation_offset",
               "same_mutation_only_in_combination", "same_position_different_change", "different_mutation",
               "same_mutation_related_organism", "locus_not_covered")
CLASS_LABELS = {
    "same_mutation": "Same mutation, same reference", "same_mutation_different_reference": "Same mutation, other reference",
    "possible_same_mutation_offset": "Possibly same, numbered differently", "same_mutation_only_in_combination": "Listed only in a combination",
    "same_position_different_change": "Same position, other change", "different_mutation": "Gene covered, mutation not listed",
    "same_mutation_related_organism": "Listed for a related organism", "locus_not_covered": "Gene not covered for the organism",
}


def statistics(data):
    records = data["records"]
    by_id = {r["record_id"]: r for r in records}
    groups = defaultdict(set)
    for r in records:
        groups[r["taxon_group"]].add(r["database"])
    single = [r for r in records if len(r["components"]) == 1]
    out = {
        "data_release": data["data_release"],
        "records": {db: sum(1 for r in records if r["database"] == db) for db in DATABASES},
        "records_total": len(records),
        "combinations": sum(1 for r in records if len(r["components"]) > 1),
        "rows_read": {db: data["sources"][db]["rows_read"] for db in DATABASES},
        "rows_excluded": {db: sum(1 for e in data["excluded"] if e["database"] == db and not e.get("partial")) for db in DATABASES},
        "organism_groups": {"total": len(groups), "two_or_more": sum(1 for v in groups.values() if len(v) >= 2),
                            "all_three": sum(1 for v in groups.values() if len(v) == 3)},
        "organism_groups_by_database": {db: sum(1 for v in groups.values() if db in v) for db in DATABASES},
        "loci": {"total": len(data["loci"]),
                 "by_databases": dict(sorted(Counter(len({x["database"] for x in l["references"]}) for l in data["loci"].values()).items()))},
        "card_tuberculosis_share": sum(1 for r in records if r["database"] == "card" and r["taxon_group"] == "mycobacterium_tuberculosis")
        / max(1, sum(1 for r in records if r["database"] == "card")),
        "self_check": {}, "pairs": {}, "shared_organisms": {},
    }
    for db in DATABASES:
        c = Counter(x["ref_check"] for r in records if r["database"] == db for x in r["components"])
        out["self_check"][db] = {"ok": c["ok"], "mismatch": c["mismatch"], "unchecked": c["unchecked"],
                                 "mismatch_share": c["mismatch"] / max(1, c["ok"] + c["mismatch"])}
    for m in data["matches"]:
        r = by_id[m["record_id"]]
        if len(r["components"]) == 1:
            out["pairs"].setdefault(f"{m['database']}>{m['other_database']}", Counter())[m["class"]] += 1
    out["pairs"] = {k: dict(sorted(v.items())) for k, v in sorted(out["pairs"].items())}
    clusters = [c for c in data["clusters"] if c["kind"] != "combination"]
    for i, a in enumerate(DATABASES):
        for b in DATABASES[i + 1:]:
            shared = {g for g, dbs in groups.items() if a in dbs and b in dbs}
            rel = [c for c in clusters if c["group"] in shared and (c["status"][a] == "listed" or c["status"][b] == "listed")]
            both = sum(1 for c in rel if c["status"][a] == "listed" and c["status"][b] == "listed")
            out["shared_organisms"][f"{a}+{b}"] = {"organism_groups": len(shared), "mutations": len(rel), "listed_by_both": both,
                                                  "share": both / max(1, len(rel))}
    out["clusters"] = {"single": len(clusters),
                       "by_databases": dict(sorted(Counter(c["databases"] for c in clusters).items()))}
    offset_pairs = {tuple(sorted((m["record_id"], h))) for m in data["matches"]
                    if m["class"] == "possible_same_mutation_offset" for h in m["matched"]}
    out["offset_pairs"] = len(offset_pairs)
    out["systematic_offsets"] = len({(o["locus_id"], tuple(sorted((o["database"], o["other_database"])))) for o in data["systematic_offsets"]})
    out["single_mutations"] = len(single)
    return out


def display(group):
    """Organism group as text: Klebsiella_pneumoniae -> Klebsiella pneumoniae."""
    text = group.replace("_", " ")
    return text[:1].upper() + text[1:]


def coverage_rows(data):
    """Single mutations (clusters) per organism group covered by at least two databases, by the databases listing them."""
    groups = defaultdict(set)
    for r in data["records"]:
        groups[r["taxon_group"]].add(r["database"])
    rows = defaultdict(Counter)
    for c in data["clusters"]:
        if c["kind"] == "combination" or len(groups[c["group"]]) < 2:
            continue
        listed = [db for db in DATABASES if c["status"][db] == "listed"]
        key = f"{LABELS[listed[0]]} only" if len(listed) == 1 else ("All three" if len(listed) == 3 else "Two databases")
        rows[c["group"]][key] += 1
    return {g: dict(v) for g, v in rows.items()}


def figure_coverage(data, out_dir, name="figure_variant_coverage"):
    import matplotlib.pyplot as plt  # noqa: F401  (optional extra `analysis`)
    from analysis import style
    rows = coverage_rows(data)
    keys = ["All three", "Two databases", "CARD only", "NCBI only", "PointFinder only"]
    colours = {"All three": "#000000", "Two databases": style.GREY, "NCBI only": "#0072B2", "CARD only": "#009E73",
               "PointFinder only": "#CC79A7"}
    order = sorted(rows, key=lambda g: (-sum(rows[g].values()), g))
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / f"source_data_{name}.csv", "w") as fh:
        fh.write("organism_group," + ",".join(k.lower().replace(" ", "_") for k in keys) + ",total\n")
        for g in order:
            fh.write(f"{g}," + ",".join(str(rows[g].get(k, 0)) for k in keys) + f",{sum(rows[g].values())}\n")
    height = 24 + 4.2 * len(order)
    fig = style.new_figure(180, height)
    ax = style.axes_mm(fig, 48, 12, 110, height - 26)
    for y, g in enumerate(order):
        total, left = sum(rows[g].values()), 0.0
        for k in keys:
            share = rows[g].get(k, 0) / total
            if share:
                ax.barh(y, share, left=left, color=colours[k], height=0.7, linewidth=0)
            left += share
        ax.text(1.01, y, f"{total:,}", va="center", ha="left", fontsize=7)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([display(g) for g in order], fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("Share of the organism's mutations, by the databases that list them")
    ax.text(1.01, -0.9, "Mutations", ha="left", va="center", fontsize=7)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colours[k]) for k in keys]
    fig.legend(handles, keys, loc="upper left", bbox_to_anchor=(46 / 180, 1 - 1 / height), frameon=False, handlelength=1,
               ncol=len(keys), columnspacing=1.2)
    style.save_figure(fig, out_dir / name)
    return fig


def figure_classes(data, out_dir, name="figure_variant_classes"):
    import matplotlib.pyplot as plt
    from analysis import style
    stats = statistics(data)
    pairs = [f"{a}>{b}" for a in DATABASES for b in DATABASES if a != b]
    colours = dict(zip(CLASS_ORDER, ["#000000", "#009E73", "#E69F00", "#F0E442", "#D55E00", "#56B4E9", "#CC79A7", "#BBBBBB"]))
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / f"source_data_{name}.csv", "w") as fh:
        fh.write("pair," + ",".join(CLASS_ORDER) + ",organism_not_covered\n")
        for p in pairs:
            c = stats["pairs"][p]
            fh.write(p + "," + ",".join(str(c.get(k, 0)) for k in CLASS_ORDER) + f",{c.get('organism_not_covered', 0)}\n")
    fig = style.new_figure(180, 82)
    ax = style.axes_mm(fig, 40, 32, 118, 44)
    for y, p in enumerate(pairs):
        c = stats["pairs"][p]
        total, left = sum(c.get(k, 0) for k in CLASS_ORDER), 0.0
        for k in CLASS_ORDER:
            share = c.get(k, 0) / total if total else 0
            if share:
                ax.barh(y, share, left=left, color=colours[k], height=0.7, linewidth=0)
            left += share
        ax.text(1.01, y, f"{total:,}", va="center", ha="left", fontsize=7)
    a_b = [p.split(">") for p in pairs]
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels([f"{LABELS[a]} against {LABELS[b]}" for a, b in a_b])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("Single mutations of organisms the other database covers")
    ax.text(1.01, -0.9, "Mutations", ha="left", va="center", fontsize=7)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colours[k]) for k in CLASS_ORDER]
    fig.legend(handles, [CLASS_LABELS[k] for k in CLASS_ORDER], loc="lower left", bbox_to_anchor=(4 / 180, 0.0),
               frameon=False, handlelength=1, fontsize=7, ncol=3, columnspacing=1.5)
    style.save_figure(fig, out_dir / name)
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--variants", default="chamredb/data/variants.json")
    parser.add_argument("--out", default="docs")
    args = parser.parse_args(argv)
    data = json.loads(Path(args.variants).read_text())
    out = Path(args.out)
    (out / "variants").mkdir(parents=True, exist_ok=True)
    stats = statistics(data)
    (out / "variants" / "statistics.json").write_text(json.dumps(stats, sort_keys=True, indent=2) + "\n")
    import matplotlib.pyplot as plt
    plt.close(figure_coverage(data, out / "images"))
    plt.close(figure_classes(data, out / "images"))
    print(f"wrote {out}/variants/statistics.json and two figures for {len(data['records'])} mutation records")


if __name__ == "__main__":
    main()
