"""The `chamredb variants` command."""
import csv
import json
import re

from rich import box
from rich.console import Console
from rich.markup import escape
from rich.table import Table

from chamredb.variants import CLASS_TEXT, DATABASES, NAMES, Variants, normalise_mutation, split_query

SOFTWARE_DATABASE = {"amrfinderplus": "ncbi", "rgi": "card", "pointfinder": "pointfinder", "resfinder": "pointfinder",
                     "resfinder 4": "pointfinder"}


def _console():
    return Console(highlight=False, soft_wrap=False)


def print_record(v, r, console):
    head = Table(box=box.SIMPLE, show_header=False, pad_edge=False)
    head.add_column(style="bold")
    head.add_column()
    head.add_row("database", escape(f"{NAMES[r['database']]}: {r['native_notation']}"))
    head.add_row("reference", escape(r["native_label"]))
    head.add_row("organism", escape(f"{r['organism_label']} (group: {r['taxon_group']})"))
    head.add_row("locus", escape(", ".join(r["locus_ids"])))
    head.add_row("on the anchor", escape(f"{v.describe(r)}  (numbering of {v.anchor_of(r)})"))
    checks = ", ".join(sorted({c["ref_check"] for c in r["components"]}))
    head.add_row("reference check", checks)
    head.add_row("drug classes", escape(", ".join(r["drug_classes"]) or "-"))
    head.add_row("phenotype", escape(r["phenotype_native"] or "-"))
    cluster = v.cluster_of.get(r["record_id"])
    if cluster:
        head.add_row("cluster", cluster["cluster_id"])
    console.print(head)
    t = Table(box=box.SIMPLE_HEAD)
    t.add_column("database")
    t.add_column("result")
    t.add_column("their notation")
    t.add_column("note")
    for db, m in v.comparison(r).items():
        result = CLASS_TEXT[m["class"]]
        if m["offset"] is not None:
            result += f" (offset {m['offset']:+d})"
        if m["distance"] is not None:
            result += f" (nearest {m['distance']} positions away)"
        their = "; ".join(f"{h['native_notation']} [{h['taxon_group']}]" if h["taxon_group"] != r["taxon_group"]
                          else h["native_notation"] for h in m["matched"][:6])
        if len(m["matched"]) > 6:
            their += f" (+{len(m['matched']) - 6} more)"
        note = m["note"]
        if m["phenotype"] and m["phenotype"] != "same_classes":
            note = "; ".join(x for x in (note, f"drug classes: {m['phenotype'].replace('_', ' ')}") if x)
        t.add_row(NAMES[db], escape(result), escape(their), escape(note))
    console.print(t)


def query_single(query, database, taxon, variants_path=None, console=None):
    console = console or _console()
    v = Variants.load(variants_path)
    console.print(f"data release: {v.release}")
    found = v.find(query, database=database, taxon=taxon)
    if not found and "," in query:
        found = v.find_combination([split_query(q) for q in query.split(",")], database=database, taxon=taxon)
    if not found:
        where = f" in {NAMES[database]}" if database else ""
        console.print(f"no mutation matching '{query}'{where}" + (f" for '{taxon}'" if taxon else ""))
        return 1
    if len({r["taxon_group"] for r in found}) > 1 and not taxon:
        console.print(f"'{query}' is listed for {len({r['taxon_group'] for r in found})} organism groups; "
                      f"use --taxon to choose one. Showing all.")
    for r in found:
        print_record(v, r, console)
    return 0


def query_locus(locus, variants_path=None, console=None):
    console = console or _console()
    v = Variants.load(variants_path)
    console.print(f"data release: {v.release}")
    ids = v.locus_ids(locus)
    if not ids:
        console.print(f"no locus matching '{locus}' (form: Organism:gene, e.g. Escherichia:gyrA)")
        return 1
    for lid in ids:
        l = v.loci[lid]
        console.print(f"[bold]{escape(lid)}[/bold]  anchor: {NAMES[l['anchor_database']]} {escape(l['anchor_accession'] or '(no sequence)')}; "
                      "references: " + escape(", ".join(f"{NAMES[x['database']]} {x.get('accession') or x['reference']}" for x in l["references"])))
        t = Table(box=box.SIMPLE_HEAD)
        t.add_column("mutation (anchor numbering)")
        for db in DATABASES:
            t.add_column(NAMES[db])
        for c in v.locus_clusters(lid):
            cells = []
            for db in DATABASES:
                status = c["status"][db]
                cells.append(escape("; ".join(v.notations(c, db))) if status == "listed" else f"[dim]{status.replace('_', ' ')}[/dim]")
            t.add_row(escape(c["hgvs"]), *cells)
        console.print(t)
    return 0


TSV_HEADER = ["input_file", "query_database", "query", "found", "record_id", "group", "loci", "hgvs_anchor", "anchor"] + \
    [f"{db}_{x}" for db in DATABASES for x in ("status", "notation")] + ["drug_classes", "data_release"]


def _rows_for(v, item):
    database, query, taxon = item["database"], item["id"], item.get("taxon")
    found = v.find(query, database=database, taxon=taxon, accession=item.get("accession"))
    if not found and item.get("accession"):
        found = v.find(query, database=database, taxon=taxon)
    if not found and "," in query:
        found = v.find_combination([split_query(q) for q in query.split(",")], database=database, taxon=taxon)
    base = [item.get("file", ""), database or "", query]
    if not found:
        return [base + ["no"] + [""] * (len(TSV_HEADER) - len(base) - 2) + [v.release]]
    rows = []
    for r in found:
        cells = []
        for db in DATABASES:
            if db == r["database"]:
                cells += ["listed", r["native_notation"]]
            else:
                m = v.comparison(r)[db]
                cells += [m["class"], "; ".join(h["native_notation"] for h in m["matched"])]
        rows.append(base + ["yes", r["record_id"], r["taxon_group"], ",".join(r["locus_ids"]), v.describe(r), v.anchor_of(r)]
                    + cells + [",".join(r["drug_classes"]), v.release])
    return rows


def query_multiple(items, outfile_path, variants_path=None):
    v = Variants.load(variants_path)
    with open(outfile_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(TSV_HEADER)
        for item in items:
            for row in _rows_for(v, item):
                w.writerow(row)


def parse_hamronization_mutations(path):
    """Mutation rows of a hAMRonization summary from AMRFinderPlus, RGI and PointFinder."""
    items = []
    for row in json.load(open(path)):
        software = (row.get("analysis_software_name") or "").lower()
        database = SOFTWARE_DATABASE.get(software)
        if database is None:
            continue
        gene = row.get("gene_symbol") or ""
        muts = [m for m in (row.get("amino_acid_mutation"), row.get("nucleotide_mutation")) if m]
        variation = (row.get("genetic_variation_type") or "").lower()
        if not muts:
            g, m = split_query(gene)
            if g and re.search(r"[A-Za-z*?]-?\d+", m or "") and (variation.endswith("variant_detected") or software == "amrfinderplus"):
                gene, muts = g, [m]
        if not muts:
            continue
        if gene and "_" in gene and split_query(gene)[0]:
            g, m = split_query(gene)
            if normalise_mutation(m) in {normalise_mutation(x) for x in muts}:
                gene = g
        for text in muts:
            parts = [p.strip() for p in text.split(",") if p.strip()]
            query = ",".join(f"{gene} {p}" for p in parts) if gene else ",".join(parts)
            items.append({"file": row.get("input_file_name", ""), "database": database, "id": query,
                          "accession": row.get("reference_accession") or None})
    return items

