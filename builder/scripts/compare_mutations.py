"""Build the mutation tables: read CARD, NCBI and PointFinder mutations, match loci, place every change on
its locus anchor, classify every record against the other databases, and write

  variants.json            all of it (records, loci, matches, clusters, exclusions, summary); bundled in the package
  variant_matches.tsv      one row per (record, other database)
  variant_comparison.tsv   one row per mutation cluster, with a status per database
  offset_review.tsv        possible_same_mutation_offset rows, for review (reviewed_equivalences.tsv)
  mutation_report.md       characterisation and comparison report

Usage (from builder/scripts): python compare_mutations.py --card data/mutations/card.json --ncbi data/mutations/ncbi
  --ncbi-proteins data/db_fastas/ncbi.protein.raw.fasta --pointfinder data/mutations/pointfinder_db
  --tables ../variants --out out_dir [--manifest data_manifest.json] [--release YYYY-MM-DD]
Exits 1 if an invariant does not hold (the outputs are still written, for inspection).
"""
import argparse
import csv
import json
import os
import sys
from pathlib import Path

from functions import mutation_compare as mc
from functions import mutation_invariants
from functions import mutation_report
from functions import mutation_sources as ms

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--card", required=True, help="card.json from the pinned CARD release")
parser.add_argument("--ncbi", required=True, help="directory with AMRProt-mutation.tsv and AMR_DNA-*.fa/.tsv")
parser.add_argument("--ncbi-proteins", required=True, help="AMRProt.fa")
parser.add_argument("--pointfinder", required=True, help="pointfinder_db checkout")
parser.add_argument("--tables", required=True, help="directory with taxa.tsv, drug_classes.tsv, gene_aliases.tsv, reviewed_equivalences.tsv")
parser.add_argument("--out", required=True)
parser.add_argument("--manifest", help="data_manifest.json, to record the input releases")
parser.add_argument("--release", default="", help="data release label written into variants.json")
args = parser.parse_args()

out = Path(args.out)
out.mkdir(parents=True, exist_ok=True)
tables = Path(args.tables)

card = ms.read_card(args.card)
ncbi = ms.read_ncbi(os.path.join(args.ncbi, "AMRProt-mutation.tsv"), args.ncbi_proteins, args.ncbi)
pointfinder = ms.read_pointfinder(args.pointfinder)
sources = {"card": card, "ncbi": ncbi, "pointfinder": pointfinder}

records = card["records"] + ncbi["records"] + pointfinder["records"]
references = {}
for s in sources.values():
    references.update(s["references"])
excluded = card["excluded"] + ncbi["excluded"] + pointfinder["excluded"]

with open(os.path.join(args.ncbi, "taxgroup.tsv")) as fh:
    amrfinder_organisms = {l.split("\t")[0] for l in fh if l.strip() and not l.startswith("#")}
mc.assign_groups(records, mc.read_table(tables / "taxa.tsv"), amrfinder_organisms)
mc.normalise_classes(records, mc.read_table(tables / "drug_classes.tsv"))
aliases = {r["label"]: r["locus_gene"] for r in mc.read_table(tables / "gene_aliases.tsv")}
loci, joined = mc.match_loci(records, references, aliases)
rrna_strands = mc.orient_rrna(records, loci, references)
maps = mc.lift_maps(loci, references)
mc.place(records, loci, references, maps)
reviewed = mc.read_table(tables / "reviewed_equivalences.tsv")
matches = mc.classify(records, loci, references, reviewed)
clusters = mc.cluster(records, matches)
records_by_id = {r["record_id"]: r for r in records}

inputs = {}
if args.manifest:
    for f in json.loads(Path(args.manifest).read_text())["files"]:
        inputs.setdefault(f["database"], {})[f["role"]] = {"release": f["release"], "url": f["url"], "sha256": f["sha256"]}

data = {
    "schema": "chamredb-variants/1",
    "data_release": args.release,
    "inputs": inputs,
    "card_version": card.get("version"),
    "settings": {"anchor_order": list(mc.DATABASES), "offset_window": mc.WINDOW,
                 "locus_sequence_identity": mc.LOCUS_SEQUENCE_IDENTITY, "locus_sequence_coverage": mc.LOCUS_SEQUENCE_COVERAGE,
                 "lift_min_identity": mc.LIFT_MIN_IDENTITY,
                 "lift_min_coverage": mc.LIFT_MIN_COVERAGE, "rrna_antisense_margin": mc.RRNA_ANTISENSE_MARGIN, "related_organism": "same family (taxa.tsv group_family)",
                 "organism_groups": "AMRFinderPlus organisms, else PointFinder species, else CARD species"},
    "sources": {db: {"rows_read": s["rows_read"], "rows_parsed": s["rows_parsed"],
                     "records": len(s["records"])} for db, s in sources.items()},
    "loci": loci,
    "loci_joined_by_sequence": joined,
    "records": records,
    "matches": matches,
    "clusters": clusters,
    "excluded": excluded,
    "rrna_strands": rrna_strands,
    "self_check_offsets": mc.self_check_offsets(records, references),
    "systematic_offsets": mc.systematic_offsets(matches, records_by_id),
}
violations = mutation_invariants.check_variants(data)
data["summary"] = mutation_report.summary(data)
(out / "variants.json").write_text(json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n")


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def describe(r):
    return " + ".join(c["hgvs"] or f"[{c['native_hgvs']}]" for c in r["components"])


write_tsv(out / "variant_matches.tsv",
          ["record_id", "database", "group", "loci", "native_notation", "hgvs_anchor", "other_database", "class",
           "matched_notation", "matched_record_ids", "offset", "distance", "phenotype", "note"],
          [[m["record_id"], m["database"], records_by_id[m["record_id"]]["taxon_group"],
            ",".join(records_by_id[m["record_id"]]["locus_ids"]), records_by_id[m["record_id"]]["native_notation"],
            describe(records_by_id[m["record_id"]]), m["other_database"], m["class"],
            "; ".join(records_by_id[h]["native_notation"] for h in m["matched"]), ";".join(m["matched"]),
            "" if m["offset"] is None else m["offset"], "" if m["distance"] is None else m["distance"],
            m["phenotype"] or "", m["note"]] for m in matches])

write_tsv(out / "variant_comparison.tsv",
          ["cluster_id", "group", "loci", "kind", "hgvs_anchor", "databases_listing"] + [f"{db}_status" for db in mc.DATABASES]
          + [f"{db}_notation" for db in mc.DATABASES] + ["drug_classes"],
          [[c["cluster_id"], c["group"], ",".join(c["loci"]), c["kind"], c["hgvs"], c["databases"]]
           + [c["status"][db] for db in mc.DATABASES]
           + ["; ".join(records_by_id[m]["native_notation"] for m in c["members"] if records_by_id[m]["database"] == db) for db in mc.DATABASES]
           + [";".join(sorted({x for m in c["members"] for x in records_by_id[m]["drug_classes"]}))] for c in clusters])

review = [m for m in matches if m["class"] == "possible_same_mutation_offset" and m["database"] < m["other_database"]]
write_tsv(out / "offset_review.tsv",
          ["record_a", "record_b", "group", "locus", "notation_a", "notation_b", "hgvs_anchor_a", "hgvs_anchor_b", "offset",
           "self_check_a", "self_check_b"],
          [[m["record_id"], h, records_by_id[m["record_id"]]["taxon_group"], records_by_id[m["record_id"]]["locus_ids"][0],
            records_by_id[m["record_id"]]["native_notation"], records_by_id[h]["native_notation"],
            describe(records_by_id[m["record_id"]]), describe(records_by_id[h]), "" if m["offset"] is None else m["offset"],
            records_by_id[m["record_id"]]["components"][0]["ref_check"], records_by_id[h]["components"][0]["ref_check"]]
           for m in review for h in m["matched"]])

(out / "mutation_report.md").write_text(mutation_report.markdown(data, violations))
print(f"records {len(records)}, loci {len(loci)}, matches {len(matches)}, clusters {len(clusters)}, excluded {len(excluded)}")
if violations:
    print("mutation invariants do not hold:", {k: v[:3] for k, v in violations.items()}, file=sys.stderr)
    sys.exit(1)
print("mutation invariants hold")
