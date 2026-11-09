"""Regenerate builder/tests/fixtures/data from the public 2021 files committed in chamredb/data.

Note: decoys are chosen using the shipped chamredb/data/graph.json (proteins with no edges there), and the
CARD raw FASTA headers are reconstructed (synthetic accession) because only parsed CARD FASTAs are committed.

Provenance only: the fixtures are committed, this script documents how they were chosen.
Standard library only. Run from the repository root: python3 builder/tests/make_fixtures.py
"""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "chamredb" / "data"
OUT = Path(__file__).resolve().parent / "fixtures" / "data"

# Designed cases (all present in the committed 2021 FASTAs):
#   dfrB1, dfrB2, dfrB3, dfrB6, tet(37): RBH between all three databases
#   card dfrB7 (ARO:3004556) / ncbi dfrB8 / resfinder dfrB8: RBH with a database name mismatch
#   resfinder dfrB7: only one-way hits (no RBH partner)
#   card DfrB9 / ncbi dfrB9: RBH with each other, only one-way hits to resfinder (no resfinder dfrB9)
#   mecI, vanU-G: RBH between card and ncbi only
#   card ARO:3002718 (215 aa) / ncbi WP_017900990.1 (214 aa): the NCBI protein is the CARD protein
#   without its first residue, so coverage is below 1 (214/215) and checkable by hand
CARD = ["ARO:3002864", "ARO:3003021", "ARO:3003022", "ARO:3003023", "ARO:3002871",
        "ARO:3004556", "ARO:3005350", "ARO:3000124", "ARO:3004253", "ARO:3002718"]
NCBI = ["WP_063844476.1", "WP_063844477.1", "WP_010890144.1", "WP_063844479.1", "WP_063855713.1",
        "WP_063844480.1", "WP_071846200.1", "WP_000369216.1", "WP_021418929.1", "WP_017900990.1"]
RESFINDER = ["dfrB1", "dfrB2", "dfrB3", "dfrB6", "tet(37)", "dfrB8", "dfrB7"]
N_DECOYS = {"card": 10, "ncbi": 8, "resfinder": 11}
MAX_ALLELES = 2


def read_fasta(path):
    records, key = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            key = line[1:]
            records[key] = []
        else:
            records[key].append(line)
    return {k: "".join(v) for k, v in records.items()}


def write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def main():
    for sub in ("db_metadata", "db_fastas"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    graph = json.load(open(SRC / "graph.json"))
    in_graph = {n["id"] for n in graph["nodes"]}

    # ---- CARD
    card_seq = {k: v for k, v in read_fasta(SRC / "db_fastas" / "card.protein.fasta").items()}
    card_meta = json.load(open(SRC / "db_metadata" / "card.metadata.json"))
    decoys = [i for i in sorted(card_seq, key=lambda i: (len(card_seq[i]), i))
              if f"card:{i}" not in in_graph and len(card_seq[i]) >= 80 and i in card_meta][: N_DECOYS["card"]]
    card_ids = CARD + decoys
    write_fasta(
        OUT / "protein_fasta_protein_homolog_model.fasta",
        # reconstructed in the current CARD header layout; the parser only reads the ARO field
        [(f"gb|FIXTURE{n:03d}.1|{i}|{card_meta[i]['name'].replace(' ', '_')}", card_seq[i])
         for n, i in enumerate(card_ids, 1)],
    )
    write_obo_subset(card_ids, OUT / "db_metadata" / "card.metadata.obo")

    # ---- NCBI
    raw = read_fasta(SRC / "db_fastas" / "ncbi.protein.raw.fasta")
    by_acc = {h.split("|")[1]: (h, s) for h, s in raw.items()}
    with open(SRC / "db_metadata" / "ncbi.metadata.tsv") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header_fields = list(reader.fieldnames or [])
        rows = list(reader)
    core = [r for r in rows if r["scope"] == "core" and r["type"] == "AMR" and r["subtype"] == "AMR"
            and r["refseq_protein_accession"] in by_acc]
    in_ncbi_graph = lambda a: f"ncbi:{a}" in in_graph
    decoy_rows = [r for r in sorted(core, key=lambda r: (len(by_acc[r["refseq_protein_accession"]][1]),
                                                         r["refseq_protein_accession"]))
                  if not in_ncbi_graph(r["refseq_protein_accession"])
                  and len(by_acc[r["refseq_protein_accession"]][1]) >= 80
                  and r["refseq_protein_accession"] not in NCBI][: N_DECOYS["ncbi"]]
    # duplicate case: several protein accessions sharing one name (empty allele, same gene_family);
    # the parsers key on allele combined with gene_family and keep only the first
    counts = {}
    for r in core:
        counts.setdefault(r["allele"] or r["gene_family"], []).append(r)
    dup = next(v for k, v in sorted(counts.items())
               if k and len(v) > 1 and all(not in_ncbi_graph(x["refseq_protein_accession"]) for x in v)
               and all(len(by_acc[x["refseq_protein_accession"]][1]) < 400 for x in v))[:3]
    # non-core case: a plus-scope row that the parsers must drop
    plus = [r for r in rows if r["scope"] != "core" and r["refseq_protein_accession"] in by_acc][:2]
    chosen_rows = [r for r in rows if r["refseq_protein_accession"] in NCBI] + decoy_rows + dup + plus
    seen, final_rows = set(), []
    for r in chosen_rows:
        key = (r["allele"], r["gene_family"], r["refseq_protein_accession"])
        if key not in seen:
            seen.add(key)
            final_rows.append(r)
    with open(OUT / "db_metadata" / "ncbi.metadata.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header_fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(final_rows)
    write_fasta(OUT / "db_fastas" / "ncbi.protein.raw.fasta",
                [by_acc[r["refseq_protein_accession"]] for r in final_rows])

    # ---- ResFinder
    nucl = read_fasta(SRC / "db_fastas" / "resfinder.nucl.fasta")
    genes = {}
    for h, s in nucl.items():
        genes.setdefault(h.split("_")[0], []).append((h, s))
    graph_genes = {i.split(":", 1)[1] for i in in_graph if i.startswith("resfinder:")}
    res_decoys = [g for g in sorted(genes, key=lambda g: (len(genes[g][0][1]), g))
                  if g not in graph_genes and len(genes[g][0][1]) >= 300 and g not in RESFINDER][: N_DECOYS["resfinder"]]
    res_genes = RESFINDER + res_decoys
    write_fasta(OUT / "db_fastas" / "resfinder.nucl.fasta",
                [rec for g in res_genes for rec in genes[g][:MAX_ALLELES]])
    meta_lines = open(SRC / "db_metadata" / "resfinder.metadata.tsv").read().splitlines()
    keep = [meta_lines[0]] + [ln for ln in meta_lines[1:] if ln.split("\t")[0].split("_")[0] in set(res_genes)]
    (OUT / "db_metadata" / "resfinder.metadata.tsv").write_text("\n".join(keep) + "\n")
    print(f"card {len(card_ids)}, ncbi {len(final_rows)} rows, resfinder {len(res_genes)} genes")


def write_obo_subset(ids, path):
    text = open(SRC / "db_metadata" / "card.metadata.obo").read()
    head, *stanzas = re.split(r"\n(?=\[(?:Term|Typedef)\])", text)
    by_id = {}
    for st in stanzas:
        m = re.search(r"^id: (\S+)", st, re.M)
        if m:
            by_id[m.group(1)] = st
    need, todo = set(), list(ids)
    while todo:
        i = todo.pop()
        if i in need or i not in by_id:
            continue
        need.add(i)
        for ref in re.findall(r"^(?:is_a|relationship): (?:\S+ )?(ARO:\d+)", by_id[i], re.M):
            todo.append(ref)
    typedefs = [st for st in stanzas if st.startswith("[Typedef]")]
    terms = [by_id[i] for i in sorted(need) if by_id[i].startswith("[Term]")]
    path.write_text(head + "\n" + "\n".join(s.rstrip("\n") + "\n" for s in typedefs + terms))


if __name__ == "__main__":
    main()
