"""Read the mutations of CARD (card.json), NCBI (AMRProt-mutation.tsv, AMR_DNA-*) and PointFinder (pointfinder_db).

Each reader returns {records, references, excluded, rows_read, rows_parsed}. A record has one change per
component (several for a combination); each component also carries gene_label, reference (a key of references) and
index (position of its first residue in the reference sequence, or None). Rows that cannot be read go to excluded
with the reason, so rows_read == rows_parsed + len(excluded).
"""
import csv
import glob
import json
import os
import re

from functions import mutation_notation as mn
from functions import validation

CARD_MODEL_TYPES = ("protein variant model", "rRNA gene variant model", "protein overexpression model")
CARD_SINGLE_PARAMS = {
    "single resistance variant", "nonsense mutation", "frameshift mutation",
    "insertion mutation from peptide sequence", "deletion mutation from peptide sequence",
    "duplication of peptide sequence", "nucleotide substitution in promoter region",
    "synonymous substitution for alternative promoter",
}
CARD_COMBINATION_PARAMS = {"multiple resistance variants"}
CARD_CODEPENDENT_PARAMS = {"co-dependent single resistance variant"}
CARD_NON_MUTATION_PARAMS = {"BLASTP bit-score", "BLASTN bit-score"}
CARD_PARAM_META = {"param_type", "param_description", "param_type_id"}

POINTFINDER_EXCLUDED_SPECIES = {"plasmodium_falciparum": "not a bacterium"}


def _record(database, record_id, notation, label, organism, components, phenotype, classes, sources):
    return {
        "database": database, "record_id": record_id, "native_notation": notation, "native_label": label,
        "organism_label": organism, "components": components, "phenotype_native": phenotype,
        "drug_classes_native": sorted(set(classes)), "sources": sorted(set(sources)),
    }


def _component(c, gene, reference, index):
    out = dict(c)
    out.update({"gene_label": gene, "reference": reference, "index": index})
    return out


def merge_duplicates(records):
    """Merge records listed more than once (per evidence source in CARD, per phenotype in PointFinder)."""
    merged = {}
    for r in records:
        m = merged.get(r["record_id"])
        if m is None:
            merged[r["record_id"]] = r
            continue
        m["sources"] = sorted(set(m["sources"]) | set(r["sources"]))
        m["drug_classes_native"] = sorted(set(m["drug_classes_native"]) | set(r["drug_classes_native"]))
        if r["phenotype_native"] and r["phenotype_native"] not in m["phenotype_native"].split("; "):
            m["phenotype_native"] = "; ".join(p for p in (m["phenotype_native"], r["phenotype_native"]) if p)
    return sorted(merged.values(), key=lambda r: r["record_id"])


def record_id(database, components):
    parts = sorted(f"{c['reference']}:{mn.hgvs(c)}" for c in components)
    return f"{database}|" + "+".join(parts)


def card_gene(model):
    """Gene symbol from the CARD short name ('Ecol_gyrA_FLO' -> 'gyrA'), else the model name."""
    short = model.get("CARD_short_name") or ""
    parts = short.split("_")
    if len(parts) >= 3:
        return "_".join(parts[1:-1])
    if len(parts) == 2:
        return parts[1]
    return short or model["model_name"]


def card_drug_classes(model):
    return [c["category_aro_name"] for c in (model.get("ARO_category") or {}).values()
            if c.get("category_aro_class_name") == "Drug Class"]


def read_card(card_json_path):
    card = json.load(open(card_json_path))
    version = card.get("_version")
    models = {k: v for k, v in card.items() if not k.startswith("_")}
    by_aro_id = {m["ARO_id"]: m for m in models.values()}
    records, references, excluded = [], {}, []
    n = {"read": 0, "parsed": 0}

    def refs_for(model):
        seqs = list((model.get("model_sequences") or {}).get("sequence", {}).values())
        if len(seqs) != 1:
            raise ValueError(f"CARD model {model['model_id']}: expected one sequence, found {len(seqs)}")
        s = seqs[0]
        tax = (s.get("NCBI_taxonomy") or {}).get("NCBI_taxonomy_name", "")
        gene = card_gene(model)
        aro = "ARO:" + model["ARO_accession"]
        out = {}
        if "protein_sequence" in s and s["protein_sequence"].get("sequence"):
            key = f"card|{aro}|protein"
            references[key] = {"database": "card", "accession": s["protein_sequence"]["accession"], "molecule": "protein",
                               "gene_label": gene, "organism_label": tax, "sequence": s["protein_sequence"]["sequence"]}
            out["protein"] = key
        if "dna_sequence" in s and s["dna_sequence"].get("sequence"):
            mol = "rna" if model["model_type"] == "rRNA gene variant model" else "dna"
            key = f"card|{aro}|{mol}"
            references[key] = {"database": "card", "accession": s["dna_sequence"]["accession"], "molecule": mol,
                               "gene_label": gene, "organism_label": tax, "sequence": s["dna_sequence"]["sequence"]}
            out[mol] = key
        key = f"card|{aro}|promoter"   # no upstream sequence in CARD: positions cannot be checked
        references[key] = {"database": "card", "accession": None, "molecule": "promoter", "gene_label": gene,
                           "organism_label": tax, "sequence": ""}
        out["promoter"] = key
        return out, gene, tax

    def comp(token, model, refs, gene):
        molecule = "rna" if model["model_type"] == "rRNA gene variant model" else "protein"
        c = mn.parse_card(token, molecule)
        ref = refs.get(c["molecule"])
        if ref is None:
            raise mn.NotationError(f"no CARD {c['molecule']} sequence for {token!r}")
        index = None if c["molecule"] == "promoter" else c["start"]
        return _component(c, gene, ref, index)

    for model in sorted(models.values(), key=lambda m: int(m["model_id"])):
        if model["model_type"] not in CARD_MODEL_TYPES:
            continue
        refs, gene, tax = refs_for(model)
        classes = card_drug_classes(model)
        aro = "ARO:" + model["ARO_accession"]
        for key, param in sorted(model.get("model_param", {}).items()):
            ptype = param.get("param_type") or ("single resistance variant" if key == "snp" else None)
            if ptype in CARD_NON_MUTATION_PARAMS:
                continue
            values = {}
            for source, entries in param.items():
                if source in CARD_PARAM_META or not isinstance(entries, dict):
                    continue
                for v in entries.values():
                    values.setdefault(v, set()).add(source)
            for value, sources in sorted(values.items()):
                n["read"] += 1
                row = f"{aro} {ptype}: {value}"
                try:
                    if ptype in CARD_SINGLE_PARAMS:
                        components = [comp(value, model, refs, gene)]
                    elif ptype in CARD_COMBINATION_PARAMS:
                        components = [comp(t, model, refs, gene) for t in value.split(",")]
                    elif ptype in CARD_CODEPENDENT_PARAMS:
                        components = []
                        for aro_id, tokens in mn.parse_card_codependent(value):
                            other = by_aro_id.get(aro_id)
                            if other is None:
                                raise mn.NotationError(f"co-dependent variant names unknown ARO id {aro_id}")
                            orefs, ogene, _ = refs_for(other)
                            components += [comp(t, other, orefs, ogene) for t in tokens]
                    else:
                        excluded.append({"database": "card", "row": row, "reason": f"parameter type not a mutation: {ptype}"})
                        continue
                except mn.NotationError as e:
                    excluded.append({"database": "card", "row": row, "reason": str(e)})
                    continue
                records.append(_record("card", record_id("card", components), value, f"{aro} {model['model_name']}",
                                       tax, components, ptype, classes, sources))
                n["parsed"] += 1
    return {"records": merge_duplicates(records), "references": references, "excluded": excluded,
            "rows_read": n["read"], "rows_parsed": n["parsed"], "version": version}


def read_fasta(path):
    seqs, key = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                key = line[1:].split()[0]
                if key in seqs:
                    key = None   # duplicate header: keep the first
                    continue
                seqs[key] = []
            elif key is not None:
                seqs[key].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def _ncbi_classes(cls, subclass):
    return [c for c in cls.split("/") if c]


def _ncbi_notes(symbol, reported, c):
    """The reported symbol if it differs, and the symbol's own number if it is not the sequence position."""
    notes = [f"reported as {reported}"] if reported != symbol else []
    written = mn.symbol_position(_ncbi_gene(symbol)[1])
    if c["molecule"] != "promoter" and written is not None and written != c["start"]:
        notes.append(f"symbol numbering {written}, sequence position {c['start']}")
    return notes


def _ncbi_gene(symbol):
    return symbol.rsplit("_", 1)


def read_ncbi(protein_mutation_tsv, protein_fasta, dna_dir):
    records, references, excluded = [], {}, []
    n = {"read": 0, "parsed": 0}
    proteins = {}
    for header, seq in read_fasta(protein_fasta).items():
        acc = header.split("|")[0]
        proteins.setdefault(acc, seq)

    with open(protein_mutation_tsv) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]
    validation.require_columns([h.lstrip("#") for h in header],
                               ["taxgroup", "accession_version", "mutation_position", "standard_mutation_symbol",
                                "reported_mutation_symbol", "class", "subclass"], protein_mutation_tsv)
    for row in rows[1:]:
        if not row or row[0].startswith("#"):
            continue
        n["read"] += 1
        taxgroup, acc, pos, symbol, reported, cls, subclass = row[:7]
        gene, mut = _ncbi_gene(symbol)
        seq = proteins.get(acc)
        if seq is None:
            excluded.append({"database": "ncbi", "row": "\t".join(row[:4]), "reason": f"protein {acc} not in AMRProt"})
            continue
        key = f"ncbi|{acc}|protein"
        references[key] = {"database": "ncbi", "accession": acc, "molecule": "protein", "gene_label": gene,
                           "organism_label": taxgroup, "sequence": seq}
        try:
            c, index = mn.parse_ncbi(mut, int(pos), "protein")
        except mn.NotationError as e:
            excluded.append({"database": "ncbi", "row": "\t".join(row[:4]), "reason": str(e)})
            continue
        components = [_component(c, gene, key, index)]
        records.append(_record("ncbi", record_id("ncbi", components), symbol, acc, taxgroup, components,
                               f"{cls} ({subclass})", _ncbi_classes(cls, subclass), _ncbi_notes(symbol, reported, c)))
        n["parsed"] += 1

    for tsv in sorted(glob.glob(os.path.join(dna_dir, "AMR_DNA-*.tsv"))):
        taxgroup = os.path.basename(tsv)[len("AMR_DNA-"):-len(".tsv")]
        seqs = read_fasta(tsv[:-4] + ".fa")
        with open(tsv) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        for row in rows[1:]:
            if not row or row[0].startswith("#"):
                continue
            n["read"] += 1
            region, pos, symbol, reported, cls, subclass = row[:6]
            acc, gene_name, rest = region.split("@")
            gene = rest.split(":")[0]
            molecule = "rna" if "ribosomal_RNA" in gene_name else "promoter" if "promoter" in gene_name else "dna"
            seq = seqs.get(region)
            if seq is None:
                excluded.append({"database": "ncbi", "row": "\t".join(row[:3]), "reason": f"{region} not in {taxgroup}.fa"})
                continue
            key = f"ncbi|{region}|{molecule}"
            references[key] = {"database": "ncbi", "accession": acc + ":" + rest.split(":")[1], "molecule": molecule,
                               "gene_label": gene, "organism_label": taxgroup, "sequence": seq}
            _, mut = _ncbi_gene(symbol)
            try:
                c, index = mn.parse_ncbi(mut, int(pos), molecule)
            except mn.NotationError as e:
                excluded.append({"database": "ncbi", "row": "\t".join(row[:3]), "reason": str(e)})
                continue
            components = [_component(c, gene, key, index)]
            records.append(_record("ncbi", record_id("ncbi", components), symbol, region, taxgroup, components,
                                   f"{cls} ({subclass})", _ncbi_classes(cls, subclass), _ncbi_notes(symbol, reported, c)))
            n["parsed"] += 1
    return {"records": merge_duplicates(records), "references": references, "excluded": excluded,
            "rows_read": n["read"], "rows_parsed": n["parsed"]}


_PF_GENE = re.compile(r"^(?P<gene>.+?)_(?P<n>\d+)_(?P<acc>[A-Z]{1,2}_?\d+(?:\.\d+)?)$")
_PF_PROMOTER = re.compile(r"^(?P<gene>.+?)-promoter-size-(?P<size>\d+)bp$")


def pointfinder_gene(accession):
    """'gyrA_1_CP073768.1' -> ('gyrA', None); 'ampC-promoter-size-53bp_1_CP037449.1' -> ('ampC', 53);
    'porA_YP_002344650.1' -> ('porA', None)."""
    m = _PF_GENE.match(accession)
    name = m.group("gene") if m else accession.split("_")[0]
    p = _PF_PROMOTER.match(name)
    if p:
        return p.group("gene"), int(p.group("size"))
    return name, None


def _translate(nucl):
    from Bio.Seq import Seq
    n = len(nucl) - len(nucl) % 3
    return str(Seq(nucl[:n]).translate())


def read_pointfinder(db_dir):
    records, references, excluded = [], {}, []
    n = {"read": 0, "parsed": 0}
    for species in sorted(os.listdir(db_dir)):
        pheno = os.path.join(db_dir, species, "phenotypes.txt")
        if not os.path.isfile(pheno):
            continue
        with open(pheno) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        validation.require_columns([h.strip().lstrip("#") for h in rows[0]],
                                   ["Gene_accession", "TypeGene", "Codon_pos", "Ref_nuc", "Ref_codon", "Res_codon",
                                    "Class", "Phenotype", "Required mutation"], pheno)
        if species in POINTFINDER_EXCLUDED_SPECIES:
            for row in rows[1:]:
                if row and not row[0].startswith("#"):
                    n["read"] += 1
                    excluded.append({"database": "pointfinder", "row": f"{species}: " + "\t".join(row[:4]),
                                     "reason": POINTFINDER_EXCLUDED_SPECIES[species]})
            continue
        fsa = {}
        for path in glob.glob(os.path.join(db_dir, species, "*.fsa")):
            for k, v in read_fasta(path).items():
                fsa[k] = v.upper()
        rna_genes = set()
        rna_path = os.path.join(db_dir, species, "RNA_genes.txt")
        if os.path.isfile(rna_path):
            rna_genes = {l.strip() for l in open(rna_path) if l.strip() and not l.startswith("#")}

        def sequence_of(accession):
            if accession in fsa:
                return fsa[accession]
            stripped = re.sub(r"-size-\d+bp", "", accession)
            return fsa.get(stripped)

        def component_for(accession, c):
            """Attach gene, reference and index to a parsed change of gene `accession`."""
            gene, size = pointfinder_gene(accession)
            file_gene = re.sub(r"_\d+_[^_]+$", "", accession)
            nucl = sequence_of(accession)
            if nucl is None:
                raise mn.NotationError(f"no sequence for {accession}")
            is_rna = file_gene in rna_genes or re.sub(r"-size-\d+bp$", "", file_gene) in rna_genes \
                or gene in ("rrs", "rrl") or gene.startswith(("16S", "23S"))
            if c["molecule"] == "dna" and is_rna:
                c = dict(c, molecule="rna")
            gene_seq = nucl[size:] if size else nucl
            if c["molecule"] == "promoter":
                key = f"pointfinder|{species}/{accession}|promoter"
                references[key] = {"database": "pointfinder", "accession": accession, "molecule": "promoter",
                                   "gene_label": gene, "organism_label": species, "sequence": nucl}
                index = (size or 0) + c["start"] + 1
            elif c["molecule"] == "protein":
                key = f"pointfinder|{species}/{accession}|protein"
                references[key] = {"database": "pointfinder", "accession": accession, "molecule": "protein",
                                   "gene_label": gene, "organism_label": species, "sequence": _translate(gene_seq)}
                index = c["start"]
            else:
                key = f"pointfinder|{species}/{accession}|{c['molecule']}"
                references[key] = {"database": "pointfinder", "accession": accession, "molecule": c["molecule"],
                                   "gene_label": gene, "organism_label": species, "sequence": gene_seq}
                index = c["start"]
            return _component(c, gene, key, index)

        def parse(accession, type_gene, pos, ref_codon, res, ref_nuc):
            _, size = pointfinder_gene(accession)
            return [component_for(accession, c)
                    for c in mn.parse_pointfinder(type_gene, pos, ref_codon, res, ref_nuc, size is not None and int(pos) < 0)]

        by_gene, type_at = {}, {}
        for row in rows[1:]:
            if row and not row[0].startswith("#"):
                by_gene.setdefault(pointfinder_gene(row[0].strip())[0], row[0].strip())
                type_at[(row[0].strip(), row[3].strip())] = row[1].strip()

        for row in rows[1:]:
            if not row or row[0].startswith("#"):
                continue
            n["read"] += 1
            row = row + [""] * (14 - len(row))
            acc = row[0].strip()
            text = f"{species}: " + "\t".join(c.strip() for c in row[:7])
            try:
                singles = parse(acc, row[1].strip(), row[3].strip(), row[5].strip(), row[6].strip(), row[4].strip())
                required = row[12].strip()
                alternatives = mn.parse_pointfinder_required(required) if required else []
            except (mn.NotationError, ValueError) as e:
                excluded.append({"database": "pointfinder", "row": text, "reason": str(e)})
                continue
            classes = [c.strip() for c in row[7].split(",") if c.strip()]
            phenotype = row[8].strip()
            gene = pointfinder_gene(acc)[0]
            if not alternatives:
                for c in singles:
                    records.append(_record("pointfinder", record_id("pointfinder", [c]), f"{gene}_{mn.compact(c)}",
                                           acc, species, [c], phenotype, classes, []))
                n["parsed"] += 1
                continue
            before, errors = len(records), []
            for c in singles:
                for group in alternatives:
                    try:
                        others = []
                        for g, ref, pos, alts in group:
                            gacc = by_gene.get(g) or by_gene.get(pointfinder_gene(g)[0])
                            if gacc is None:
                                raise mn.NotationError(f"required mutation in unknown gene {g}")
                            # one combination per listed residue
                            type_gene = type_at.get((gacc, str(pos)), "AA")
                            others.append([parse(gacc, type_gene, str(pos), ref, a, "")[0] for a in alts])
                    except mn.NotationError as e:
                        errors.append(str(e))
                        continue
                    for combo in _product(others):
                        components = [c] + combo
                        notation = "+".join(f"{x['gene_label']}_{mn.compact(x)}" for x in components)
                        records.append(_record("pointfinder", record_id("pointfinder", components), notation, acc,
                                               species, components, phenotype, classes, [f"required: {required}"]))
            if len(records) > before:
                n["parsed"] += 1
                for e in sorted(set(errors)):   # some alternatives of the row were kept, others not
                    excluded.append({"database": "pointfinder", "row": text + f" (required {required})",
                                     "reason": e, "partial": True})
            else:
                excluded.append({"database": "pointfinder", "row": text + f" (required {required})",
                                 "reason": "; ".join(sorted(set(errors))) or "no combination could be built"})
    return {"records": merge_duplicates(records), "references": references, "excluded": excluded,
            "rows_read": n["read"], "rows_parsed": n["parsed"]}


def _product(lists):
    out = [[]]
    for options in lists:
        out = [prefix + [o] for prefix in out for o in options]
    return out
