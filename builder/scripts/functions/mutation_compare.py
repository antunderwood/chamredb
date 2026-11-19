"""Compare the mutations of CARD, NCBI and PointFinder.

assign_groups and normalise_classes map organism labels and drug classes with the tables in builder/variants.
match_loci joins references of one organism and molecule that share a gene name or a near-identical sequence, and
picks the anchor of each locus (CARD, else NCBI, else PointFinder). place puts every change on its anchor, classify
compares every record with the other two databases, and cluster joins the records that are the same mutation.
"""
import csv
import hashlib
from collections import Counter, defaultdict

from functions import liftover
from functions import mutation_notation as mn

DATABASES = ("card", "ncbi", "pointfinder")   # also the anchor preference
WINDOW = 10                                    # positions; possible_same_mutation_offset
LOCUS_SEQUENCE_IDENTITY = 0.9                  # join differently named references of one organism
LOCUS_SEQUENCE_COVERAGE = 0.9
KMER_FILTER = 0.3
LIFT_MIN_IDENTITY = 0.5                        # below this a reference is not lifted onto the anchor
LIFT_MIN_COVERAGE = 0.5                        # share of the shorter sequence that must be aligned to lift
RRNA_ANTISENSE_MARGIN = 3                      # more sites must fit the antisense strand to read positions along it
WILDCARD_KINDS = ("any_substitution",)         # 'any residue' changes match specific ones but do not join clusters

SINGLE_CLASSES = (
    "same_mutation", "same_mutation_different_reference", "possible_same_mutation_offset",
    "same_mutation_only_in_combination", "same_position_different_change", "different_mutation",
    "same_mutation_related_organism", "locus_not_covered", "organism_not_covered",
)
COMBINATION_CLASSES = (
    "same_combination", "partial_combination", "components_listed_separately", "combination_not_listed",
    "locus_not_covered", "organism_not_covered",
)
CLUSTER_CLASSES = {"same_mutation", "same_mutation_different_reference", "same_combination"}
# One order per record type: locus_not_covered and organism_not_covered are in both lists at different ranks.
PRIORITY = {c: i for i, c in enumerate(SINGLE_CLASSES)}
COMBINATION_PRIORITY = {c: i for i, c in enumerate(COMBINATION_CLASSES)}


def _liftable(aln):
    return aln is not None and aln["identity"] >= LIFT_MIN_IDENTITY and aln["coverage"] >= LIFT_MIN_COVERAGE


def read_table(path):
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def organism_key(record):
    label = record["organism_label"]
    return " ".join(label.split()[:2]) if record["database"] == "card" else label


def assign_groups(records, taxa_rows, amrfinder_organisms=None):
    """Set the organism group of each record from taxa.tsv. AMRFinderPlus groups must be in its taxgroup.tsv."""
    taxa = {(r["database"], r["label"]): r for r in taxa_rows}
    unknown = sorted({(r["database"], organism_key(r)) for r in records} - set(taxa))
    if unknown:
        raise ValueError(f"organism labels missing from taxa.tsv: {unknown}")
    if amrfinder_organisms is not None:
        bad = sorted({(t["database"], t["label"], t["group"]) for t in taxa_rows
                      if t["group_source"] == "amrfinderplus" and t["group"] not in amrfinder_organisms})
        if bad:
            raise ValueError(f"taxa.tsv groups that are not AMRFinderPlus organisms: {bad}")
    for r in records:
        t = taxa[(r["database"], organism_key(r))]
        r["taxon_group"], r["taxon_family"], r["taxon_source"] = t["group"], t["group_family"], t["group_source"]


def normalise_classes(records, class_rows):
    table = {(r["database"], r["term"]): r["class"] for r in class_rows}
    unknown = sorted({(r["database"], t) for r in records for t in r["drug_classes_native"]} - set(table))
    if unknown:
        raise ValueError(f"drug class terms missing from drug_classes.tsv: {unknown}")
    for r in records:
        r["drug_classes"] = sorted({table[(r["database"], t)] for t in r["drug_classes_native"]} - {""})


def gene_key(label, aliases):
    return aliases.get(label, label).lower()


def sequence_sha256(seq):
    return hashlib.sha256(seq.upper().replace("U", "T").encode()).hexdigest()


def _seq(references, ref):
    """Sequence of a reference, or None (CARD promoters have none)."""
    entry = references.get(ref)
    return entry["sequence"] if entry and entry["sequence"] else None


def _numbered_seq(references, ref):
    """The sequence of a reference on the strand its database numbers it along: the stored one, unless orient_rrna
    found the positions to count along the other strand."""
    seq = _seq(references, ref)
    entry = references.get(ref) or {}
    if seq and entry.get("stored_strand", 1) * entry.get("numbering_strand", 1) == -1:
        return liftover.revcomp(seq)
    return seq


def _anchor_seq(references, locus):
    """The anchor sequence of a locus; for rRNA always the sense strand, whichever strand the reference stores."""
    seq = _seq(references, locus["anchor_reference"])
    if seq and references[locus["anchor_reference"]].get("stored_strand", 1) == -1:
        return liftover.revcomp(seq)
    return seq


def _kind(molecule):
    return "protein" if molecule == "protein" else "nucleotide"


class _UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def match_loci(records, references, aliases):
    """Return {locus_id: locus} and set component['locus_id'] on every component."""
    use = defaultdict(int)          # (group, molecule, reference) -> number of changes
    for r in records:
        for c in r["components"]:
            use[(r["taxon_group"], c["molecule"], c["reference"])] += 1
    uf = _UnionFind()
    by_name = {}
    for node in sorted(use):
        group, molecule, ref = node
        name = (group, molecule, gene_key(references_gene(references, ref), aliases))
        uf.find(node)
        if name in by_name:
            uf.union(node, by_name[name])
        else:
            by_name[name] = node
    joined_by_sequence = []
    by_gm = defaultdict(list)
    for node in sorted(use):
        if node[1] != "promoter" and _seq(references, node[2]):
            by_gm[node[:2]].append(node)
    for (group, molecule), nodes in sorted(by_gm.items()):
        for i, a in enumerate(nodes):
            for b in nodes[i + 1:]:
                if uf.find(a) == uf.find(b) or references[a[2]]["database"] == references[b[2]]["database"]:
                    continue
                sa, sb = references[a[2]]["sequence"], references[b[2]]["sequence"]
                if liftover.kmer_similarity(sa, sb) < KMER_FILTER:
                    continue
                aln = liftover.align(sa, sb, _kind(molecule))
                if aln["identity"] >= LOCUS_SEQUENCE_IDENTITY and aln["coverage"] >= LOCUS_SEQUENCE_COVERAGE:
                    uf.union(a, b)
                    joined_by_sequence.append({"group": group, "molecule": molecule, "a": a[2], "b": b[2],
                                               "identity": round(aln["identity"], 4), "coverage": round(aln["coverage"], 4)})
    members = defaultdict(list)
    for node in sorted(use):
        members[uf.find(node)].append(node)

    loci, locus_of = {}, {}
    for root, nodes in sorted(members.items()):
        group, molecule = root[0], root[1]

        def rank(node):
            ref = node[2]
            db = references[ref]["database"]
            has_seq = _seq(references, ref) is not None
            return (not has_seq, DATABASES.index(db), -use[node], ref)
        anchor = sorted(nodes, key=rank)[0][2]
        anchor_ref = references[anchor] if _seq(references, anchor) else None
        gene = references_gene(references, anchor)
        display = aliases.get(gene, gene)
        locus_id = f"{group}:{display}:{molecule}"
        if locus_id in loci:   # two loci with one display name (e.g. two unrelated references named alike)
            n = 2
            while f"{locus_id}#{n}" in loci:
                n += 1
            locus_id = f"{locus_id}#{n}"
        refs = []
        for node in nodes:
            ref = node[2]
            entry = {"reference": ref, "database": ref.split("|")[0], "changes": use[node],
                     "gene_label": references_gene(references, ref)}
            if _seq(references, ref):
                entry["accession"] = references[ref]["accession"]
                entry["sha256"] = sequence_sha256(references[ref]["sequence"])
                entry["length"] = len(references[ref]["sequence"])
            refs.append(entry)
            locus_of[node] = locus_id
        loci[locus_id] = {
            "locus_id": locus_id, "group": group, "gene": display, "molecule": molecule,
            "anchor_reference": anchor, "anchor_database": anchor.split("|")[0],
            "anchor_accession": anchor_ref["accession"] if anchor_ref else None,
            "anchor_sha256": sequence_sha256(anchor_ref["sequence"]) if anchor_ref else None,
            "references": refs,
        }
    for r in records:
        for c in r["components"]:
            c["locus_id"] = locus_of[(r["taxon_group"], c["molecule"], c["reference"])]
    return loci, joined_by_sequence


def references_gene(references, ref):
    return references[ref]["gene_label"]


STRAND_NAMES = {1: "sense", -1: "antisense", 0: "unknown"}


def orient_rrna(records, loci, references):
    """Find the strand every rRNA reference is stored on and the strand its database numbers it along.

    Databases store some rRNA genes as the reverse complement of the rRNA, and number them either along the
    stored sequence or along the rRNA. The stored strand comes from conserved rRNA stretches. Positions are read
    along the sense strand, as the literature numbers rRNA, unless the stated original residues fit the antisense
    strand at RRNA_ANTISENSE_MARGIN more sites than the sense strand. Sets stored_strand and numbering_strand on
    the references and the locus entries, and returns one row per reference for the report.
    """
    sites = defaultdict(set)
    for r in records:
        for c in r["components"]:
            if c["molecule"] == "rna" and c["ref"] and c["index"] is not None:
                sites[c["reference"]].add((c["index"], c["ref"].upper().replace("U", "T")))
    rows = []
    for locus in loci.values():
        if locus["molecule"] != "rna":
            continue
        for entry in locus["references"]:
            ref = entry["reference"]
            seq = _seq(references, ref)
            if not seq:
                continue
            stored = liftover.rrna_strand(seq)
            sense = liftover.revcomp(seq) if stored == -1 else seq
            antisense = liftover.revcomp(sense)
            fit = {name: sum(_residues_at(strand, i, len(x)) == x for i, x in sites[ref])
                   for name, strand in (("sense", sense), ("antisense", antisense))}
            if stored == 0:      # positions are read on the sequence as stored
                numbering = 0
            else:
                numbering = -1 if fit["antisense"] - fit["sense"] >= RRNA_ANTISENSE_MARGIN else 1
            references[ref]["stored_strand"], references[ref]["numbering_strand"] = stored, numbering
            entry["stored_strand"], entry["numbering_strand"] = STRAND_NAMES[stored], STRAND_NAMES[numbering]
            rows.append({"locus_id": locus["locus_id"], "database": entry["database"], "reference": ref,
                         "stored_strand": STRAND_NAMES[stored], "numbering_strand": STRAND_NAMES[numbering],
                         "sites": len(sites[ref]), "sites_fit_sense": fit["sense"] if stored else None,
                         "sites_fit_antisense": fit["antisense"] if stored else None})
    return sorted(rows, key=lambda x: (x["locus_id"], x["database"], x["reference"]))


def lift_maps(loci, references):
    """{reference: alignment to its locus anchor} for every reference with a sequence (not promoters)."""
    maps = {}
    for locus in loci.values():
        if locus["molecule"] == "promoter" or not _seq(references, locus["anchor_reference"]):
            continue
        anchor_seq = _anchor_seq(references, locus)
        for entry in locus["references"]:
            ref = entry["reference"]
            if not _seq(references, ref):
                continue
            aln = liftover.align(anchor_seq, _numbered_seq(references, ref), _kind(locus["molecule"]))
            entry["identity_to_anchor"] = round(aln["identity"], 4)
            entry["coverage_to_anchor"] = round(aln["coverage"], 4)
            entry["strand_to_anchor"] = aln["strand"]
            entry["lifted"] = _liftable(aln)
            maps[ref] = aln
    return maps


def _residues_at(seq, index, n):
    if index is None or index < 1 or index - 1 + n > len(seq):
        return None
    return seq[index - 1:index - 1 + n].upper().replace("U", "T")


def self_check(c, references):
    """The source's stated reference residues against its own sequence (on the strand its positions are read
    along): ok, mismatch or unchecked."""
    if not c["ref"] or c["index"] is None or not _seq(references, c["reference"]):
        return "unchecked"
    got = _residues_at(_numbered_seq(references, c["reference"]), c["index"], len(c["ref"]))
    return "ok" if got == c["ref"].upper().replace("U", "T") else "mismatch"


def on_anchor(c):
    """The change as written on its anchor: anchor positions, and residues complemented if the reference is stored
    on the other strand. The residue named beside an insertion is on its other side then, so it is dropped."""
    out = dict(c, start=c["anchor_start"], end=c["anchor_end"])
    if c.get("anchor_strand", 1) == -1:
        out["ref"] = "" if c["kind"] == "insertion" else liftover.revcomp(c["ref"])
        out["alt"] = liftover.revcomp(c["alt"])
    return out


def place(records, loci, references, maps):
    for r in records:
        for c in r["components"]:
            locus = loci[c["locus_id"]]
            seq = _seq(references, c["reference"])
            c["reference_sha256"] = sequence_sha256(seq) if seq else None
            c["ref_check"] = self_check(c, references)
            c["native_hgvs"] = mn.hgvs(c)
            c["anchor_strand"] = 1
            if c["molecule"] == "promoter":
                c["anchor_start"], c["anchor_end"] = c["start"], c["end"]
                c["anchor_ref_check"] = "unchecked"
            else:
                aln = maps.get(c["reference"])
                if not _liftable(aln):
                    c["anchor_start"], c["anchor_end"] = None, None
                else:
                    c["anchor_start"], c["anchor_end"] = liftover.lift(aln["map"], c["start"], c["end"])
                    c["anchor_strand"] = aln["strand"]
                anchor_seq = _anchor_seq(references, locus)
                ref = on_anchor(c)["ref"] if c["anchor_start"] is not None else ""
                if not ref or anchor_seq is None:
                    c["anchor_ref_check"] = "unchecked"
                else:
                    got = _residues_at(anchor_seq, c["anchor_start"], 1)
                    c["anchor_ref_check"] = "ok" if got == ref[0].upper().replace("U", "T") else "mismatch"
            c["hgvs"] = mn.hgvs(on_anchor(c)) if c["anchor_start"] is not None else None
            c["native_compact"] = mn.compact(c)
            c["anchor_compact"] = mn.compact(on_anchor(c)) if c["anchor_start"] is not None else None
        r["locus_ids"] = sorted({c["locus_id"] for c in r["components"]})
        r["kind"] = "combination" if len(r["components"]) > 1 else r["components"][0]["kind"]


def canonical(c, references, loci):
    """(start, end, key) on the anchor; a duplication is written as an insertion after the duplicated residues."""
    if c["anchor_start"] is None:
        return None
    s, e = c["anchor_start"], c["anchor_end"]
    if c["kind"] == "duplication":
        anchor = _anchor_seq(references, loci[c["locus_id"]])
        if anchor is None:
            return (s, e, ("duplication", e - s))
        return (e, e + 1, ("insertion", anchor[s - 1:e].upper()))
    return (s, e, mn.change_key(on_anchor(c)))


def _keys_match(k1, k2):
    sub = ("substitution", "nonsense", "any_substitution")
    if k1[0] == "any_substitution" and k2[0] in sub or k2[0] == "any_substitution" and k1[0] in sub:
        return True
    return k1 == k2


def compare_changes(c1, c2, references, loci):
    """Class of change c2 relative to c1 (same locus), with the anchor distance where defined."""
    p1, p2 = canonical(c1, references, loci), canonical(c2, references, loci)
    same_native = (c1["reference_sha256"] is not None and c1["reference_sha256"] == c2["reference_sha256"]
                   and c1["start"] == c2["start"] and c1["end"] == c2["end"])
    k1 = p1[2] if p1 else mn.change_key(c1)
    k2 = p2[2] if p2 else mn.change_key(c2)
    if same_native and _keys_match(k1, k2):
        return "same_mutation", 0
    same_change = _keys_match(k1, k2) and _refs_agree(c1, c2)
    # the two databases write the same thing (number and residues) although their sequences do not put it at
    # one place: compared as written, since references stored on opposite strands complement the residues
    same_notation = (c1["start"] == c2["start"] and _keys_match(mn.change_key(c1), mn.change_key(c2))
                     and _refs_agree(dict(c1, anchor_start=None), c2))
    if p1 and p2:
        delta = p2[0] - p1[0]
        if p1[:2] == p2[:2] and _keys_match(k1, k2):
            return "same_mutation_different_reference", 0
        if same_notation or (same_change and abs(delta) <= WINDOW and _numbering_suspect(c1, c2)):
            return "possible_same_mutation_offset", delta
        if p1[0] == p2[0]:
            return "same_position_different_change", 0
        return "different_mutation", delta
    if same_notation:
        return "possible_same_mutation_offset", None
    return "different_mutation", None


def _refs_agree(c1, c2):
    """The first original residues agree, on the anchor strand where both changes are placed."""
    if c1["anchor_start"] is not None and c2["anchor_start"] is not None:
        c1, c2 = on_anchor(c1), on_anchor(c2)
    a, b = c1["ref"][:1].upper().replace("U", "T"), c2["ref"][:1].upper().replace("U", "T")
    return not a or not b or a == b


def _numbering_suspect(c1, c2):
    """True if a reference check fails for either change, or for promoter changes."""
    if c1["molecule"] == "promoter":
        return True
    return "mismatch" in (c1["ref_check"], c2["ref_check"], c1["anchor_ref_check"], c2["anchor_ref_check"])


def _phenotype(r1, r2):
    a, b = set(r1["drug_classes"]), set(r2["drug_classes"])
    if not a or not b:
        return "unknown"
    if a == b:
        return "same_classes"
    return "shared_class" if a & b else "no_shared_class"


class Index:
    def __init__(self, records, loci, references):
        self.records = records
        self.loci = loci
        self.references = references
        self.singles = defaultdict(list)       # (db, locus) -> records
        self.combos = defaultdict(list)        # (db, group) -> records
        self.groups = defaultdict(int)         # (db, group) -> records
        self.loci_by_gene = defaultdict(list)  # (family, gene key, molecule) -> loci
        for r in records:
            self.groups[(r["database"], r["taxon_group"])] += 1
            if len(r["components"]) == 1:
                self.singles[(r["database"], r["components"][0]["locus_id"])].append(r)
            else:
                self.combos[(r["database"], r["taxon_group"])].append(r)
        self.covered = defaultdict(set)        # db -> loci with any change of db
        for r in records:
            for c in r["components"]:
                self.covered[r["database"]].add(c["locus_id"])
        self.family_of = {r["taxon_group"]: r["taxon_family"] for r in records}
        for locus in loci.values():
            fam = self.family_of.get(locus["group"], "")
            if fam:
                self.loci_by_gene[(fam, locus["gene"].lower(), locus["molecule"])].append(locus["locus_id"])
        self._cross = {}

    def cross_map(self, from_locus, to_locus):
        """Alignment of the anchor of from_locus onto the anchor of to_locus (related organisms)."""
        key = (from_locus, to_locus)
        if key not in self._cross:
            a = _anchor_seq(self.references, self.loci[to_locus])
            b = _anchor_seq(self.references, self.loci[from_locus])
            if a is None or b is None or self.loci[to_locus]["molecule"] == "promoter":
                self._cross[key] = None
            else:
                aln = liftover.align(a, b, _kind(self.loci[to_locus]["molecule"]))
                self._cross[key] = aln if _liftable(aln) else None
        return self._cross[key]


def classify_single(r, other_db, ix):
    c = r["components"][0]
    candidates = ix.singles.get((other_db, c["locus_id"]), [])
    found = _classify_against_singles(c, candidates, ix) if candidates else None
    if found is None or PRIORITY[found[0]] > PRIORITY["same_mutation_only_in_combination"]:
        in_combos = [s for s in ix.combos.get((other_db, r["taxon_group"]), []) if _component_found(c, [s], ix)]
        if in_combos:
            return "same_mutation_only_in_combination", in_combos, None, "listed only together with other mutations"
    if found is not None:
        return found
    if c["locus_id"] in ix.covered[other_db]:   # covered only through combinations of other_db
        return "different_mutation", [], None, "locus covered by combinations only"
    related = _related(c, r, other_db, ix)
    if related:
        return related
    if ix.groups.get((other_db, r["taxon_group"])):
        return "locus_not_covered", [], None, ""
    return "organism_not_covered", [], None, ""


def _classify_against_singles(c, candidates, ix):
    if candidates:
        best_cls, best, best_delta = None, [], None
        for s in candidates:
            cls, delta = compare_changes(c, s["components"][0], ix.references, ix.loci)
            if best_cls is None or PRIORITY[cls] < PRIORITY[best_cls]:
                best_cls, best, best_delta = cls, [s], delta
            elif cls == best_cls:
                if cls == "different_mutation":
                    if delta is not None and (best_delta is None or abs(delta) < abs(best_delta)):
                        best, best_delta = [s], delta
                else:
                    best.append(s)
        if best_cls in ("same_mutation", "same_mutation_different_reference"):
            # all records that are the same change, so the relation is symmetric
            best = [s for s in candidates if compare_changes(c, s["components"][0], ix.references, ix.loci)[0]
                    in ("same_mutation", "same_mutation_different_reference")]
        return best_cls, best, best_delta, ""
    return None


def _related(c, r, other_db, ix):
    if not r["taxon_family"] or c["anchor_start"] is None:
        return None
    locus = ix.loci[c["locus_id"]]
    for other_locus in ix.loci_by_gene.get((r["taxon_family"], locus["gene"].lower(), locus["molecule"]), []):
        if other_locus == c["locus_id"] or not ix.singles.get((other_db, other_locus)):
            continue
        aln = ix.cross_map(c["locus_id"], other_locus)
        if aln is None:
            continue
        s_, e_ = liftover.lift(aln["map"], c["anchor_start"], c["anchor_end"])
        if s_ is None:
            continue
        key = canonical(c, ix.references, ix.loci)[2]
        if aln["strand"] == -1:   # the two anchors are stored on opposite strands
            if key[0] == "duplication":
                continue
            key = key[:1] + tuple(liftover.revcomp(x) if isinstance(x, str) else x for x in key[1:])
        hits = []
        for s in ix.singles[(other_db, other_locus)]:
            p = canonical(s["components"][0], ix.references, ix.loci)
            if p and p[0] == s_ and p[1] == e_ and _keys_match(key, p[2]):
                hits.append(s)
        if hits:
            return "same_mutation_related_organism", hits, None, f"in {ix.loci[other_locus]['group']}"
    return None


def _component_found(c, records, ix):
    for s in records:
        for sc in s["components"]:
            if sc["locus_id"] == c["locus_id"]:
                cls, _ = compare_changes(c, sc, ix.references, ix.loci)
                if cls in ("same_mutation", "same_mutation_different_reference"):
                    return True
    return False


def classify_combination(r, other_db, ix):
    comps = r["components"]
    best_cls, best, shared = None, [], 0
    for s in ix.combos.get((other_db, r["taxon_group"]), []):
        n = sum(_component_found(c, [s], ix) for c in comps)
        if n == 0:
            continue
        cls = "same_combination" if n == len(comps) == len(s["components"]) else "partial_combination"
        if best_cls is None or COMBINATION_PRIORITY[cls] < COMBINATION_PRIORITY[best_cls] or (cls == best_cls and n > shared):
            best_cls, best, shared = cls, [s], n
        elif cls == best_cls and n == shared:
            best.append(s)
    if best_cls:
        return best_cls, best, None, f"{shared} of {len(comps)} components shared" if best_cls == "partial_combination" else ""
    singles = [s for c in comps for s in ix.singles.get((other_db, c["locus_id"]), [])]
    if all(_component_found(c, singles, ix) for c in comps):
        return "components_listed_separately", [], None, ""
    if any(c["locus_id"] in ix.covered[other_db] for c in comps):
        return "combination_not_listed", [], None, ""
    if ix.groups.get((other_db, r["taxon_group"])):
        return "locus_not_covered", [], None, ""
    return "organism_not_covered", [], None, ""


def classify(records, loci, references, reviewed=()):
    """One match row per (record, other database)."""
    ix = Index(records, loci, references)
    decisions = {}
    for row in reviewed:
        decisions[(row["record_a"], row["record_b"])] = row
        decisions[(row["record_b"], row["record_a"])] = row
    matches = []
    for r in records:
        for other in DATABASES:
            if other == r["database"]:
                continue
            if len(r["components"]) == 1:
                cls, hits, delta, note = classify_single(r, other, ix)
            else:
                cls, hits, delta, note = classify_combination(r, other, ix)
            if cls == "possible_same_mutation_offset":
                for h in hits:
                    d = decisions.get((r["record_id"], h["record_id"]))
                    if d and d["decision"] == "same":
                        cls, note = "same_mutation_different_reference", f"reviewed equivalence ({d['reviewer']}, {d['date']})"
                        hits = [h]
                        break
            matches.append({
                "record_id": r["record_id"], "database": r["database"], "other_database": other, "class": cls,
                "matched": [h["record_id"] for h in hits], "offset": delta if cls == "possible_same_mutation_offset" else None,
                "distance": abs(delta) if cls == "different_mutation" and delta is not None else None,
                "phenotype": _phenotype(r, hits[0]) if hits else None, "note": note,
            })
    return matches


def cluster(records, matches):
    by_id = {r["record_id"]: r for r in records}

    def wildcards(rid):
        return sum(1 for c in by_id[rid]["components"] if c["kind"] in WILDCARD_KINDS)

    uf = _UnionFind()
    for r in records:
        uf.find(r["record_id"])
    for m in matches:
        if m["class"] in CLUSTER_CLASSES:
            for h in m["matched"]:
                # S315Var matches S315T and S315N, but those two are not one mutation: an 'any residue' record
                # joins only records that are 'any residue' too
                if wildcards(m["record_id"]) == wildcards(h):
                    uf.union(m["record_id"], h)
    groups = defaultdict(list)
    for r in records:
        groups[uf.find(r["record_id"])].append(r)
    by_record = defaultdict(list)
    for m in matches:
        by_record[m["record_id"]].append(m)
    clusters = []
    for members in groups.values():
        members.sort(key=lambda r: (DATABASES.index(r["database"]), r["record_id"]))
        lead = members[0]
        status = {}
        for db in DATABASES:
            mine = [r for r in members if r["database"] == db]
            if mine:
                status[db] = "listed"
                continue
            classes = [m["class"] for r in members for m in by_record[r["record_id"]] if m["other_database"] == db]
            order = COMBINATION_PRIORITY if len(lead["components"]) > 1 else PRIORITY
            status[db] = min(classes, key=lambda c: order[c]) if classes else "organism_not_covered"
        clusters.append({
            "cluster_id": None, "group": lead["taxon_group"], "loci": lead["locus_ids"], "kind": lead["kind"],
            "hgvs": " + ".join(c["hgvs"] or f"[{c['native_hgvs']} not placed]" for c in lead["components"]),
            "members": [r["record_id"] for r in members], "status": status,
            "databases": sum(1 for v in status.values() if v == "listed"),
        })
    clusters.sort(key=lambda c: (c["group"], c["loci"], c["hgvs"], c["members"]))
    for i, c in enumerate(clusters, 1):
        c["cluster_id"] = f"V{i:05d}"
    return clusters


def systematic_offsets(matches, records_by_id, min_count=3):
    """The same anchor offset between two databases at one locus, seen at least min_count times."""
    seen = Counter()
    for m in matches:
        if m["class"] == "possible_same_mutation_offset" and m["offset"] is not None:
            r = records_by_id[m["record_id"]]
            seen[(r["components"][0]["locus_id"], m["database"], m["other_database"], m["offset"])] += 1
    return [{"locus_id": k[0], "database": k[1], "other_database": k[2], "offset": k[3], "count": n}
            for k, n in sorted(seen.items()) if n >= min_count]


def self_check_offsets(records, references, min_count=3, max_shift=150):
    """References whose self-check mismatches are explained by one constant shift (>= min_count changes)."""
    by_ref = defaultdict(list)
    for r in records:
        for c in r["components"]:
            if c["ref"] and c["index"] is not None and _seq(references, c["reference"]) and c["kind"] in (
                    "substitution", "nonsense", "any_substitution", "synonymous", "frameshift"):
                by_ref[c["reference"]].append(c)
    out = []
    for ref, cs in sorted(by_ref.items()):
        bad = [c for c in cs if c["ref_check"] == "mismatch"]
        if len(bad) < min_count:
            continue
        seq = _numbered_seq(references, ref).upper().replace("U", "T")

        def ok(c, d):
            i = c["index"] - 1 + d
            return 0 <= i < len(seq) and seq[i] == c["ref"].upper().replace("U", "T")
        best = max(range(-max_shift, max_shift + 1), key=lambda d: (sum(ok(c, d) for c in cs), -abs(d)))
        fixed = sum(ok(c, best) for c in bad)
        if best != 0 and fixed >= min_count:
            out.append({"reference": ref, "gene": references[ref]["gene_label"], "organism": references[ref]["organism_label"],
                        "changes": len(cs), "mismatches": len(bad), "shift": best, "mismatches_explained": fixed})
    return out
