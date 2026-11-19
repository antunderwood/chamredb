"""Look up mutations in the bundled comparison of CARD, NCBI and PointFinder (data/variants.json, pinned by sha256 in
data/graph_release.json). Positions are in the numbering of each locus's anchor (CARD, else NCBI, else PointFinder).
"""
import hashlib
import json
import re
from pathlib import Path

from chamredb import graph_loader

DATABASES = ("card", "ncbi", "pointfinder")
NAMES = {"ncbi": "NCBI", "card": "CARD", "pointfinder": "PointFinder"}
CLASS_TEXT = {
    "listed": "listed",
    "same_mutation": "same mutation (same reference sequence)",
    "same_mutation_different_reference": "same mutation, different reference sequence",
    "possible_same_mutation_offset": "possibly the same mutation, numbered differently (review)",
    "same_mutation_only_in_combination": "listed only in combination with other mutations",
    "same_position_different_change": "same position, different change",
    "different_mutation": "gene covered, mutation not listed",
    "same_mutation_related_organism": "same mutation listed in a related organism",
    "locus_not_covered": "no mutations listed in this gene for this organism",
    "organism_not_covered": "no mutations listed for this organism",
    "same_combination": "same combination",
    "partial_combination": "combination sharing some mutations",
    "components_listed_separately": "each mutation listed, but not as this combination",
    "combination_not_listed": "combination not listed",
}
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "TER": "*", "STOP": "*", "XAA": "?", "VAR": "?"}
CUSTOM_RELEASE = graph_loader.CUSTOM_RELEASE


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_variants(path=None, pin_file=None):
    """(file, data release). A given path is used as is; the bundled file is checked against its pin."""
    if path:
        return Path(path), CUSTOM_RELEASE
    pin_path = Path(pin_file or graph_loader.PACKAGE_DATA / "graph_release.json")
    pin = json.loads(pin_path.read_text())
    entry = pin.get("variants")
    if entry is None:
        raise ValueError(f"{pin_path} pins no variants file")
    file = pin_path.parent / entry["file"]
    if _sha256(file) != entry["sha256"]:
        raise ValueError(f"{file} sha256 does not match graph_release.json")
    return file, entry.get("data_release") or pin["data_release"]


# Three-letter residues are read in the case HGVS writes them (Ser, Ter), wherever they are (p.Arg177_Thr178insIleGly),
# so runs of one-letter residues that spell one (TERY, ASNG) are left alone. Any case is read in a plain substitution.
_THREE = re.compile("|".join(k.capitalize() for k in sorted(AA3, key=len, reverse=True)))
_THREE_SUBSTITUTION = re.compile(r"^([A-Za-z]{3})(-?\d+)([A-Za-z]{3,4}|[*=?])$")
_INDEL = re.compile(r"^[A-Z*]*(-?\d+)(?:_[A-Z*]*(-?\d+))?(del|dup|fs|ins|delins)([A-Z*?]*)$")
_POSITION = r"[A-Za-z*]{0,3}-?\d+"
_RANGE = r"(?:[pncg]\.)?" + _POSITION + "_" + _POSITION + r"(?:delins|del|dup|ins)[A-Za-z*?]*"
_GENE_AND_RANGE = re.compile(r"^(?P<gene>.+?)[ :_](?P<mut>" + _RANGE + r")$")


def normalise_mutation(text):
    """Short form of a mutation: p.Ser83Leu, p.S83L and S83L give S83L; Q10Ter gives Q10*; g2032t, n.2032G>T give
    G2032T. Indels are written without the residues beside the positions, which the databases give or leave out:
    p.Asn437_Asn438del and N437_N438del give 437_438del, p.Arg177_Thr178insIle gives 177_178insI, G129fs stays."""
    t = text.strip().replace(" ", "")
    t = re.sub(r"^[pncg]\.", "", t)
    m = _THREE_SUBSTITUTION.match(t)
    if m and m.group(1).upper() in AA3 and (m.group(3) in "*=?" or m.group(3).upper() in AA3):
        t = AA3[m.group(1).upper()] + m.group(2) + AA3.get(m.group(3).upper(), m.group(3))
    t = _THREE.sub(lambda m: AA3[m.group(0).upper()], t)
    m = re.match(r"^([A-Z])(-?\d+)=$", t)                               # p.Leu197= : unchanged residue
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(1)}"
    m = re.match(r"^(-?\d+)=$", t)
    if m:
        return t
    m = _INDEL.match(t)
    if m and m.group(3) != "fs":
        start, end, kind, new = m.groups()
        if kind in ("del", "dup") and new:
            return t
        return (start if end in (None, start) else f"{start}_{end}") + kind + new
    m = re.match(r"^(-?\d+)([ACGTUacgtu])>([ACGTUacgtuN])$", t)       # 2059A>G, -15C>T
    if m:
        alt = "?" if m.group(3).upper() == "N" else m.group(3).upper().replace("U", "T")
        return f"{m.group(2).upper().replace('U', 'T')}{m.group(1)}{alt}"
    m = re.match(r"^([acgtu])(-?\d+)([acgtu])$", t)                    # CARD rRNA g2032t
    if m:
        return f"{m.group(1).upper()}{m.group(2)}{m.group(3).upper()}".replace("U", "T")
    m = re.match(r"^(-?\d+)([acgt])>([acgt])$", t)
    if m:
        return f"{m.group(2).upper()}{m.group(1)}{m.group(3).upper()}"
    return t


def split_query(query):
    """'gyrA_S83L', 'gyrA S83L', 'gyrA:p.Ser83Leu' -> ('gyrA', 'S83L'); 'S83L' -> (None, 'S83L'). A mutation over a
    range keeps both ends: 'rpoB_N437_N438del' -> ('rpoB', '437_438del'), 'D63_S67del' -> (None, '63_67del')."""
    q = query.strip()
    if re.fullmatch(_RANGE, q):
        return None, normalise_mutation(q)
    m = _GENE_AND_RANGE.match(q)
    if m:
        return m.group("gene").strip(), normalise_mutation(m.group("mut"))
    for sep in (" ", ":", "_"):
        if sep in q:
            gene, mut = q.rsplit(sep, 1)
            if re.search(r"\d", mut):
                return gene.strip(), normalise_mutation(mut)
    return None, normalise_mutation(q)


def _label(text):
    return text.lower().replace("_", " ").strip()


class Variants:
    def __init__(self, data, release):
        self.data = data
        self.release = release
        self.records = data["records"]
        self.by_id = {r["record_id"]: r for r in self.records}
        self.loci = data["loci"]
        self.matches = {}
        for m in data["matches"]:
            self.matches.setdefault(m["record_id"], {})[m["other_database"]] = m
        self._form_cache = {}
        self.cluster_of = {}
        for c in data["clusters"]:
            for rid in c["members"]:
                self.cluster_of[rid] = c

    @classmethod
    def load(cls, path=None):
        file, release = resolve_variants(path)
        return cls(json.loads(Path(file).read_text()), release)

    def _gene_names(self, c):
        locus = self.loci[c["locus_id"]]
        names = {_label(c["gene_label"]), _label(locus["gene"])}
        if locus["gene"].endswith(" rRNA"):
            names.add(_label(locus["gene"][:-5]))
        return names

    def _forms(self, r, c):
        """Short forms a component can be asked for: its own and its anchor notation, compact and HGVS."""
        key = (r["record_id"], id(c))
        if key not in self._form_cache:
            forms = {c["native_compact"], c.get("anchor_compact"), c.get("native_hgvs"), c.get("hgvs")}
            if len(r["components"]) == 1:
                forms.add(split_query(r["native_notation"])[1])
            self._form_cache[key] = {normalise_mutation(f) for f in forms if f} | {f for f in forms if f}
        return self._form_cache[key]

    def _taxon_ok(self, r, taxon):
        if not taxon:
            return True
        t = _label(taxon)
        return t in _label(r["taxon_group"]) or t in _label(r["organism_label"]) or _label(r["taxon_group"]) in t

    def find(self, query, database=None, taxon=None, accession=None):
        """Records matching a record id, a database's own notation, or gene and mutation in a common form."""
        if query in self.by_id:
            return [self.by_id[query]]
        gene, mut = split_query(query)
        found = []
        for r in self.records:
            if database and r["database"] != database:
                continue
            if not self._taxon_ok(r, taxon):
                continue
            if accession and accession not in r["native_label"] and not any(
                    accession in (self.loci[c["locus_id"]]["anchor_accession"] or "") for c in r["components"]):
                continue
            if r["native_notation"] == query.strip():
                found.append(r)
                continue
            if len(r["components"]) != 1:
                continue
            c = r["components"][0]
            if mut in self._forms(r, c) and (gene is None or _label(gene) in self._gene_names(c)):
                found.append(r)
        return sorted(found, key=lambda r: (DATABASES.index(r["database"]), r["taxon_group"], r["record_id"]))

    def find_combination(self, mutations, database=None, taxon=None):
        """Combination records whose components are exactly the given (gene, mutation) pairs."""
        wanted = {(_label(g) if g else None, m) for g, m in mutations}
        out = []
        for r in self.records:
            if len(r["components"]) < 2 or (database and r["database"] != database) or not self._taxon_ok(r, taxon):
                continue
            have = set()
            for c in r["components"]:
                for g, m in wanted:
                    if m in self._forms(r, c) and (g is None or g in self._gene_names(c)):
                        have.add((g, m))
            if have == wanted and len(r["components"]) == len(wanted):
                out.append(r)
        return out

    def locus_ids(self, query):
        """Locus ids matching e.g. Escherichia:gyrA or 'Neisseria gonorrhoeae:23S'."""
        q = [_label(x) for x in query.split(":")]
        out = []
        for lid, locus in self.loci.items():
            group = _label(locus["group"])
            gene = locus["gene"].lower()
            if len(q) < 2 or q[0] not in group or q[1] not in (gene, gene.replace(" rrna", "")):
                continue
            if len(q) == 2 or q[2] == locus["molecule"]:
                out.append(lid)
        return sorted(out)

    def comparison(self, r):
        """Per other database: class, matched records, offset, note."""
        out = {}
        for db in DATABASES:
            if db == r["database"]:
                continue
            m = self.matches[r["record_id"]][db]
            out[db] = {"class": m["class"], "matched": [self.by_id[h] for h in m["matched"]],
                       "offset": m["offset"], "distance": m["distance"], "phenotype": m["phenotype"], "note": m["note"]}
        return out

    def describe(self, r):
        return " + ".join(c["hgvs"] or f"{c['native_hgvs']} (not placed on the anchor)" for c in r["components"])

    def anchor_of(self, r):
        loci = [self.loci[c["locus_id"]] for c in r["components"]]
        return "; ".join(f"{NAMES[l['anchor_database']]} {l['anchor_accession'] or '(no sequence)'}" for l in loci)

    def locus_clusters(self, locus_id):
        rows = [c for c in self.data["clusters"] if locus_id in c["loci"]]

        def pos(c):
            r = self.by_id[c["members"][0]]
            comp = next(x for x in r["components"] if x["locus_id"] == locus_id)
            return (len(c["loci"]) > 1 or c["kind"] == "combination",
                    comp["anchor_start"] if comp["anchor_start"] is not None else 10 ** 9, c["hgvs"])
        return sorted(rows, key=pos)

    def notations(self, cluster, database):
        return [self.by_id[m]["native_notation"] for m in cluster["members"] if self.by_id[m]["database"] == database]
