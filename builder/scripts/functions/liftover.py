"""Global alignment of two reference sequences and the position map between them. End gaps are free, so a
shorter reference (a truncated gene) still maps onto the anchor. Nucleotide references are not all stored on the
same strand (rRNA genes in particular), so a nucleotide sequence is aligned as its reverse complement when that
shares more with the anchor; the result then has strand -1 and its changes must be complemented (revcomp).
"""
from functools import lru_cache

from Bio import Align
from Bio.Align import substitution_matrices


@lru_cache(maxsize=None)
def _aligner(kind):
    a = Align.PairwiseAligner()
    a.mode = "global"
    if kind == "protein":
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score, a.extend_gap_score = -10.0, -0.5
    else:
        a.match_score, a.mismatch_score = 2.0, -3.0
        a.open_gap_score, a.extend_gap_score = -5.0, -2.0
    if "end_insertion_score" in dir(a):   # Biopython >= 1.85 names
        a.end_insertion_score = a.end_deletion_score = 0.0
    else:                                 # Biopython 1.81 (builder environment)
        a.target_end_gap_score = a.query_end_gap_score = 0.0
    return a


def _clean(seq, kind):
    s = seq.upper().replace("*", "") if kind == "protein" else seq.upper().replace("U", "T")
    return "".join(ch for ch in s if ch.isalpha())


_COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def revcomp(seq):
    """Reverse complement; characters that are not nucleotides ('?' for any residue) are kept."""
    return seq.translate(_COMPLEMENT)[::-1]


# Stretches of 16S and 23S rRNA that are conserved across bacteria, archaea and organelles, as written on the
# sense strand (the rRNA itself). They tell which strand of the gene a reference sequence is stored on.
RRNA_SENSE_MOTIFS = (
    "ACGGGAGGCAGCAG", "GCCAGCAGCCGCGGTAA", "GCCAGCCGCCGCGGTAA", "CACACCGCCCGTC", "AGTCGTAACAAGG", "ATTAGATACCC",
    "CAACGAGCGCAACCC",                                                                             # 16S
    "CGGTGGATGCC", "GTACGAGAGGACC", "GATAACAGGCTGAT", "CGGAAAGACCCC", "CCTCGATGTCGG", "AAGGTTAAGGAAG",   # 23S
)


def rrna_strand(seq):
    """Strand an rRNA gene sequence is stored on: 1 sense (reads as the rRNA), -1 antisense (its reverse
    complement), 0 if the conserved stretches do not tell."""
    s = _clean(seq, "nucleotide")
    rc = revcomp(s)
    sense = sum(m in s for m in RRNA_SENSE_MOTIFS)
    antisense = sum(m in rc for m in RRNA_SENSE_MOTIFS)
    return (sense > antisense) - (antisense > sense)


def _shared_kmers(a, b, k=11):
    kb = {b[i:i + k] for i in range(len(b) - k + 1)}
    return sum(1 for i in range(len(a) - k + 1) if a[i:i + k] in kb)


def align(anchor, other, kind):
    """Align other onto anchor: identity, coverage of the shorter sequence, {other position: anchor position} and
    strand (-1 if other was aligned as its reverse complement; positions in the map are still those of other)."""
    if kind not in ("protein", "nucleotide"):
        raise ValueError(kind)
    a, b = _clean(anchor, kind), _clean(other, kind)
    if not a or not b:
        return {"identity": 0.0, "coverage": 0.0, "map": {}, "strand": 1}
    if a == b:
        return {"identity": 1.0, "coverage": 1.0, "map": {i: i for i in range(1, len(a) + 1)}, "strand": 1}
    strand = 1
    if kind == "nucleotide":
        rc = revcomp(b)
        if _shared_kmers(a, rc) > _shared_kmers(a, b):
            b, strand = rc, -1
    aln = _aligner(kind).align(a, b)[0]
    pos_map, same, cols = {}, 0, 0
    for (ta, tb), (qa, qb) in zip(*aln.aligned):
        ta, tb, qa = int(ta), int(tb), int(qa)   # numpy integers in recent Biopython
        for k in range(tb - ta):
            q = qa + k + 1
            pos_map[q if strand == 1 else len(b) - q + 1] = ta + k + 1
            cols += 1
            same += a[ta + k] == b[qa + k]
    return {"identity": same / cols if cols else 0.0, "coverage": cols / min(len(a), len(b)), "map": pos_map,
            "strand": strand}


def kmer_similarity(a, b, k=4):
    """Share of k-mers of the shorter sequence found in the other; a cheap filter before aligning."""
    a, b = a.upper(), b.upper()
    ka = {a[i:i + k] for i in range(len(a) - k + 1)}
    kb = {b[i:i + k] for i in range(len(b) - k + 1)}
    if not ka or not kb:
        return 0.0
    return len(ka & kb) / min(len(ka), len(kb))


def lift(pos_map, start, end):
    """Anchor (start, end) of a change at (start, end) of the other sequence, or (None, None) if either end falls
    in a gap of the alignment."""
    s, e = pos_map.get(start), pos_map.get(end)
    if s is None or e is None:
        return None, None
    return min(s, e), max(s, e)   # a reverse-strand map turns the range around
