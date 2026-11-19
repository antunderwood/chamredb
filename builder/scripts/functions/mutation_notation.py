"""Parse CARD, NCBI and PointFinder mutation notation into one form, and write it back as HGVS.

A change is a dict with molecule (protein, rna, dna, promoter), kind (substitution, any_substitution, nonsense,
frameshift, deletion, insertion, duplication, delins, synonymous), start, end, ref and alt. Positions are 1-based on
the source's own reference; promoter positions are negative offsets from the gene start, with no position 0.
An insertion has start = residue before it and end = start + 1. alt is '*' for a stop and '?' for any residue.
"""
import re

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His",
    "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "Ter", "?": "Xaa", "U": "Sec", "X": "Xaa",
}
AA1 = {v.upper(): k for k, v in AA3.items() if k not in "?X"}
AA = "ACDEFGHIKLMNPQRSTVWY"
NUC = "ACGTU"


class NotationError(ValueError):
    pass


def change(molecule, kind, start, end, ref="", alt=""):
    return {"molecule": molecule, "kind": kind, "start": int(start), "end": int(end), "ref": ref, "alt": alt}


def _aa(token):
    """One-letter amino acid from 'S', 'Ser' or 'Ter'."""
    if len(token) == 1 and (token in AA or token == "*"):
        return token
    t = token.upper()
    if t in AA1:
        return AA1[t]
    if t == "TER" or t == "STOP":
        return "*"
    raise NotationError(f"not an amino acid: {token!r}")


def substitution(molecule, pos, ref, alt):
    """Classify a one-position change. alt '*' is nonsense, alt '?' any substitution, ref == alt synonymous."""
    if alt == "*":
        kind = "nonsense"
    elif alt == "?":
        kind = "any_substitution"
    elif alt == ref:
        kind = "synonymous"
    else:
        kind = "substitution"
    return change(molecule, kind, pos, pos, ref, alt)


_CARD_AA_SUB = re.compile(r"^([A-Z])(\d+)([A-Z]|Ter|Var)$")
_CARD_NUC_SUB = re.compile(r"^([acgtu])(\d+)([acgtu])$")
_CARD_PROMOTER = re.compile(r"^(-\d+)([acgt])>([acgt]|Var)$")
_CARD_FS = re.compile(r"^([A-Z]|Ter)(\d+)fs$")
_CARD_DEL = re.compile(r"^([A-Z])(\d+)(?:_([A-Z])(\d+))?del$")
_CARD_INS = re.compile(r"^([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$")
_CARD_DUP = re.compile(r"^([A-Z])(\d+)(?:_([A-Z])(\d+))?dup$")


def parse_card(token, molecule):
    """One CARD mutation, e.g. S83L, C52Ter, R447Var, G129fs, D63_S67del, R97_G98insA, P627dup, g2032t, -11a>g."""
    t = token.strip()
    m = _CARD_PROMOTER.match(t)
    if m:
        alt = "?" if m.group(3) == "Var" else m.group(3).upper()
        return substitution("promoter", int(m.group(1)), m.group(2).upper(), alt)
    m = _CARD_NUC_SUB.match(t)
    if m:
        mol = "rna" if molecule == "rna" else "dna"
        return substitution(mol, int(m.group(2)), m.group(1).upper(), m.group(3).upper())
    if molecule == "rna":
        raise NotationError(f"CARD rRNA mutation not understood: {token!r}")
    m = _CARD_AA_SUB.match(t)
    if m:
        alt = {"Ter": "*", "Var": "?"}.get(m.group(3), m.group(3))
        if m.group(1) not in AA or (alt not in AA and alt not in "*?"):
            raise NotationError(f"CARD mutation with a residue that is not an amino acid: {token!r}")
        return substitution("protein", int(m.group(2)), m.group(1), alt)
    m = _CARD_FS.match(t)
    if m:
        return change("protein", "frameshift", m.group(2), m.group(2), _aa(m.group(1)))
    m = _CARD_DEL.match(t)
    if m:
        start = int(m.group(2))
        end = int(m.group(4)) if m.group(4) else start
        ref = m.group(1) if not m.group(3) else ""
        if m.group(3) and end == start + 1:
            ref = m.group(1) + m.group(3)
        return change("protein", "deletion", start, end, ref)
    m = _CARD_INS.match(t)
    if m:
        if int(m.group(4)) != int(m.group(2)) + 1:
            raise NotationError(f"CARD insertion between non-adjacent residues: {token!r}")
        return change("protein", "insertion", m.group(2), m.group(4), m.group(1) + m.group(3), m.group(5))
    m = _CARD_DUP.match(t)
    if m:
        start = int(m.group(2))
        end = int(m.group(4)) if m.group(4) else start
        return change("protein", "duplication", start, end, m.group(1) if start == end else "")
    raise NotationError(f"CARD mutation not understood: {token!r}")


def parse_card_codependent(value):
    """'45612,N87K,V199A+45613,D481E' -> [('45612', ['N87K', 'V199A']), ('45613', ['D481E'])] (ARO ids)."""
    parts = []
    for block in value.split("+"):
        fields = [f.strip() for f in block.split(",") if f.strip()]
        if len(fields) < 2 or not fields[0].isdigit():
            raise NotationError(f"CARD co-dependent variant not understood: {value!r}")
        parts.append((fields[0], fields[1:]))
    return parts


_NCBI_FS = re.compile(r"^([A-Z])(\d+)([A-Z])?fsTer(\d+)(?:(ins|del)(\d+))?$")
_NCBI_SIMPLE = re.compile(r"^([A-Z]+)(-?\d+)([A-Z]+|Ter|del)$")


def parse_ncbi(symbol, position, molecule):
    """NCBI symbol without the gene (S83L, Q163Ter, ASLD153del, G585GTG, A2074G, T-35G) at mutation_position.
    Returns (change, index in the reference sequence)."""
    m = _NCBI_FS.match(symbol)
    if m:
        if molecule != "protein":
            raise NotationError(f"NCBI frameshift outside a protein: {symbol!r}")
        return change("protein", "frameshift", position, position, m.group(1)), position
    m = _NCBI_SIMPLE.match(symbol)
    if not m:
        raise NotationError(f"NCBI mutation not understood: {symbol!r}")
    ref, sympos, alt = m.group(1), int(m.group(2)), m.group(3)
    # the symbol can use a conventional numbering; mutation_position is the position in the sequence
    if molecule == "promoter" or sympos < 0:
        start = sympos
        molecule = "promoter"
    else:
        start = position
    if alt == "Ter":
        return change(molecule, "nonsense", start, start, ref[0], "*") if len(ref) == 1 else \
            change(molecule, "delins", start, start + len(ref) - 1, ref, "*"), position
    if alt == "del":
        return change(molecule, "deletion", start, _advance(start, len(ref) - 1), ref), position
    if len(ref) == 1 and len(alt) == 1:
        return substitution(molecule, start, ref, alt), position
    if len(alt) > len(ref) and alt.startswith(ref):
        before = _advance(start, len(ref) - 1)
        # ref is the residue before the insertion, so the index is that of the last residue of the symbol's ref
        return change(molecule, "insertion", before, _advance(before, 1), ref[-1], alt[len(ref):]), position + len(ref) - 1
    return change(molecule, "delins", start, _advance(start, len(ref) - 1), ref, alt), position


def symbol_position(symbol):
    """The number written in an NCBI symbol ('S80I' -> 80), or None."""
    m = re.match(r"^[A-Z]+(-?\d+)", symbol)
    return int(m.group(1)) if m else None


def _advance(pos, n):
    """pos + n on a coordinate without 0 (promoter offsets run ..., -2, -1, then 1 is the first gene base)."""
    out = pos + n
    if pos < 0 <= out:
        out += 1
    return out


def parse_pointfinder(type_gene, codon_pos, ref_codon, res_codons, ref_nuc, promoter):
    """One PointFinder phenotypes.txt row, as one change per residue listed in Res_codon.

    A negative position is upstream of the gene, so it is a nucleotide change of the promoter whatever TypeGene
    says (PointFinder types some promoter rows AA). A row 'AA, position 1, Ref_nuc N, N -> *' is PointFinder's
    way of listing any premature stop codon in the gene: it has no position and is not read as p.Asn1Ter."""
    pos = int(codon_pos)
    if promoter or pos < 0:
        molecule = "promoter"
    elif type_gene == "AA":
        molecule = "protein"
    else:
        molecule = "dna"
    if molecule == "protein" and ref_nuc.strip().upper() == "N" and ref_codon == "N" and res_codons.strip() == "*":
        raise NotationError("any premature stop codon in the gene (no position)")
    out = []
    alts = [a.strip() for a in res_codons.split(",") if a.strip()]
    if not alts:
        raise NotationError(f"PointFinder row with no Res_codon at {codon_pos}")
    for alt in alts:
        if ref_codon == "ins":
            out.append(change(molecule, "insertion", pos, _advance(pos, 1), "", alt.upper()))
        elif ref_codon == "del":
            out.append(change(molecule, "deletion", pos, _advance(pos, len(alt) - 1), alt.upper()))
        else:
            ref = ref_codon.upper() if molecule != "protein" else ref_codon
            if molecule != "protein" and len(ref) == 1 and len(alt) == 1:
                out.append(substitution(molecule, pos, ref, alt.upper()))
            elif molecule == "protein":
                out.append(substitution("protein", pos, ref, _aa(alt) if alt != "*" else "*"))
            else:
                out.append(change(molecule, "delins", pos, _advance(pos, len(ref) - 1), ref, alt.upper()))
    return out


_PF_REQ = re.compile(r"^(?P<gene>.+)_(?P<ref>[A-Za-z*]+)(?P<pos>-?\d+)(?P<alts>[A-Za-z*]+(?:\.[A-Za-z*]+)*)$")


def parse_pointfinder_required(value):
    """'gyrA_S83L.W.A.V,parC_S80I;gyrA_D87N' -> list of alternatives, each a list of (gene, ref, pos, [alts])."""
    groups = []
    for alternative in value.split(";"):
        group = []
        for item in alternative.split(","):
            item = item.strip()
            if not item:
                continue
            m = _PF_REQ.match(item)
            if not m:
                raise NotationError(f"PointFinder required mutation not understood: {item!r}")
            group.append((m.group("gene"), m.group("ref"), int(m.group("pos")), m.group("alts").split(".")))
        if group:
            groups.append(group)
    return groups


def hgvs(c, start=None, end=None):
    """HGVS-like text of a change, on (start, end) if given (anchor coordinates) else its own."""
    s = c["start"] if start is None else start
    e = c["end"] if end is None else end
    k, ref, alt = c["kind"], c["ref"], c["alt"]
    if c["molecule"] == "protein":
        r1 = AA3.get(ref[:1], "") if ref else ""
        r2 = AA3.get(ref[-1:], "") if ref else ""
        if k in ("substitution", "nonsense", "any_substitution", "synonymous"):
            return f"p.{AA3[ref]}{s}{'=' if k == 'synonymous' else AA3[alt]}"
        if k == "frameshift":
            return f"p.{r1}{s}fs"
        if k == "deletion":
            return f"p.{r1}{s}del" if s == e else f"p.{r1}{s}_{AA3.get(ref[-1:], '') if len(ref) == e - s + 1 else ''}{e}del"
        if k == "duplication":
            return f"p.{r1}{s}dup" if s == e else f"p.{s}_{e}dup"
        if k == "insertion":
            return f"p.{r1}{s}_{r2}{e}ins{''.join(AA3.get(a, a) for a in alt)}"
        new = "".join(AA3.get(a, a) for a in alt)
        return f"p.{r1}{s}delins{new}" if s == e else f"p.{s}_{e}delins{new}"
    prefix = {"rna": "n.", "dna": "c.", "promoter": "c."}[c["molecule"]]
    if k in ("substitution", "nonsense", "any_substitution", "synonymous"):
        return f"{prefix}{s}{ref}>{'N' if alt == '?' else alt}" if k != "synonymous" else f"{prefix}{s}="
    if k == "deletion":
        return f"{prefix}{s}del" if s == e else f"{prefix}{s}_{e}del"
    if k == "insertion":
        return f"{prefix}{s}_{e}ins{alt}"
    if k == "duplication":
        return f"{prefix}{s}_{e}dup"
    if k == "frameshift":
        return f"{prefix}{s}fs"
    return f"{prefix}{s}_{e}delins{alt}"


def compact(c):
    """Short notation of a change in its own coordinates, e.g. S83L, Q10*, G129fs, 63_67del, 131_132insSG, C-42T."""
    k, s, e, ref, alt = c["kind"], c["start"], c["end"], c["ref"], c["alt"]
    if k in ("substitution", "nonsense", "any_substitution", "synonymous"):
        return f"{ref}{s}{alt}"
    if k == "frameshift":
        return f"{ref}{s}fs"
    if k == "deletion":
        return f"{ref or ''}{s}del" if s == e else f"{s}_{e}del"
    if k == "insertion":
        return f"{s}_{e}ins{alt}"
    if k == "duplication":
        return f"{s}_{e}dup"
    return f"{s}_{e}delins{alt}"


def change_key(c):
    """Kind and new residues of a change; with equal positions, equal keys are the same change."""
    kind = c["kind"]
    if kind in ("frameshift",):
        return (kind,)
    if kind == "deletion":
        return (kind, c["end"] - c["start"])
    return (kind, c["alt"])
