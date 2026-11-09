"""Testable invariants of a chamredb graph (spec #1, predicates 1 to 4). Standard library only.

check_graph returns {violation_name: [examples]}; an empty dict means all graph predicates hold.
Names start with the predicate: p1 (ids and endpoints), p2 (RBH reciprocal), p3 (identity and coverage range).
"""
DATABASES = {"card", "ncbi", "resfinder"}


def check_graph(graph):
    out = {}

    def add(name, example):
        out.setdefault(name, []).append(example)

    ids = {n["id"] for n in graph["nodes"]}
    for i in sorted(ids):
        db, sep, rest = i.partition(":")
        if not sep or not rest:
            add("p1_node_id_without_database", i)
        elif db not in DATABASES:
            add("p1_unknown_database", i)
    types = {}
    for e in graph["links"]:
        s, t = e["source"], e["target"]
        types[(s, t)] = e["type"]
        for end in (s, t):
            if end not in ids:
                add("p1_edge_endpoint_not_a_node", (s, t))
        if s.split(":")[0] == t.split(":")[0]:
            add("p1_edge_within_one_database", (s, t))
        if not 0 < e["identity"] <= 1:
            add("p3_identity_out_of_range", (s, t, e["identity"]))
        if not 0 < e["coverage"] <= 1:
            add("p3_coverage_out_of_range", (s, t, e["coverage"]))
    for (s, t), ty in types.items():
        if ty == "RBH" and types.get((t, s)) != "RBH":
            add("p2_rbh_not_reciprocal", (s, t))
    return out


def duplicate_fasta_ids(path):
    """Ids (first word after '>') that occur more than once, sorted (predicate 4)."""
    seen, dup = set(), set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                i = line[1:].split()[0]
                (dup if i in seen else seen).add(i)
    return sorted(dup)
