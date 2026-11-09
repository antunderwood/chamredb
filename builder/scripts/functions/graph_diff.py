"""Diff of two chamredb graphs for release review (ticket #7). Standard library only.

diff() extends graph_compare.compare with ID continuity: for every node id in the old graph, is it kept,
reassigned (same id, different name) or absent (with the new ids in the same database that carry its name).
"""
from graph_compare import TOLERANCE, compare


def _database(node_id):
    return node_id.split(":", 1)[0]


def _continuity(old, new):
    old_nodes = {n["id"]: n for n in old["nodes"]}
    new_nodes = {n["id"]: n for n in new["nodes"]}
    by_name = {}
    for i, n in new_nodes.items():
        by_name.setdefault((_database(i), n.get("name")), []).append(i)
    kept, reassigned, absent = [], [], []
    by_db = {}
    for i in sorted(old_nodes):
        counts = by_db.setdefault(_database(i), {"kept": 0, "reassigned": 0, "absent": 0})
        if i not in new_nodes:
            absent.append((i, sorted(by_name.get((_database(i), old_nodes[i].get("name")), []))))
            counts["absent"] += 1
        elif old_nodes[i].get("name") != new_nodes[i].get("name"):
            reassigned.append((i, old_nodes[i].get("name"), new_nodes[i].get("name")))
            counts["reassigned"] += 1
        else:
            kept.append(i)
            counts["kept"] += 1
    return {"kept": kept, "reassigned": reassigned, "absent": absent, "by_database": by_db}


def diff(old, new, tol=TOLERANCE):
    d = compare(old, new, tol)
    d["continuity"] = _continuity(old, new)
    return d


def report(d, limit=20):
    c = d["continuity"]
    lines = ["# Graph diff", "",
             f"- Nodes: old {d['shipped_nodes']}, new {d['rebuilt_nodes']}",
             f"- Edges: old {d['shipped_edges']}, new {d['rebuilt_edges']}",
             f"- Nodes gained: {len(d['nodes_gained'])}", f"- Nodes lost: {len(d['nodes_lost'])}",
             f"- Edges gained: {len(d['edges_gained'])}", f"- Edges lost: {len(d['edges_lost'])}",
             f"- Edges with changed type: {len(d['type_changed'])}",
             f"- Edges with identity changed by more than {TOLERANCE}: {len(d['identity_off'])}",
             f"- Edges with coverage changed by more than {TOLERANCE}: {len(d['coverage_off'])}",
             f"- Node attribute differences: {len(d['attribute_diffs'])}", "",
             "## ID continuity (ids in the old graph)", "",
             f"- IDs kept: {len(c['kept'])}",
             f"- IDs absent in the new graph: {len(c['absent'])}",
             f"- IDs reassigned (same id, different name): {len(c['reassigned'])}", "",
             "| database | kept | reassigned | absent |", "| --- | --- | --- | --- |"]
    lines += [f"| {db} | {n['kept']} | {n['reassigned']} | {n['absent']} |" for db, n in sorted(c["by_database"].items())]
    if c["absent"]:
        lines += ["", f"### Absent ids (first {limit}; successor = same name in the same database)"]
        lines += [f"- {i} -> {', '.join(s) if s else '(no successor)'}" for i, s in c["absent"][:limit]]
    if c["reassigned"]:
        lines += ["", f"### Reassigned ids (first {limit})"]
        lines += [f"- {i}: {a} -> {b}" for i, a, b in c["reassigned"][:limit]]
    for key, label in (("type_changed", "Edges with changed type"), ("identity_off", "Identity changes"),
                       ("coverage_off", "Coverage changes"), ("nodes_lost", "Nodes lost"), ("nodes_gained", "Nodes gained"),
                       ("edges_lost", "Edges lost"), ("edges_gained", "Edges gained"), ("attribute_diffs", "Node attribute differences")):
        if d[key]:
            lines += ["", f"### {label} (first {limit})"] + [f"- {x}" for x in d[key][:limit]]
    return "\n".join(lines) + "\n"
