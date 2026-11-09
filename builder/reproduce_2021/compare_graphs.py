"""Compare a rebuilt graph with the shipped one (ADR 0003 step 1). Standard library only.

Pass criterion (spec #1): identical node set, edge set and edge types; identity and coverage within
TOLERANCE. Node attribute differences are reported as findings but do not fail the comparison.
Usage: python compare_graphs.py SHIPPED.json REBUILT.json [REPORT.md]   (exit 0 if ok, 1 if not)
"""
import json
import sys

TOLERANCE = 1e-3


def _norm(value):
    return sorted(value) if isinstance(value, list) and all(isinstance(v, (str, int, float)) for v in value) else value


def compare(shipped, rebuilt, tol=TOLERANCE):
    sn = {n["id"]: n for n in shipped["nodes"]}
    rn = {n["id"]: n for n in rebuilt["nodes"]}
    se = {(e["source"], e["target"]): e for e in shipped["links"]}
    re_ = {(e["source"], e["target"]): e for e in rebuilt["links"]}
    r = {
        "shipped_nodes": len(sn), "rebuilt_nodes": len(rn), "shipped_edges": len(se), "rebuilt_edges": len(re_),
        "nodes_lost": sorted(sn.keys() - rn.keys()), "nodes_gained": sorted(rn.keys() - sn.keys()),
        "edges_lost": sorted(se.keys() - re_.keys()), "edges_gained": sorted(re_.keys() - se.keys()),
        "type_changed": [], "identity_off": [], "coverage_off": [], "attribute_diffs": [], "findings": [],
    }
    for key in sorted(se.keys() & re_.keys()):
        a, b = se[key], re_[key]
        if a["type"] != b["type"]:
            r["type_changed"].append((key, a["type"], b["type"]))
        if abs(a["identity"] - b["identity"]) > tol:
            r["identity_off"].append((key, a["identity"], b["identity"]))
        if abs(a["coverage"] - b["coverage"]) > tol:
            r["coverage_off"].append((key, a["coverage"], b["coverage"]))
    for nid in sorted(sn.keys() & rn.keys()):
        for field in sorted(sn[nid].keys() | rn[nid].keys()):
            if field != "id" and _norm(sn[nid].get(field)) != _norm(rn[nid].get(field)):
                r["attribute_diffs"].append((nid, field, sn[nid].get(field), rn[nid].get(field)))
    r["ok"] = not any(r[k] for k in ("nodes_lost", "nodes_gained", "edges_lost", "edges_gained",
                                     "type_changed", "identity_off", "coverage_off"))
    if r["attribute_diffs"]:
        r["findings"].append(f"{len(r['attribute_diffs'])} node attribute differences (reported, not failing)")
    return r


def report(r, limit=20):
    lines = [f"# 2021 graph reproduction: {'MATCH' if r['ok'] else 'MISMATCH'}", "",
             f"- Nodes: shipped {r['shipped_nodes']}, rebuilt {r['rebuilt_nodes']}",
             f"- Edges: shipped {r['shipped_edges']}, rebuilt {r['rebuilt_edges']}"]
    for key, label in (("nodes_lost", "Nodes only in shipped"), ("nodes_gained", "Nodes only in rebuilt"),
                       ("edges_lost", "Edges only in shipped"), ("edges_gained", "Edges only in rebuilt"),
                       ("type_changed", "Edges with changed type"), ("identity_off", f"Identity off by more than {TOLERANCE}"),
                       ("coverage_off", f"Coverage off by more than {TOLERANCE}"), ("attribute_diffs", "Node attribute differences")):
        lines.append(f"- {label}: {len(r[key])}")
    for key in ("nodes_lost", "nodes_gained", "edges_lost", "edges_gained", "type_changed", "identity_off", "coverage_off", "attribute_diffs"):
        if r[key]:
            lines += ["", f"## {key} (first {limit})"] + [f"- {x}" for x in r[key][:limit]]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    result = compare(json.load(open(sys.argv[1])), json.load(open(sys.argv[2])))
    text = report(result)
    print(text)
    if len(sys.argv) > 3:
        open(sys.argv[3], "w").write(text)
    sys.exit(0 if result["ok"] else 1)
