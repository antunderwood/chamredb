"""Plan a graph release (ticket #8, ADR 0005). Standard library only.

  plan_release.py needs-build --previous PREV_MANIFEST.json|none --inputs INPUTS.json
      prints `build=true|false` and `reasons=...` lines (for GITHUB_OUTPUT)
  plan_release.py notes --tag TAG --manifest MANIFEST.json --report REPORT.md --out NOTES.md
"""
import argparse
import datetime
import json
import sys


def tag_for(date):
    return f"graph-{date:%Y.%m.%d}"


def _short(commit):
    return commit[:7]


def releases_from_manifest(manifest):
    out = {}
    for f in manifest["files"]:
        if f["database"] == "card":
            out["card"] = f["card_release"]
        elif f["database"] == "ncbi":
            out["ncbi"] = f["release"]
        elif f["database"] == "resfinder":
            out["resfinder"] = _short(f["release"])
    return out


def releases_from_inputs(inputs):
    return {"card": inputs["card"]["release"], "ncbi": inputs["ncbi"]["release"],
            "resfinder": _short(inputs["resfinder"]["commit"])}


def needs_build(previous_manifest, inputs):
    """(build?, reasons): build when there is no previous release or any database release differs."""
    if previous_manifest is None:
        return True, ["no previous release"]
    old, new = releases_from_manifest(previous_manifest), releases_from_inputs(inputs)
    reasons = [f"{db}: {old.get(db)} -> {new[db]}" for db in ("card", "ncbi", "resfinder") if old.get(db) != new[db]]
    return bool(reasons), reasons


def release_notes(tag, manifest, report):
    r = releases_from_manifest(manifest)
    return "\n".join([
        f"# {tag}", "",
        "Graph of Source database matches (RBH and one-way hits) built from:", "",
        f"- CARD {r['card']} (data and ontology from the same release)",
        f"- NCBI {r['ncbi']}",
        f"- ResFinder db commit {r['resfinder']}", "",
        "Assets: `graph.json`, `data_manifest.json` (URL, release, retrieval date and sha256 per input file), `diff_report.md`, and the three database metadata files (`card`, `ncbi`, `resfinder` `.metadata.json`) that go with the graph.", "",
        "## Changes since the previous release", "", report.strip(), ""])


def main(argv=None):
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    n = sub.add_parser("needs-build")
    n.add_argument("--previous", required=True)
    n.add_argument("--inputs", required=True)
    o = sub.add_parser("notes")
    o.add_argument("--tag", required=True)
    o.add_argument("--manifest", required=True)
    o.add_argument("--report", required=True)
    o.add_argument("--out", required=True)
    args = p.parse_args(argv)
    if args.cmd == "needs-build":
        prev = None if args.previous == "none" else json.load(open(args.previous))
        build, reasons = needs_build(prev, json.load(open(args.inputs)))
        print(f"build={'true' if build else 'false'}")
        print("reasons=" + "; ".join(reasons))
    else:
        open(args.out, "w").write(release_notes(args.tag, json.load(open(args.manifest)), open(args.report).read()))


if __name__ == "__main__":
    main()
