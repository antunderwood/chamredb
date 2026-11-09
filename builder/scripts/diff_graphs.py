"""Compare two chamredb graph JSON files and write a release-review report.

Usage: python diff_graphs.py OLD.json NEW.json [--out report.md] [--json report.json]
Informational: always exits 0. Use scripts/functions/graph_compare.py for a pass/fail comparison.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "functions"))
from graph_diff import diff, report  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("old")
parser.add_argument("new")
parser.add_argument("--out", help="write the Markdown report here (default: print)")
parser.add_argument("--json", help="also write the full diff as JSON here")
args = parser.parse_args()

result = diff(json.load(open(args.old)), json.load(open(args.new)))
text = report(result)
if args.out:
    open(args.out, "w").write(text)
else:
    print(text)
if args.json:
    json.dump(result, open(args.json, "w"), indent=1, default=list)
