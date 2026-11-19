"""Write the mutation part of a release diff report.

Usage: python diff_variants.py OLD_variants.json NEW_variants.json --out variants_diff.md
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from functions import variants_diff  # noqa: E402

p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
p.add_argument("old")
p.add_argument("new")
p.add_argument("--out", required=True)
args = p.parse_args()
d = variants_diff.diff(json.loads(Path(args.old).read_text()), json.loads(Path(args.new).read_text()))
Path(args.out).write_text(variants_diff.markdown(d))
print(f"wrote {args.out}: {len(d['added'])} added, {len(d['removed'])} removed, {len(d['anchor_changes'])} anchor changes")
