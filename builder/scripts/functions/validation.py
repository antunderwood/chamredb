"""Strict input checks for the builder parsers (spec #1, predicates 4 to 7). Standard library only.

Every check raises ValueError with a message naming the source and the problem, so a changed upstream
format stops the build instead of producing plausible but wrong data.
"""
import hashlib
import json
from pathlib import Path


def require_columns(columns, required, source):
    missing = [c for c in required if c not in set(columns)]
    if missing:
        raise ValueError(f"{source}: missing required columns {missing}; found {sorted(columns)}")


def load_bounds(path):
    """{'card': {'min': int, 'max': int}, ...} from a JSON file, or {} when no file is given."""
    return json.loads(Path(path).read_text()) if path else {}


def check_count(n, bounds, what):
    """bounds is {'min': .., 'max': ..} (either may be missing) or None (no check)."""
    if not bounds:
        return
    lo, hi = bounds.get("min"), bounds.get("max")
    if (lo is not None and n < lo) or (hi is not None and n > hi):
        raise ValueError(f"{what}: {n} outside the stated bounds [{lo}, {hi}]")


def check_unique(ids, what):
    seen, dup = set(), set()
    for i in ids:
        (dup if i in seen else seen).add(i)
    if dup:
        raise ValueError(f"{what}: {len(dup)} duplicate ids, e.g. {sorted(dup)[:3]}")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_manifest(manifest, data_dir=None):
    """Predicate 7: the CARD sequences and ontology declare one CARD release (`card_release`; the ontology
    is a separate download with no release tag, so the pairing is declared, with a `note`). With data_dir,
    also check that every listed file exists with the recorded sha256."""
    files = manifest["files"]
    card = {f["role"]: f.get("card_release") for f in files if f["database"] == "card"}
    if not card or None in card.values():
        raise ValueError(f"CARD manifest entries must all declare card_release: {card}")
    if len(set(card.values())) > 1:
        raise ValueError(f"CARD files come from different releases: {card}")
    for f in files:
        for key in ("database", "role", "url", "release", "retrieved", "sha256"):
            if not f.get(key):
                raise ValueError(f"manifest entry {f.get('path')}: missing {key}")
        if data_dir is not None:
            path = Path(data_dir) / f["path"]
            if not path.is_file():
                raise ValueError(f"manifest lists {f['path']} but it does not exist")
            if sha256_file(path) != f["sha256"]:
                raise ValueError(f"{f['path']}: sha256 differs from the manifest")
