"""Find, download, verify and cache the graph (ADR 0004). Standard library only.

The package pins a graph release in data/graph_release.json: data_release, sha256 and, once a release asset
exists, its url. With no url the bundled graph.json is used. A pinned url is downloaded once into the cache
directory, verified against the sha256, and reused offline afterwards.
"""
import hashlib
import json
import os
import urllib.request
from pathlib import Path

PACKAGE_DATA = Path(__file__).resolve().parent / "data"
CUSTOM_RELEASE = "unknown (custom graph)"


def default_cache_dir():
    """CHAMREDB_CACHE_DIR, else $XDG_CACHE_HOME/chamredb, else ~/.cache/chamredb (same on Linux and macOS)."""
    if os.environ.get("CHAMREDB_CACHE_DIR"):
        return Path(os.environ["CHAMREDB_CACHE_DIR"])
    return Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "chamredb"


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_graph(path=None, pin_file=None, cache_dir=None):
    """Return (graph file path, data release label). An explicit path wins and is not verified."""
    if path:
        return Path(path), CUSTOM_RELEASE
    pin = json.loads(Path(pin_file or PACKAGE_DATA / "graph_release.json").read_text())
    if not pin.get("url"):
        bundled = PACKAGE_DATA / "graph.json"
        if _sha256(bundled) != pin["sha256"]:
            raise ValueError(f"bundled graph.json sha256 does not match graph_release.json ({pin['data_release']})")
        return bundled, pin["data_release"]
    cached = Path(cache_dir or default_cache_dir()) / pin["sha256"] / "graph.json"
    if cached.is_file() and _sha256(cached) == pin["sha256"]:
        return cached, pin["data_release"]
    with urllib.request.urlopen(pin["url"]) as resp:
        data = resp.read()
    if hashlib.sha256(data).hexdigest() != pin["sha256"]:
        raise ValueError(f"downloaded graph from {pin['url']} has the wrong sha256; not cached")
    cached.parent.mkdir(parents=True, exist_ok=True)
    cached.write_bytes(data)
    return cached, pin["data_release"]


def resolve_metadata(database, graph_path=None, pin_file=None):
    """Metadata dict of one database from the files pinned next to the graph (ticket #33).

    The pin lists `metadata: {database: {file, sha256}}`, with `file` relative to the pin file. Returns None for a
    custom graph (its metadata is unknown) or when the pin lists no file for the database; a file whose sha256
    differs from the pin is an error.
    """
    if graph_path:
        return None
    pin_path = Path(pin_file or PACKAGE_DATA / "graph_release.json")
    entry = json.loads(pin_path.read_text()).get("metadata", {}).get(database)
    if entry is None:
        return None
    path = pin_path.parent / entry["file"]
    if _sha256(path) != entry["sha256"]:
        raise ValueError(f"{path} sha256 does not match the pin ({database} metadata)")
    return json.loads(path.read_text())

