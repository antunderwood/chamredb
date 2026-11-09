"""Resolve the current upstream releases into a pinned inputs file (ticket #8). Standard library plus git.

Usage: python resolve_inputs.py OUT.json
Downloads the current CARD data and ontology tarballs (checked to be the same release), the NCBI AMRFinderPlus
database files (versioned directory, not `latest`) and the ResFinder db HEAD commit, hashes them, and writes the
same schema as builder/refresh/inputs_*.json, so the build is reproducible from the file.
"""
import hashlib
import json
import re
import subprocess
import sys
import tarfile
import io
import urllib.request

CARD = "https://card.mcmaster.ca"
NCBI = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database"
RESFINDER = "https://bitbucket.org/genomicepidemiology/resfinder_db"


def parse_card_version(card_json_bytes):
    version = json.loads(card_json_bytes).get("_version")
    if not version:
        raise ValueError("card.json has no _version")
    return version


def pick_ncbi_dir(version, listing_html, exists):
    """Newest major directory (for example 4.2) that holds the database version."""
    dirs = re.findall(r'href="(\d+\.\d+)/"', listing_html)
    for d in sorted(dirs, key=lambda x: tuple(map(int, x.split("."))), reverse=True):
        if exists(d, version):
            return d
    raise ValueError(f"no database directory holds NCBI version {version}")


def _get(url):
    with urllib.request.urlopen(url) as resp:
        return resp.read()


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _member(tar_bytes, name):
    with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:bz2") as tar:
        return tar.extractfile(tar.getmember(name)).read()


def resolve():
    latest_data = _get(f"{CARD}/latest/data")
    version = parse_card_version(_member(latest_data, "./card.json"))
    data_url = f"{CARD}/download/0/broadstreet-v{version}.tar.bz2"
    ont_url = f"{CARD}/download/5/ontology-v{version}.tar.bz2"
    data, ontology = _get(data_url), _get(ont_url)
    if _sha(data) != _sha(latest_data):
        raise ValueError(f"{data_url} differs from {CARD}/latest/data")
    if _sha(ontology) != _sha(_get(f"{CARD}/latest/ontology")):
        raise ValueError(f"{ont_url} differs from {CARD}/latest/ontology: ontology and data are not the same release")
    obo = _member(ontology, "./aro.obo")

    ncbi_version = _get(f"{NCBI}/latest/version.txt").decode().strip()
    listing = _get(f"{NCBI}/").decode()

    def exists(d, v):
        try:
            urllib.request.urlopen(urllib.request.Request(f"{NCBI}/{d}/{v}/version.txt", method="HEAD")).close()
            return True
        except Exception:
            return False

    ncbi_dir = pick_ncbi_dir(ncbi_version, listing, exists)
    base = f"{NCBI}/{ncbi_dir}/{ncbi_version}/"
    files = {n: _sha(_get(base + n)) for n in ("ReferenceGeneCatalog.txt", "AMRProt.fa", "AMR_CDS.fa")}

    head = subprocess.run(["git", "ls-remote", RESFINDER, "HEAD"], check=True, capture_output=True, text=True).stdout.split()[0]
    return {
        "purpose": "Current upstream releases resolved by builder/release/resolve_inputs.py",
        "card": {"release": version, "url": data_url, "sha256": _sha(data), "member": "./protein_fasta_protein_homolog_model.fasta"},
        "ontology": {"release": version, "url": ont_url, "sha256": _sha(ontology), "member": "./aro.obo", "member_sha256": _sha(obo)},
        "ncbi": {"release": f"AMRFinderPlus database {ncbi_dir}, {ncbi_version}", "base_url": base, "files": files},
        "resfinder": {"release": f"resfinder_db commit {head[:7]}", "git_url": RESFINDER, "commit": head},
    }


if __name__ == "__main__":
    out = resolve()
    json.dump(out, open(sys.argv[1], "w"), indent=2)
    print(f"resolved CARD {out['card']['release']}, {out['ncbi']['release']}, ResFinder {out['resfinder']['commit'][:7]}")
