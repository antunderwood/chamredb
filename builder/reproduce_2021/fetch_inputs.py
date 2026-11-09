"""Download the archived 2021 inputs, verify sha256, and lay them out as the builder expects.
Standard library plus git. Usage: python fetch_inputs.py <data_dir>   (data_dir = <build>/data)
"""
import hashlib
import json
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent


def get(url, sha256, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as resp:
        data = resp.read()
    got = hashlib.sha256(data).hexdigest()
    if got != sha256:
        sys.exit(f"sha256 mismatch for {url}: expected {sha256}, got {got}")
    dest.write_bytes(data)


def main(data):
    m = json.loads((HERE / "inputs.json").read_text())
    data.mkdir(parents=True, exist_ok=True)
    tmp = data / "_download"
    get(m["card"]["url"], m["card"]["sha256"], tmp / "card.tar.bz2")
    with tarfile.open(tmp / "card.tar.bz2") as tar:
        member = tar.getmember(m["card"]["member"])
        (data / "protein_fasta_protein_homolog_model.fasta").write_bytes(tar.extractfile(member).read())
    get(m["aro"]["url"], m["aro"]["sha256"], data / "db_metadata" / "card.metadata.obo")
    n = m["ncbi"]
    for name, dest in (("ReferenceGeneCatalog.txt", "db_metadata/ncbi.metadata.tsv"),
                       ("AMRProt", "db_fastas/ncbi.protein.raw.fasta"), ("AMR_CDS", "db_fastas/ncbi.nucl.raw.fasta")):
        get(n["base_url"] + name, n["files"][name], data / dest)
    r, repo = m["resfinder"], tmp / "resfinder_db"
    subprocess.run(["git", "clone", "-q", r["git_url"], str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", r["commit"]], check=True)
    (data / "db_metadata").mkdir(parents=True, exist_ok=True)
    (data / "db_metadata" / "resfinder.metadata.tsv").write_bytes((repo / "phenotypes.txt").read_bytes())
    # same as the builder's `cat resfinder_db/*.fsa` (all.fsa did not exist at this commit)
    (data / "db_fastas").mkdir(parents=True, exist_ok=True)
    with open(data / "db_fastas" / "resfinder.nucl.fasta", "wb") as out:
        for f in sorted(repo.glob("*.fsa")):
            out.write(f.read_bytes())


if __name__ == "__main__":
    main(Path(sys.argv[1]))
