"""Download the pinned inputs, verify sha256, and lay them out as the builder expects.
Standard library plus git. Usage: python fetch_inputs.py <data_dir> [inputs.json]   (data_dir = <build>/data;
inputs.json defaults to the 2021 inputs next to this file)
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


def main(data, inputs=None):
    m = json.loads(Path(inputs or HERE / "inputs.json").read_text())
    data.mkdir(parents=True, exist_ok=True)
    tmp = data / "_download"
    get(m["card"]["url"], m["card"]["sha256"], tmp / "card.tar.bz2")
    with tarfile.open(tmp / "card.tar.bz2") as tar:
        member = tar.getmember(m["card"]["member"])
        (data / "protein_fasta_protein_homolog_model.fasta").write_bytes(tar.extractfile(member).read())
        if m["card"].get("variants_member"):   # card.json, for the mutation comparison
            (data / "mutations").mkdir(parents=True, exist_ok=True)
            (data / "mutations" / "card.json").write_bytes(tar.extractfile(tar.getmember(m["card"]["variants_member"])).read())
    o = m["ontology"]
    get(o["url"], o["sha256"], tmp / "ontology.tar.bz2")
    with tarfile.open(tmp / "ontology.tar.bz2") as tar:
        obo = tar.extractfile(tar.getmember(o["member"])).read()
    if hashlib.sha256(obo).hexdigest() != o["member_sha256"]:
        sys.exit(f"sha256 mismatch for {o['member']} in {o['url']}")
    (data / "db_metadata").mkdir(parents=True, exist_ok=True)
    (data / "db_metadata" / "card.metadata.obo").write_bytes(obo)
    n = m["ncbi"]
    # file names changed between releases (AMRProt in 2021, AMRProt.fa in 2026); the role is the name prefix
    for name, sha in n["files"].items():
        dest = {"ReferenceGeneCatalog": "db_metadata/ncbi.metadata.tsv", "AMRProt": "db_fastas/ncbi.protein.raw.fasta",
                "AMR_CDS": "db_fastas/ncbi.nucl.raw.fasta"}[name.split(".")[0]]
        get(n["base_url"] + name, sha, data / dest)
    r, repo = m["resfinder"], tmp / "resfinder_db"
    subprocess.run(["git", "clone", "-q", r["git_url"], str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", r["commit"]], check=True)
    (data / "db_metadata").mkdir(parents=True, exist_ok=True)
    (data / "db_metadata" / "resfinder.metadata.tsv").write_bytes((repo / "phenotypes.txt").read_bytes())
    # concatenate the per-class files; all.fsa (present from 2022, a copy of all the others) must never be included:
    # it would double every record (ticket #14). When it exists it is used only to cross-check the per-class union.
    (data / "db_fastas").mkdir(parents=True, exist_ok=True)
    per_class = [f for f in sorted(repo.glob("*.fsa")) if f.name != "all.fsa"]
    with open(data / "db_fastas" / "resfinder.nucl.fasta", "wb") as out:
        for f in per_class:
            out.write(f.read_bytes())
    all_fsa = repo / "all.fsa"
    if all_fsa.is_file():
        ids = lambda paths: {line[1:].split()[0] for p in paths for line in open(p) if line.startswith(">")}
        if ids([all_fsa]) != ids(per_class):
            sys.exit("resfinder all.fsa does not match the union of the per-class files")

    fetch_mutation_inputs(data, m, tmp)


def fetch_mutation_inputs(data, m, tmp):
    """NCBI point-mutation files and pointfinder_db, when the inputs file lists them."""
    nm = m.get("ncbi_mutations")
    if nm:
        for name, sha in nm["files"].items():
            get(nm["base_url"] + name, sha, data / "mutations" / "ncbi" / name)
    pf = m.get("pointfinder")
    if pf:
        repo = tmp / "pointfinder_db"
        subprocess.run(["git", "clone", "-q", pf["git_url"], str(repo)], check=True)
        subprocess.run(["git", "-C", str(repo), "checkout", "-q", pf["commit"]], check=True)
        dest = data / "mutations" / "pointfinder_db"
        lines = []
        for f in sorted(p for p in repo.rglob("*") if p.is_file() and ".git" not in p.relative_to(repo).parts):
            rel = f.relative_to(repo)
            (dest / rel).parent.mkdir(parents=True, exist_ok=True)
            content = f.read_bytes()
            (dest / rel).write_bytes(content)
            lines.append(f"{hashlib.sha256(content).hexdigest()}  {rel.as_posix()}\n")
        (data / "mutations" / "pointfinder_db.sha256").write_text("".join(lines))


if __name__ == "__main__":
    main(Path(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else None)
