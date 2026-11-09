"""Write data_manifest.json: for every input file, its source URL, release, retrieval date and sha256.

Usage: python write_manifest.py --inputs inputs.json --data data --out data_manifest.json [--retrieved YYYY-MM-DD]
inputs.json has the keys card, ontology, ncbi and resfinder (see builder/reproduce_2021/inputs.json). The
manifest is checked (CARD release agreement, every file present with its hash) before it is written.
"""
import argparse
import datetime
import json
from pathlib import Path

from functions import validation

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--inputs", required=True)
parser.add_argument("--data", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--retrieved", default=datetime.datetime.now(datetime.timezone.utc).date().isoformat())
args = parser.parse_args()

src = json.loads(Path(args.inputs).read_text())
data = Path(args.data)
card_release = src["card"]["release"]


def entry(db, role, path, url, release, **extra):
    return {"database": db, "role": role, "path": path, "url": url, "release": release,
            "retrieved": args.retrieved, "sha256": validation.sha256_file(data / path), **extra}


ncbi = src["ncbi"]


def _ncbi_name(prefix):
    return next(name for name in ncbi["files"] if name.split(".")[0] == prefix)


rf = src["resfinder"]
files = [
    entry("card", "sequences", "protein_fasta_protein_homolog_model.fasta", src["card"]["url"], card_release,
          card_release=card_release, note="extracted from the release tarball; sha256 is of the extracted file"),
    entry("card", "ontology", "db_metadata/card.metadata.obo", src["ontology"]["url"], src["ontology"]["release"],
          card_release=card_release, note="extracted from the CARD ontology tarball of the same release; sha256 is of the extracted aro.obo"),
    entry("ncbi", "catalogue", "db_metadata/ncbi.metadata.tsv", ncbi["base_url"] + _ncbi_name("ReferenceGeneCatalog"), ncbi["release"]),
    entry("ncbi", "proteins", "db_fastas/ncbi.protein.raw.fasta", ncbi["base_url"] + _ncbi_name("AMRProt"), ncbi["release"]),
    entry("ncbi", "nucleotides", "db_fastas/ncbi.nucl.raw.fasta", ncbi["base_url"] + _ncbi_name("AMR_CDS"), ncbi["release"]),
    entry("resfinder", "phenotypes", "db_metadata/resfinder.metadata.tsv", rf["git_url"], rf["commit"]),
    entry("resfinder", "nucleotides", "db_fastas/resfinder.nucl.fasta", rf["git_url"], rf["commit"],
          note="concatenation of the per-class *.fsa files at this commit"),
]
manifest = {"files": files}
validation.check_manifest(manifest, data)
Path(args.out).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
print(f"wrote {args.out} ({len(files)} files)")
