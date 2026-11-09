import pandas as pd
import os
import json
import argparse
import tempfile
from functions import validation

parser = argparse.ArgumentParser(description='Process input/output files')
parser.add_argument('--meta', type=str, help='input metadata file name')
parser.add_argument('--nucl', type=str, required=True, help='resfinder nucleotide FASTA; the allele that comes first for each gene is the one whose sequence is used')
parser.add_argument('--output', type=str, help='output file name')

args = parser.parse_args()

if args.meta:
	metadata_file = args.meta
else:
	metadata_file =  os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"..", "data", "db_metadata",
		"resfinder.metadata.tsv"
	)

if args.output:
	out_path = args.output
else:
	out_path =  os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"..", "data", "db_metadata",
		"resfinder.metadata.json"
	)

with open(metadata_file) as metadata_fh:
	with tempfile.NamedTemporaryFile('w') as metadata_parsed:
		for line in metadata_fh:
			# get rid of bad lines
			if line.count('\t') > 6:
				line = line.replace('\t\t', '\t')
			metadata_parsed.write(line)
		metadata_parsed.flush()
		metadata_df = pd.read_csv(metadata_parsed.name, sep = "\t", on_bad_lines="warn")
		validation.require_columns(metadata_df.columns, ['Gene_accession no.', 'Phenotype'], metadata_file)
		metadata_df = metadata_df[['Gene_accession no.', 'Phenotype']]

		# the phenotype of a gene is the phenotype of the allele whose sequence is used (the first record of the gene in the
		# nucleotide FASTA, as in resfinder_orfs_to_proteins.py); before ticket #14 the last row of the table won
		used_allele = {}
		with open(args.nucl) as nucl_fh:
			for line in nucl_fh:
				if line.startswith('>'):
					allele = line[1:].split()[0]
					used_allele.setdefault(allele.split("_")[0], allele)
		phenotype_by_allele = dict(zip(metadata_df['Gene_accession no.'], metadata_df['Phenotype']))
		metadata = {}
		for _, row in metadata_df.iterrows():
			name = row['Gene_accession no.'].split("_")[0]
			if name not in used_allele:
				metadata.setdefault(name, {'name': name, 'phenotype': f"confers resistance to {row['Phenotype']}"})
		for name, allele in used_allele.items():
			if allele not in phenotype_by_allele:
				raise ValueError(f"{metadata_file}: no phenotype row for {allele}, the allele whose sequence is used for {name}")
			metadata[name] = {
				'name': name,
				'phenotype': f"confers resistance to {phenotype_by_allele[allele]}"
			}

		with open(out_path, "w") as out_file:
			out_file.write(
				json.dumps(metadata, sort_keys=True, indent=2)
    )



