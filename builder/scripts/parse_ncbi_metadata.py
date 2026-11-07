# %%
import pandas as pd
import json
import os
import argparse

# %%
parser = argparse.ArgumentParser(description='Process input/output files')
parser.add_argument('--meta', type=str, help='input metadata file name')
parser.add_argument('--output', type=str, help='output file name')

args = parser.parse_args()

if args.meta:
	metadata_file = args.meta
else:
	metadata_file =  os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"..", "data", "db_metadata",
		"ncbi.metadata.tsv"
	)

if args.output:
	out_path = args.output
else:
	out_path =  os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"..", "data", "db_metadata",
		"ncbi.metadata.json"
	)

# %%
metadata_df = pd.read_csv(metadata_file, sep="\t")

# %%
metadata_df = metadata_df.query('scope == "core" and subtype !="POINT"' )

metadata_df = metadata_df.sort_values(by=['gene_family', 'allele'])

metadata_df = metadata_df.filter(['product_name', 'allele', 'gene_family', 'subclass', 'refseq_protein_accession', 'refseq_nucleotide_accession'])

# %%
# combine allele and gene family and drop duplicates
metadata_df = (metadata_df
    .assign(allele = metadata_df['allele'].combine_first(metadata_df['gene_family']))
    .drop(columns=['gene_family'])
    # .drop_duplicates(subset=['allele'])
)

# %%
metadata_deduplicated = metadata_df.drop_duplicates(subset=['allele'])

# %%
# aggregated metadata
aggregated_refseq_protein_accessions = (
    metadata_df.groupby(['allele'])['refseq_protein_accession']
    .apply(list)
    .reset_index()
    .rename(columns={"refseq_protein_accession": "aggregated_refseq_protein_accession"})
)
duplicated_refseq_protein_accessions =  (
    aggregated_refseq_protein_accessions['aggregated_refseq_protein_accession']
    .apply(lambda x: x[1:])
)

duplicated_refseq_protein_accessions = (
    aggregated_refseq_protein_accessions
    .assign(duplicated_refseq_protein_accessions = duplicated_refseq_protein_accessions)
    .drop(columns=['aggregated_refseq_protein_accession'])
)
# %%
# aggregated metdata
aggregated_refseq_nucleotide_accessions = (
    metadata_df
    .groupby(['allele'])['refseq_nucleotide_accession']
    .apply(list)
    .reset_index()
    .rename(columns={"refseq_nucleotide_accession": "aggregated_refseq_nucleotide_accessions"})
)

duplicated_refseq_nucleotide_accessions =  (
    aggregated_refseq_nucleotide_accessions['aggregated_refseq_nucleotide_accessions']
    .apply(lambda x: x[1:])
)

duplicated_refseq_nucleotide_accessions = (
    aggregated_refseq_nucleotide_accessions
    .assign(duplicated_refseq_nucleotide_accessions = duplicated_refseq_nucleotide_accessions)
    .drop(columns=['aggregated_refseq_nucleotide_accessions'])
)

# %%
# join together dataframes with deduplicated and duplicates
metadata_deduplicated = (
    metadata_deduplicated
    .join(duplicated_refseq_protein_accessions.set_index('allele'), on='allele')
    .join(duplicated_refseq_nucleotide_accessions.set_index('allele'), on='allele')
)

# %%
metadata = {}
for _, row in metadata_deduplicated.iterrows():
    id = row['refseq_protein_accession']
    metadata[id] = {
        'name': row['allele'],
        'alternative_id': {
            'nucleotide_accession': row['refseq_nucleotide_accession']
        },
        'duplicate_allele_ids': {
            'protein_accessions': row['duplicated_refseq_protein_accessions'],
            'nucleotide_accessions': row['duplicated_refseq_nucleotide_accessions']
        },
        'product': row['product_name'],
        'phenotype': f"confers resistance to subclass {row['subclass']}"
    }

with open(out_path, "w") as out_file:
    out_file.write(
        json.dumps(metadata, sort_keys=True, indent=2)
    )
# %%
