import os
import argparse
from functions import validation
from Bio import SeqIO
from pronto.relationship import Relationship

parser = argparse.ArgumentParser(description='Process input/output files')
parser.add_argument('--input', type=str, help='input file name')
parser.add_argument('--output', type=str, help='output file name')

parser.add_argument('--bounds', type=str, help='JSON file with min and max protein counts per database (optional)')

args = parser.parse_args()

if args.input:
	in_path = args.input
else:
	in_path = os.path.join(
		        os.path.dirname(os.path.abspath(__file__)),
			    "..", "data", "card_db",
				"protein_fasta_protein_homolog_model.fasta"
	)

if args.output:
	out_path = args.output
else:
	out_path = os.path.join(
		       os.path.dirname(os.path.abspath(__file__)),
			    "..", "data", "db_fastas",
				"card.protein.fasta"
	)

with open(in_path) as in_file, open(out_path, "w") as out_file:
    aro_ids = []
    for record in SeqIO.parse(in_file, "fasta"):
        fields = record.id.split("|")
        if len(fields) < 4:
            raise ValueError(f"{in_path}: unexpected CARD header (need gb|accession|ARO|name): {record.id}")
        aro_id, description = (fields[2], fields[3])
        aro_ids.append(aro_id)
        out_file.write(f">{aro_id}\n{str(record.seq)}\n")

validation.check_unique(aro_ids, 'card proteins')
validation.check_count(len(aro_ids), validation.load_bounds(args.bounds).get('card'), 'card proteins')



