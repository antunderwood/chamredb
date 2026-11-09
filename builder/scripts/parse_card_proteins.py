import os
import argparse
from Bio import SeqIO
from pronto.relationship import Relationship

parser = argparse.ArgumentParser(description='Process input/output files')
parser.add_argument('--input', type=str, help='input file name')
parser.add_argument('--output', type=str, help='output file name')

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
    for record in SeqIO.parse(in_file, "fasta"):

        aro_id, description = (record.id.split("|")[2],record.id.split("|")[3])
        out_file.write(f">{aro_id}\n{str(record.seq)}\n")


