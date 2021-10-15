from typing import Sequence
import orfipy_core 
import re
from Bio import SeqIO
from Bio.Seq import Seq

with open("/data/mAMRite/db_fastas/resfinder.nucl.fasta") as in_file, open("/data/mAMRite/db_fastas/resfinder.protein.fasta", "w") as out_file:
    genes_added = []
    for record in SeqIO.parse(in_file, "fasta"):
        id = record.id.split("_")[0]
        if id in genes_added:
            continue
        genes_added.append(id)
        sequence = str(record.seq)
        max_length = 0
        orf_start = None
        orf_end = None
        orf_description = None
        orf_frame = None
        # find best ORF
        for start,stop,strand,description in orfipy_core.orfs(sequence, partial3 = True):
            if stop - start > max_length:
                max_length = stop - start
                orf_start = start
                orf_end = stop
                orf_frame = int(re.search('ORF_frame=([-\d]+)', description).group(1))
        if (orf_frame < 1):
            protein_sequence = record.seq[orf_start:orf_end].reverse_complement().translate()
        else:
            protein_sequence = record.seq[orf_start:orf_end].translate()

        out_file.write(f">{id}\n{protein_sequence}\n")

