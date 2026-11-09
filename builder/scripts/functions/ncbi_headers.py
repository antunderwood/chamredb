"""AMRProt FASTA header layouts (ticket #13). Standard library only."""


def parse_header(header):
    """Return (refseq protein accession, allele name) from an AMRProt record id.

    11 fields (2021): gi|accession|1|1|allele|...    10 fields (2026): accession|1|1|allele|...
    """
    fields = header.split("|")
    if len(fields) == 11:
        return fields[1], fields[4]
    if len(fields) == 10:
        return fields[0], fields[3]
    raise ValueError(f"unexpected AMRProt header with {len(fields)} fields: {header[:80]}")
