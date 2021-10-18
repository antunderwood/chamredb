# Python dependecies pre-load
import  os
import  sys
import  re
import  pandas as pd
from    os.path import join as join

shell.executable("bash")

DATA       = 'data'

rule all:
	input:
		# CARD
		join(DATA, 'db_metadata', 'card.metadata.obo'),
		join(DATA, 'protein_fasta_protein_homolog_model.fasta'),
		join(DATA, 'nucleotide_fasta_protein_homolog_model.fasta'),
		join(DATA, 'db_fastas', 'card.protein.fasta'),
		join(DATA, 'db_metadata', 'card.metadata.json'),
		# NCBI
		join(DATA, 'db_metadata', 'ncbi.metadata.tsv'),
		join(DATA, 'db_fastas', 'ncbi.protein.raw.fasta'),
		join(DATA, 'db_fastas', 'ncbi.nucl.raw.fasta'),
		join(DATA, 'db_fastas', 'ncbi.protein.fasta'),
		join(DATA, 'db_metadata', 'ncbi.metadata.json'),
		# Resfinder
		join(DATA, 'db_metadata', 'resfinder.metadata.tsv'),
		join(DATA, 'db_fastas', 'resfinder.nucl.fasta'),
		join(DATA, 'db_fastas', 'resfinder.protein.fasta'),
		join(DATA, 'db_metadata', 'resfinder.metadata.json'),
		# RBH
		expand(join(DATA, 'mmseqs2_results', '{db1}_vs_{db2}.rbh.tsv'), db1='resfinder', db2=['card','ncbi']),
		expand(join(DATA, 'mmseqs2_results', '{db1}_vs_{db2}.rbh.tsv'), db1='card', db2=['resfinder','ncbi']),
		expand(join(DATA, 'mmseqs2_results', '{db1}_vs_{db2}.rbh.tsv'), db1='ncbi', db2=['card','resfinder']),
		expand(join(DATA, 'mmseqs2_DBs', '{db1}.protein'), db1=['resfinder','card','ncbi']),
		expand(join(DATA, 'mmseqs2_search_DBs', '{db1}_vs_{db2}.search.index'), db1='resfinder', db2=['card','ncbi']),
		expand(join(DATA, 'mmseqs2_search_DBs', '{db1}_vs_{db2}.search.index'), db1='card', db2=['resfinder','ncbi']),
		expand(join(DATA, 'mmseqs2_search_DBs', '{db1}_vs_{db2}.search.index'), db1='ncbi', db2=['card','resfinder']),
		expand(join(DATA, 'mmseqs2_results', '{db1}_vs_{db2}.search.tsv'), db1='resfinder', db2=['card','ncbi']),
		expand(join(DATA, 'mmseqs2_results', '{db1}_vs_{db2}.search.tsv'), db1='card', db2=['resfinder','ncbi']),
		expand(join(DATA, 'mmseqs2_results', '{db1}_vs_{db2}.search.tsv'), db1='ncbi', db2=['card','resfinder']),

rule card_download:
	output:
		obo  = join(DATA, 'db_metadata', 'card.metadata.obo'),
		data = join(DATA, 'card_db', 'data'),
		prot = join(DATA, 'protein_fasta_protein_homolog_model.fasta'),
		nucl = join(DATA, 'nucleotide_fasta_protein_homolog_model.fasta')
	shell:
		"""
		wget https://github.com/arpcard/aro/raw/master/src/ontology/aro.obo -O {output.obo}
		wget https://card.mcmaster.ca/latest/data -O {output.data}
		tar -xvf {output.data} ./protein_fasta_protein_homolog_model.fasta
		tar -xvf {output.data} ./nucleotide_fasta_protein_homolog_model.fasta
		mv ./protein_fasta_protein_homolog_model.fasta {output.prot}
		mv ./nucleotide_fasta_protein_homolog_model.fasta {output.nucl}
		"""

rule ncbi_download:
	output:
		meta = join(DATA, 'db_metadata', 'ncbi.metadata.tsv'),
		prot = join(DATA, 'db_fastas', 'ncbi.protein.raw.fasta'),
		nucl = join(DATA, 'db_fastas', 'ncbi.nucl.raw.fasta')
	shell:
		"""
		wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt -O {output.meta}
		wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt -O {output.prot}
		wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS -O {output.nucl}
		"""

rule resfinder_download:
	params:
		outdir = join(DATA, 'resfinder_db')
	output:
		meta = join(DATA, 'db_metadata', 'resfinder.metadata.tsv'),
		nucl = join(DATA, 'db_fastas', 'resfinder.nucl.fasta')
	shell:
		"""
		wget https://bitbucket.org/genomicepidemiology/resfinder_db/raw/master/phenotypes.txt -O {output.meta}
		git clone https://bitbucket.org/genomicepidemiology/resfinder_db {params.outdir}
		cat {params.outdir}/*.fsa > {output.nucl}
		"""

rule card_parse_proteins:
	input:
		join(DATA, 'protein_fasta_protein_homolog_model.fasta')
	output:
		join(DATA, 'db_fastas', 'card.protein.fasta')
	shell:
		"""
		python scripts/parse_card_proteins.py --input {input} --output {output}
		"""

rule card_parse_metadata:
	input:
		join(DATA, 'db_metadata', 'card.metadata.obo')
	output:
		join(DATA, 'db_metadata', 'card.metadata.json')
	shell:
		"""
		python scripts/parse_card_metadata.py --input {input} --output {output}
		"""

rule ncbi_parse_proteins:
	input:
		meta = join(DATA, 'db_metadata', 'ncbi.metadata.tsv'),
		prot = join(DATA, 'db_fastas', 'ncbi.protein.raw.fasta')
	output:
		join(DATA, 'db_fastas', 'ncbi.protein.fasta')
	shell:
		"""
		python scripts/parse_ncbi_proteins.py --meta {input.meta} --input {input.prot} --output {output}
		"""

rule ncbi_parse_metdata:
	input:
		join(DATA, 'db_metadata', 'ncbi.metadata.tsv'),
	output:
		join(DATA, 'db_metadata', 'ncbi.metadata.json')
	shell:
		"""
		python scripts/parse_ncbi_metadata.py --meta {input} --output {output}
		"""

rule resfinder_parse_proteins:
	input:
		join(DATA, 'db_fastas', 'resfinder.nucl.fasta')
	output:
		join(DATA, 'db_fastas', 'resfinder.protein.fasta')
	shell:
		"""
		python scripts/resfinder_orfs_to_proteins.py --nucl {input} --protein {output}
		"""

rule resfinder_parse_metadata:
	input:
		join(DATA, 'db_metadata', 'resfinder.metadata.tsv')
	output:
		join(DATA, 'db_metadata', 'resfinder.metadata.json')
	shell:
		"""
		python scripts/parse_resfinder_metadata.py --meta {input} --output {output}
		"""

rule mmseqs_easy_rbh:
	params:
		sensitivity = 7.5
	input:
		db1 = join(DATA, 'db_fastas', '{db1}.protein.fasta'),
		db2 = join(DATA, 'db_fastas', '{db2}.protein.fasta')
	output:
		join(DATA, 'mmseqs2_results', '{db1}_vs_{db2}.rbh.tsv')
	shell:
		"""
		mmseqs easy-rbh {input.db1} {input.db2} {output} tmp -s {params.sensitivity}
		"""

rule mmseqs_createdb:
	input:
		join(DATA, 'db_fastas', '{db1}.protein.fasta'),
	output:
		join(DATA, 'mmseqs2_DBs', '{db1}.protein')
	shell:
		"""
		mmseqs createdb {input} {output}
		mmseqs createindex {output} tmp
		"""

rule mmseqs_search:
	params:
		sensitivity = 7.0,
		max_accept = 3,
		stem = join(DATA, 'mmseqs2_search_DBs', '{db1}_vs_{db2}.search')
	input:
		db1 = join(DATA, 'mmseqs2_DBs', '{db1}.protein'),
		db2 = join(DATA, 'mmseqs2_DBs', '{db2}.protein'),
	output:
		join(DATA, 'mmseqs2_search_DBs', '{db1}_vs_{db2}.search.index')
	shell:
		"""
		mmseqs search {input.db1} {input.db2} {params.stem} tmp -s {params.sensitivity} --max-accept {params.max_accept}
		"""

rule mmseqs_convertalis:
	params:
		stem = join(DATA, 'mmseqs2_search_DBs', '{db1}_vs_{db2}.search')
	input:
		db1 = join(DATA, 'mmseqs2_DBs', '{db1}.protein'),
		db2 = join(DATA, 'mmseqs2_DBs', '{db2}.protein'),
		search = join(DATA, 'mmseqs2_search_DBs', '{db1}_vs_{db2}.search.index')
	output:
		join(DATA, 'mmseqs2_results', '{db1}_vs_{db2}.search.tsv')
	shell:
		"""
		mmseqs convertalis {input.db1} {input.db2} {params.stem} {output}
		"""
