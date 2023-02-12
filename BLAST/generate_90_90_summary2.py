import sys
import os
import argparse
import itertools
import math
import pandas as pd
import numpy as np
from Bio import SeqIO

out_dir = "/lustre/scratch/npaulat/RayLib-Masking/external_ht"
blast_dir = "/lustre/scratch/npaulat/RayLib-Masking/external_ht/blast_out"
taxa_files = "/lustre/scratch/npaulat/RayLib-Masking/external_ht/metatheria"
summary_file = os.path.join(out_dir, "nonmammal_ht_3_90_sp_hits_summary.csv")
te_list_file = "/lustre/scratch/npaulat/RayLib-Masking/external_ht/ht_te_list"
library = "/lustre/scratch/npaulat/RayLib-Masking/te_fastas/final_mammal_library_reduced.fa"

with open(te_list_file) as g:
	te_list = list(line for line in g.read().splitlines() if line)
#add keys to te_dict with te_list
te_dict = dict.fromkeys(te_list, None)

te_dict2 = {}
#te_lib = SeqIO.to_dict(SeqIO.parse(library, "fasta")
te_lib = SeqIO.index(library, "fasta")
for te in te_list:
	te_len = len(te_lib[te].seq)
	te_dict2[te] = te_len

## Make pythonic wc -l function
def _make_gen(reader):
	b = reader(1024 * 1024)
	while b:
		yield b
		b = reader(1024*1024)

def rawgencount(filename):
	f = open(filename, 'rb')
	f_gen = _make_gen(f.raw.read)
	return sum(buf.count(b'\n') for buf in f_gen)

def hit_counter(te, blast_out, taxa_out, te_dict, te_dict2):
	species_dict = {}
	with open(taxa_out) as f:
		file_lines = list(line for line in f.read().splitlines() if line)
	
	
	for line in file_lines:
		taxon_id = line.split("|")[0].strip()
		species_id = line.split("|")[1].strip()
		taxon_genus = line.split("|")[3].strip()
		taxon_family = line.split("|")[4].strip()
		taxon_order = line.split("|")[5].strip()
		taxon_class = line.split("|")[6].strip()
		taxon_phylum = line.split("|")[7].strip()
		taxon_kingdom = line.split("|")[8].strip()
		taxon_domain = line.split("|")[9].strip()
		species_dict[taxon_id] = ((species_id, taxon_genus, taxon_family, taxon_order, taxon_class, taxon_phylum, taxon_kingdom, taxon_domain))
	hits = pd.read_csv(blast_out, sep='\t', header=None)
	hits[15] = hits[15].astype('str')
	#check if taxon id is 12 or 13
	#taxa_list = hits[12].unique().astype(str).tolist()
	taxa_list = hits[15].unique().tolist()
	
	for taxa in taxa_list:
		if taxa in species_dict:
			#hit_count = (hits[12].values == taxa).sum()
			taxa_hits = hits.loc[hits[15] == taxa]
			hit_count = len(taxa_hits)
			#hit_count_80 = hit_count - hit_count_90
			try:
				te_dict[te].append((str(hit_count),species_dict[taxa][0], species_dict[taxa][1], species_dict[taxa][2], species_dict[taxa][3], species_dict[taxa][4], species_dict[taxa][5], species_dict[taxa][6], species_dict[taxa][7], taxa))
			except AttributeError:
				te_dict[te] = [(str(hit_count),species_dict[taxa][0], species_dict[taxa][1], species_dict[taxa][2], species_dict[taxa][3], species_dict[taxa][4], species_dict[taxa][5], species_dict[taxa][6], species_dict[taxa][7], taxa)]
	
	return te_dict

#BLAST custom output 6 format: qseqid sseqid pident qcovs qcovhsp length mismatch gapopen gaps qstart qend sstart send evalue bitscore staxids sscinames scomnames
#want staxids -> field 16 (15 if index starts at 0)
#CraTho-1.145	gi|588472287|ref|NW_006533190.1|	92.030	527	39	2	38	563	85003	85527	0.0	737	176946	N/A	N/A
#.taxa files: tax_id	species_name	unique_name(?)	genus	family	order	class	phylum	kingdom	domain
#	9430	|	Desmodus rotundus	|		|	Desmodus	|	Phyllostomidae	|	Chiroptera	|	Mammalia	|	Chordata	|	Metazoa	|	Eukaryota	|

for te in te_dict:
	file = os.path.join(blast_dir, te) + "_90_90.out"
	taxa_file = os.path.join(taxa_files, te) + "_90_nonEutherian.out"
	num_hits = rawgencount(file)
	num_non_mamm = rawgencount(taxa_file)
	if num_non_mamm > 0:
		con_len = te_dict2[te]
		te_dict = hit_counter(te, file, taxa_file, te_dict, con_len)
	else:
		te_dict[te] = [("N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A")]

with open(summary_file, 'w+') as g:
	file_header = ['Element', '90_Hits', 'Name', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Domain', 'Tax_ID']
	g.write(",".join(str(x) for x in file_header) + "\n")
	for key, values in te_dict.items():
		for tuple_value in values:
			csv_row = [key] + list(tuple_value)
			g.write(",".join(str(x) for x in csv_row) + "\n")
