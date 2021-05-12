import sys
import os
import argparse
import itertools
import math
from Bio import SeqIO

blast_dir = "/lustre/scratch/npaulat/yin_yang/blast_80"
summary_file = os.path.join(blast_dir, "potential_ht_3_80_sp_hits_summary.csv")
library = "/lustre/scratch/npaulat/yin_yang/te_fastas/final_mammal_library_reduced.fa"
## Make TE list w/ TE consensus length
te_list_file = "/lustre/scratch/npaulat/yin_yang/potential_ht_tes_list_reduced"
te_dict = {}
with open(te_list_file) as g:
	te_list = list(line for line in g.read().splitlines() if line)

#te_lib = SeqIO.to_dict(SeqIO.parse(library, "fasta")
te_lib = SeqIO.index(library, "fasta")
for te in te_list:
	te_len = len(te_lib[te].seq)
	te_dict[te] = te_len

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

def hit_counter(blast_out, te_dict, te):
	good_hit_count = 0
	with open(blast_out) as f:
		file_lines = list(line for line in f.read().splitlines() if line)
	for line in file_lines:
		#if default BLAST fmt 6:
		#hit_len = int(line.split("\t")[3]) - int(line.split("\t")[5])
		hit_len = int(line.split("\t")[3]) - int(line.split("\t")[6])
		if hit_len >= math.floor(te_dict[te]*0.8):
			if hit_len >= 80:
				good_hit_count += 1
	return good_hit_count

#BLAST output 6 format: query_id, subject_id, perc_identity, alignment_length, mismatches, gaps, gap_opens, q_start, q_end, s_start, s_end, e_value, bit_score
#AmnHarb1	KN195036.1	81.092	238	41	4	4	714	948	76701	76937	1.18e-45	191
#Actual hit length will be alignment_length - gaps (# gaps)
#/lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blast_formatter -archive /lustre/scratch/npaulat/yin_yang/blast_80/AmnHarb1/AmnHarb1_TolMat.asn -out /lustre/scratch/npaulat/yin_yang/blast_80/AmnHarb1/AmnHarb1_TolMat.out -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore"

#need to exlude 4 bat1k species b/c not published (gSor) or no closely related species with curated de novo TE library
exclude_list = ["gSor", "sBil", "nTum", "hMon"]

summary_list = {}
for te in te_dict:
	#for each file in subdirectory:
	#this assumes no subdirs have been added to TE blast subdirs from blast runs
	#assumes only BLAST .out files in each subdirectory
	subdirectory = os.path.join(blast_dir, te)
	files = [os.path.join(subdirectory, file) for file in os.listdir(subdirectory) if file.endswith('.out')]
	for file in files:
		species_name = os.path.basename(file).split("_")[-1].split(".out")[0]
		if species_name not in exclude_list:
			num_hits = rawgencount(file)
			if num_hits > 0:
				good_hit_count = hit_counter(file, te_dict, te)
			else:
				good_hit_count = 0
			try:
				summary_list[te].append((species_name, good_hit_count))
			except KeyError:
				summary_list[te] = [(species_name, good_hit_count)]

with open(summary_file, 'a+') as g:
	for key, values in summary_list.items():
		for tuple_value in values:
			csv_row = [key] + list(tuple_value)
			g.write(",".join(str(x) for x in csv_row) + "\n")
