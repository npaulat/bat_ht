import sys
import os
import argparse
import itertools
import math
from Bio import SeqIO

blast_dir = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast"
#summary_file = os.path.join(blast_dir, "mamm_ht_417_90_filtered_hits_summary_check.csv")
#summary_file = os.path.join(blast_dir, "mamm_ht_225_90_filtered_hits_summary_check.csv")
summary_file = os.path.join(blast_dir, "mamm_ht_227_90_filtered_hits_summary_check2.csv")
#summary_file = os.path.join(blast_dir, "bat_ht_417_90_filtered_hits_summary_check.csv")
#summary_file = os.path.join(blast_dir, "bat_ht_225_90_filtered_hits_summary.csv")
library = "/lustre/scratch/npaulat/RayLib-Masking/te_fastas/final_mammal_library_reduced.fa"
## Make TE list w/ TE consensus length
#te_list_file = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/ht_te_417_list"
#te_list_file = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/ht_te_225_list"
te_list_file = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/ht_te_227_list"
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

# def hit_counter(blast_out, te_dict, te):
	# good_hit_count = 0
	# with open(blast_out) as f:
		# file_lines = list(line for line in f.read().splitlines() if line)
	# for line in file_lines:
		# if line.split("\t")[1] == te:
		# #if default BLAST fmt 6:
		# #hit_len = int(line.split("\t")[3]) - int(line.split("\t")[5])
			# hit_len = int(line.split("\t")[3]) - int(line.split("\t")[6])
			# if hit_len >= math.floor(te_dict[te]*0.9) and hit_len <= math.floor(te_dict[te]*1.1):
				# if hit_len >= 90:
					# good_hit_count += 1
	# return good_hit_count

def hit_counter(blast_out, te):
	good_hit_count = 0
	with open(blast_out) as f:
		file_lines = list(line for line in f.read().splitlines() if line)
	for line in file_lines:
		if line.split("\t")[1] == te:
			good_hit_count += 1
	return good_hit_count


#BLAST output 6 format: sseqid qseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore
#scaffold_m16_p_14	Mariner2_mMol	93.297	2715	147	22	35	10	2699	53737582	53734878	0.0	3975

summary_list = {}
files = [os.path.join(blast_dir, file) for file in os.listdir(blast_dir) if file.endswith('filtered_blast90.out')]
for te in te_dict:
	#subdirectory = os.path.join(blast_dir, te)
	#files = [os.path.join(subdirectory, file) for file in os.listdir(subdirectory) if file.endswith('filtered_blast90.out')]
	for file in files:
		species_name = os.path.basename(file).split("_")[0]
		num_hits = rawgencount(file)
		if num_hits > 0:
			#good_hit_count = hit_counter(file, te_dict, te)
			good_hit_count = hit_counter(file,  te)
		else:
			good_hit_count = 0
		try:
			summary_list[te].append((species_name, good_hit_count))
		except KeyError:
			summary_list[te] = [(species_name, good_hit_count)]
		#summary_list.append(species_name + "\t" + numhits)
	#sum_file = os.path.join(subdirectory, "summary_file.txt")

with open(summary_file, 'a+') as g:
	for key, values in summary_list.items():
		for tuple_value in values:
			csv_row = [key] + list(tuple_value)
			g.write(",".join(str(x) for x in csv_row) + "\n")


