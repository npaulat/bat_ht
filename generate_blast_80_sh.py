import sys
import os
import argparse
import itertools
import pandas as pd
import numpy as np

OUTDIR = "/lustre/scratch/npaulat/yin_yang/blast_80/"
INPUT1 = "/lustre/scratch/npaulat/yin_yang/final_sp_TE_heatmap_min100_DNA_RC_max100sp.csv"
INPUT2 = "/lustre/scratch/npaulat/yin_yang/heatmap_new_headers"
OUTPUT1 = os.path.join(OUTDIR, "blast_all_TEs_all_mammals_mkdb.sh")
OUTPUT2 = os.path.join(OUTDIR, "blast_all_TEs_all_mammals.sh")
GENOME_DIR = "/lustre/scratch/npaulat/yin_yang/masked_assemblies/"
TE_DIR = "/lustre/scratch/npaulat/yin_yang/te_fastas/"

with open(INPUT2) as g:
	NEW_HEADERS = list(line for line in (l.strip() for l in g) if line)

DF = pd.read_csv(INPUT1, sep="\t", header=0)
DF.columns = NEW_HEADERS

TE_DICT = {}
for row in range(len(DF)):
	for column in range(263):
		if DF.iloc[row, column] == 1:
			try:
				TE_DICT[DF.iloc[row, 0]].append(DF.columns[column])
			except KeyError:
				TE_DICT[DF.iloc[row, 0]] = [DF.columns[column]]

#DF.drop(columns=['TE_class'], inplace=True)
SP_LIST = []
DF.loc['col_total'] = DF.sum(0, numeric_only=True)
total_index = len(DF) - 1
for column in range(1, 263):
	if DF.iloc[total_index, column] > 0:
		SP_LIST.append(DF.columns[column])

#BLAST output 6 format: query_id, subject_id, perc_identity, alignment_length, mismatches, gap_opens, q_start, q_end, s_start, s_end, e_value, bit_score
#Also want gaps
#-outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore"

with open(OUTPUT1, "w+") as file1:
	with open(OUTPUT2, "w+") as file2:
		for SPECIES in SP_LIST:
			GENOME = GENOME_DIR + SPECIES + ".fa"
			file1.write("/lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/makeblastdb -dbtype nucl -in " + GENOME + " -out " + SPECIES + "\n")
		for TE in TE_DICT:
			TE_SUBDIR = os.path.join(OUTDIR, TE)
			if not os.path.isdir(TE_SUBDIR):
				os.makedirs(TE_SUBDIR)
			for SPECIES in TE_DICT[TE]:
				TE_FA = TE_DIR + TE + ".fa"
				OUT_FILE1 = TE_SUBDIR + "/" + TE + "_" + SPECIES + ".asn"
				OUT_FILE2 = TE_SUBDIR + "/" + TE + "_" + SPECIES + ".out"
				file2.write('/lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blastn -query ' + TE_FA + ' -db ' + SPECIES + ' -perc_identity 80 -qcov_hsp_perc 80 -out ' + OUT_FILE1 + ' -outfmt 11\n')
				file2.write('/lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blast_formatter -archive ' + OUT_FILE1 + ' -out ' + OUT_FILE2 + ' -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore"\n')
