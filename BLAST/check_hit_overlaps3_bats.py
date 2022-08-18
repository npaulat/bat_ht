import sys
import os
import argparse
import itertools
import random
import math
from Bio import SeqIO
import pandas as pd
import numpy as np

### Version 1, created 26 April 2022 by Nicole S Paulat ###
### Modified from rm_myotis_dups_assess.py (teava) ###

### To check for overlapping/duplicate calls of a TE insertion across multiple BLAST searches
### First exclude any hits <90 bp or those <90% or >110% the TE consensus length
### Select best hit from a cluster of overlapping hits likely representing the same insertion based on:
### 1) % Sequence Identity 
### 2) E-value
### 3) Bitscore
### 4) Length
### 5) Random if absolutely identical hits
###
### Output:
### 1) Table of unique hits (including best hit from overlapping hits across TE families)
### 2) Table of overlapping hits between TE families
###
### Input format is custom BLAST 6 .out: QSEQID SSEQID PIDENT LENGTH MISMATCH GAPOPEN GAPS QSTART QEND SSTART SEND EVALUE BITSCORE

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Filtering of BLAST duplicate/overlapping calls by species", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input file of list of TEs
	#parser.add_argument('-i', '--input', help='list of TEs to check', required=True)
	#Give file of list of taxon abbreviations
	parser.add_argument('-t', '--taxon', help='taxon abbreviation input (prefix used in the BLAST .out files)', required=True)
	#Argument of directory containing the TE BLAST result subdirectories (need full path)
	parser.add_argument('-d', '--directory', type=str, help='Location of directory of the query BLAST files', default=".")
	#Argument of the output directory (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file', default=".")
	
	args = parser.parse_args()
	#TE_LIST_FILE = args.input
	SPECIES = args.taxon
	DIR = args.directory
	OUTDIR = args.outdir
	
	#argument sanity checks
	if args.taxon:
		print('The taxon of interest is ' + SPECIES +'.')
	else:
		sys.exit('You must provide a taxon.')
	
	#return TE_LIST_FILE, SPECIES, DIR, OUTDIR
	return SPECIES, DIR, OUTDIR

SPECIES, DIR, OUTDIR = get_args()

if DIR == ".":
	DIR = os.getcwd()
if OUTDIR == ".":
	OUTDIR = os.getcwd()

## Set up access to input file w/ path
#INPUT = os.path.join(DIR, TE_LIST_FILE)
#BLAST_DIR = "/lustre/scratch/npaulat/RayLib-Masking/blast_90_mamm"
BLAST_DIR = "/lustre/scratch/npaulat/RayLib-Masking/blast2_90"
#with open(SPECIES_LIST_FILE) as g:
#	SP_LIST = list(line for line in g.read().splitlines() if line)

## Make TE list object
TE_LIST_FILE = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/ht_te_417_list"
#TE_LIST_FILE = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/ht_te_225_list"
with open(TE_LIST_FILE) as g:
	TE_LIST = list(line for line in g.read().splitlines() if line)

# Make TE size dictionary
LIBRARY = "/lustre/scratch/npaulat/RayLib-Masking/te_fastas/final_mammal_library_reduced.fa"
TE_DICT = {}
TE_LIB = SeqIO.index(LIBRARY, "fasta")
for TE in TE_LIST:
	TE_LEN = len(TE_LIB[TE].seq)
	TE_DICT[TE] = TE_LEN

## Set output directory
OUTDIR = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast"
#set up access to ouput files w/ path
#BASENAME = os.path.basename(SPECIES).split("_dups")[0]
#BASE_OUT = os.path.join(OUTDIR, BASENAME)
BASE_OUT = os.path.join(OUTDIR, SPECIES)
OUTPUT = BASE_OUT + "_filtered_blast90.out"
OUTPUT2 = BASE_OUT + "_overlaps_blast90.out"


##### Function for pythonic wc -l function #####
def _make_gen(reader):
	b = reader(1024 * 1024)
	while b:
		yield b
		b = reader(1024*1024)

##### Function to check if file is empty #####
def rawgencount(filename):
	f = open(filename, 'rb')
	f_gen = _make_gen(f.raw.read)
	return sum(buf.count(b'\n') for buf in f_gen)

##### Make function which checks if hit is >=90% and <= 110% consensus len, and >= 90 bp
def hit_check(ungapped_length, cons_len):
	if ungapped_length >= math.floor(cons_len*0.9) and ungapped_length <= math.floor(cons_len*1.1):
		if ungapped_length >= 90:
			return True
		else:
			return False
	else:
		return False

##### Function to read in input file into a dictionary #####
def read_blast(SPECIES, TELIST, LEN_DICT):
	HIT_DICT = {}
	for TE in TELIST:
		FILE_NAME = BLAST_DIR + "/" + TE + "/" + TE + "_" + SPECIES + "_90.out"
		#FILE_NAME = BLAST_DIR + "/" + TE + "/" + TE + "_" + SPECIES + ".out"
		if os.path.isfile(FILE_NAME):
			NUM_HITS = rawgencount(FILE_NAME)
			if NUM_HITS > 0:
				TE_CONS = LEN_DICT[TE]
				for line in open(FILE_NAME):
					line = line.split('\t')
					QSEQID = line[0]
					SSEQID = line[1]
					PIDENT = float(line[2])
					GLENGTH = int(line[3])
					MISMATCH = int(line[4])
					GAPOPEN = int(line[5])
					GAPS = int(line[6])
					QSTART = int(line[7])
					QEND = int(line[8])
					SSTART = int(line[9])
					SEND = int(line[10])
					EVALUE = float(line[11])
					BITSCORE = int(line[12])
					LENGTH = int(GLENGTH - GAPS)
					#if the scaffold is already in the dictionary, add base pos only
					#if the scaffold is not in the dictionary, add entire entry
					if hit_check(LENGTH, TE_CONS):
						try:
							HIT_DICT[SSEQID].append((QSEQID, PIDENT, GLENGTH, MISMATCH, GAPOPEN, GAPS, QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE, LENGTH))
						except KeyError:
							HIT_DICT[SSEQID] = [(QSEQID, PIDENT, GLENGTH, MISMATCH, GAPOPEN, GAPS, QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE, LENGTH)]
		else:
			print('Warning: ' + FILE_NAME + ' does not exist.')
	return HIT_DICT

##### Make function defining what is within each hit to put on the output list #####
def blast_call(HIT_DICT,SCAFF,INSERTION1):
	return SCAFF,HIT_DICT[SCAFF][INSERTION1][0],HIT_DICT[SCAFF][INSERTION1][1],HIT_DICT[SCAFF][INSERTION1][2],HIT_DICT[SCAFF][INSERTION1][3],HIT_DICT[SCAFF][INSERTION1][4],HIT_DICT[SCAFF][INSERTION1][5],HIT_DICT[SCAFF][INSERTION1][6],HIT_DICT[SCAFF][INSERTION1][7],HIT_DICT[SCAFF][INSERTION1][8],HIT_DICT[SCAFF][INSERTION1][9],HIT_DICT[SCAFF][INSERTION1][10],HIT_DICT[SCAFF][INSERTION1][11]


# Make dictionary from input file

#HIT_DICT = read_blast("dups_test2.txt")
print('Made first dictionary.')

#for SPECIES in SP_LIST:
HIT_DICT = read_blast(SPECIES, TE_LIST, TE_DICT)
print('Made {} hit dictionary.'.format(SPECIES))
# Initialize output list (List of unique (non-overlapping) BLAST hits)
UNIQUE_LIST = []
# Initialize secondary output list (list of overlapping BLAST hits)
OVERLAP_LIST = []

# Compare lines to each other
for SCAFF in HIT_DICT:
	print('Looking at scaffold ' + SCAFF + '.')
	# Create lists per each cluster
	# Where if same scaffold and overlap in calls
	# Append to "cluster" listdups
	# Initialize list of lists of candidate clusters
	INSERTION_CANDIDATES = {x:[] for x in range(len(HIT_DICT[SCAFF]))}
	# Initialize list of insertions that are close to each other
	ALL_CLOSE = []
	# For every unique combination of calls on a given scaffold
	for INSERTION1, INSERTION2 in itertools.combinations(range(len(HIT_DICT[SCAFF])),r=2):
		#define closeness range as an insertion's ends
		a = HIT_DICT[SCAFF][INSERTION2][8]
		b = HIT_DICT[SCAFF][INSERTION2][9]
		c = HIT_DICT[SCAFF][INSERTION1][8]
		d = HIT_DICT[SCAFF][INSERTION1][9]
		# If overlapping, add to "close" list, and list of lists of candidates
		if a <= (HIT_DICT[SCAFF][INSERTION1][8] or HIT_DICT[SCAFF][INSERTION1][9]) <= b or c <= (HIT_DICT[SCAFF][INSERTION2][8] or HIT_DICT[SCAFF][INSERTION2][9]) <= d:
			INSERTION_CANDIDATES[INSERTION1].append(INSERTION2)
			ALL_CLOSE.append(INSERTION2)
	ALL_CLOSE = set(ALL_CLOSE)
	#print(INSERTION_CANDIDATES)
	# For insertion on "close" list, delete it from the list of lists of candidates (removes calls already in a cluster)
	for INSERTION in ALL_CLOSE:
		del INSERTION_CANDIDATES[INSERTION]
	# For each insertion on the list of candidate cluster lists
	for INSERTION in INSERTION_CANDIDATES:
		#print(INSERTION, INSERTION_CANDIDATES[INSERTION])
		# If insertion has no insertions clustering with it (none w/n "close" range), and not on "close" list, send to output
		if len(INSERTION_CANDIDATES[INSERTION]) == 0:
			if INSERTION not in ALL_CLOSE:
				BLAST_CALL = blast_call(HIT_DICT,SCAFF,INSERTION)
				UNIQUE_LIST.append(BLAST_CALL)
		# Otherwise, it has a cluster of "close" calls, make further filtering decisions
		else: 
			# Recall the actual line data for all insertion calls in a cluster
			INSERTION_CLUSTER = [HIT_DICT[SCAFF][INSERTION]] + [HIT_DICT[SCAFF][INSERTION2] for INSERTION2 in INSERTION_CANDIDATES[INSERTION]]
			# Double-check that there is more than one insertion on a cluster list
			if len(INSERTION_CLUSTER) > 1:
				# Make dataframe of the insertion cluster data
				df = pd.DataFrame(INSERTION_CLUSTER)
				# If any hits overlap:
				if all(x[8] <= (z[8] or z[9]) <= x[9] or z[8] <= (x[8] or x[9]) <= z[9] for x,z in itertools.combinations(INSERTION_CLUSTER,r=2)):
					# Orientation check?
					# Get hit subset with highest BITSCORE (max(x[11]))
					df_subset_MAXBIT = df[df[11] == df[11].max()].copy()
					#if more than one hit has the max BITSCORE value, make new subset with max LENGTH (max(x[2]-x[5]))
					if len(df_subset_MAXBIT) > 1:
						#df_subset_MAXLEN = df_subset_MAXBIT.copy()
						#df_subset_MAXLEN[12] = df_subset_MAXLEN[2] - df_subset_MAXLEN[5]
						#df_subset_MAXLEN2 = df_subset_MAXLEN[df_subset_MAXLEN[12] == df_subset_MAXLEN[12].max()].copy()
						df_subset_MAXLEN2 = df_subset_MAXBIT[df_subset_MAXBIT[2] == df_subset_MAXBIT[2].max()].copy()
						#if multiple hits have max BITSCORE, max LENGTH, pick max %SEQID (min(x[10])
						if len(df_subset_MAXLEN2) > 1:
							#df_subset_MAXLEN2.drop(df_subset_MAXLEN2.columns[12], axis=1, inplace=True)
							df_subset_MAXID = df_subset_MAXLEN2[df_subset_MAXLEN2[1] == df_subset_MAXLEN2[1].max()].copy()
							#if multiple hits have max BITSCORE, max LENGTH, max %SEQID, pick min EVALUE (min(x[10])
							if len(df_subset_MAXID) > 1:
								df_subset_MINEVAL = df_subset_MAXID[df_subset_MAXID[10] == df_subset_MAXID[10].min()].copy()
								#print(df_subset_MAXLEN)
								#if multiple hits have max BITSCORE, max LENGTH, max %SEQID, min EVALUE, pick one at random
								if len(df_subset_MINEVAL) > 1:
									#UNIQUE_LIST.append([SCAFF]+ list(df_subset_MINEVAL.sample(1,axis=0).iloc[0]))
									SAMPLE = df_subset_MINEVAL.sample(1,axis=0).copy()
									SAMPLE.drop(SAMPLE.columns[12], axis=1, inplace=True)
									UNIQUE_LIST.append([SCAFF]+ list(SAMPLE.iloc[0]))
								#if single hit with max BITSCORE, max LENGTH, max %SEQID, min EVALUE, append to output list
								else:
									df_subset_MINEVAL.drop(df_subset_MINEVAL.columns[12], axis=1, inplace=True)
									UNIQUE_LIST.append([SCAFF]+ list(df_subset_MINEVAL.iloc[0]))
							#if single hit with max BITSCORE, max LENGTH, and max %SEQID, append to output list
							else:
								df_subset_MAXID.drop(df_subset_MAXID.columns[12], axis=1, inplace=True)
								UNIQUE_LIST.append([SCAFF]+ list(df_subset_MAXID.iloc[0]))
						#if single hit with max BITSCORE, and max LENGTH, append to output list
						else:
							df_subset_MAXLEN2.drop(df_subset_MAXLEN2.columns[12], axis=1, inplace=True)
							UNIQUE_LIST.append([SCAFF]+ list(df_subset_MAXLEN2.iloc[0]))
					#if single hit with max BITSCORE, append to output list
					else:
						df_subset_MAXBIT.drop(df_subset_MAXBIT.columns[12], axis=1, inplace=True)
						UNIQUE_LIST.append([SCAFF]+ list(df_subset_MAXBIT.iloc[0]))
			# Append cluster of overlapping hits to OVERLAP_LIST
			for x in INSERTION_CLUSTER:
				line = list(x)
				del line[-1]
				#OVERLAP_LIST.append([SCAFF]+ list(x))
				OVERLAP_LIST.append([SCAFF]+ list(line))


# write the final output list to the output file as tab-delimited lines
UNIQUE_HITS = pd.DataFrame(UNIQUE_LIST)
UNIQUE_HITS.drop_duplicates(keep='first', inplace=True)
UNIQUE_HITS.to_csv(OUTPUT, sep="\t", header=False, index=False)

OVERLAP_HITS = pd.DataFrame(OVERLAP_LIST)
OVERLAP_HITS.drop_duplicates(keep='first', inplace=True)
OVERLAP_HITS.to_csv(OUTPUT2, sep="\t", header=False, index=False)

