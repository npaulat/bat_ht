import sys
import os
import argparse
import itertools
import random
import pandas as pd
import numpy as np
import pathlib

### Version 1, created 2 May 2022 by Nicole S Paulat ###

### To separate filtered blast result files of each species by TE, slightly reformatted by swapping column 1 and 2's order, and then creating BED files for each of those
###Input format is custom BLAST 6 .out: SSEQID	QSEQID	PIDENT LENGTH MISMATCH GAPOPEN GAPS QSTART QEND SSTART SEND EVALUE BITSCORE
### Output:
###		1. More standard custom BLAST 6 .out for each taxon, each TE: QSEQID	SSEQID	PIDENT LENGTH MISMATCH GAPOPEN GAPS QSTART QEND SSTART SEND EVALUE BITSCORE
###		2. From output 1 format, BED files with SSEQID	SSTART SEND

## separate filtered_blast files by TE and swap column 1 and 2 positions
## need this file for making extracts for calculating TE div/age
## then take those files and make separate reduced files of just TE name, scaffold, start, stop for the ortholog comparisons

## Make TE list object
TE_LIST_FILE = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/ht_te_221_list"
with open(TE_LIST_FILE) as g:
	TE_LIST = list(line for line in g.read().splitlines() if line)


TE_TYPE_LIST_FILE = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/te_type_list"
with open(TE_TYPE_LIST_FILE) as g:
	TYPE_LIST = list(line for line in g.read().splitlines() if line)


## Make species list object
#TAXON_LIST_FILE = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/taxon_abbrev_list"
#TAXON_LIST_FILE = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast/taxon_abbrev_list3"
#with open(TAXON_LIST_FILE) as g:
#	SPECIES_LIST = list(line for line in g.read().splitlines() if line)
SPECIES_LIST = ['eFus2', 'aPal2']


## Set input file directory
INPUT_DIR = "/lustre/scratch/npaulat/RayLib-Masking/filter_blast"

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

##### Function to read in input file into TE-separated output files #####
def read_blast(SPECIES, TELIST, TYPELIST, BLAST_DIR):
	HEADERS=['SSEQID', 'QSEQID', 'PIDENT', 'LENGTH', 'MISMATCH', 'GAPOPEN', 'GAPS', 'QSTART', 'QEND', 'SSTART', 'SEND', 'EVALUE', 'BITSCORE']
	OUTDIR1 = BLAST_DIR + "/filtered_beds/"
	OUTFILE1 = OUTDIR1 + "/" + SPECIES + "_blast90_flank50.bed"
	FILE = BLAST_DIR + "/" + SPECIES  + "_filtered_blast90.out"
	if os.path.isfile(FILE):
		NUM_HITS = rawgencount(FILE)
		if NUM_HITS > 0:
			DATA = pd.read_csv(FILE, sep="\t", names=HEADERS)
			DATA.drop(columns=['PIDENT', 'LENGTH', 'MISMATCH', 'GAPOPEN', 'GAPS', 'QSTART', 'QEND', 'EVALUE', 'BITSCORE'], inplace=True)
			MASK = DATA['SEND'] < DATA['SSTART']
			DATA.loc[MASK, ['SSTART','SEND']] = DATA.loc[MASK, ['SEND','SSTART']].values
			DATA.sort_values(by=['SSEQID', 'SSTART'], ascending=[False, True], inplace=True)
			DATA['NAME'] = DATA['SSEQID'].map(str) + ':' + DATA['SSTART'].map(str) + '-' + DATA['SEND'].map(str)
			DATA['TE'] = 'x'
			for TE in TELIST:
				SUB = TE + "#"
				TYPE = next((i for i in TYPELIST if SUB in i), None)
				DATA.loc[DATA['QSEQID'] == TE, 'TE'] = TYPE
			DATA2 = DATA.loc[DATA['TE'] != 'x'].copy()
			DATA2.drop(columns=['QSEQID'], inplace=True)
			DATA2['SSTART'] = DATA['SSTART'] - 50
			DATA2['SEND'] = DATA['SEND'] + 50
			DATA2.to_csv(OUTFILE1, sep="\t", header=False, index=False)
	for TE in TELIST:
		SUB = TE + "#"
		TYPE = next((i for i in TYPELIST if SUB in i), None)
		OUTDIR2 = BLAST_DIR + "/filtered_blast90/" + TE + "/"
		pathlib.Path(OUTDIR2).mkdir(exist_ok=True)
		OUTFILE2 = OUTDIR2 + "/" + SPECIES + "_" + TE + "_blast90.out"
		if os.path.isfile(FILE):
			NUM_HITS = rawgencount(FILE)
			if NUM_HITS > 0:
				DATA = pd.read_csv(FILE, sep="\t", names=HEADERS)
				DATA2 = DATA[DATA['QSEQID'] == TE].copy()
				if len(DATA2) > 0:
					DATA2.insert(0, 'QSEQID', DATA2.pop('QSEQID'))
					DATA2.drop(columns=['GAPS'], inplace=True)
					DATA2.sort_values(by=['SSEQID', 'SSTART'], ascending=[False, True], inplace=True)
					DATA2.to_csv(OUTFILE2, sep="\t", header=False, index=False)

OUTDIR1 = INPUT_DIR + "/filtered_beds/"
pathlib.Path(OUTDIR1).mkdir(exist_ok=True)
OUTDIR2 = INPUT_DIR + "/filtered_blast90/"
pathlib.Path(OUTDIR2).mkdir(exist_ok=True)
for SPECIES in SPECIES_LIST:
	read_blast(SPECIES, TE_LIST, TYPE_LIST, INPUT_DIR)

