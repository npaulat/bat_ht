import sys
import os
import argparse
import itertools
import pandas as pd
import numpy as np

### Version 1.1, created 12 July 2021 by Nicole S Paulat ###

#set up access to input file w/ path
OUTDIR = os.getcwd()
INPUT = "/lustre/scratch/npaulat/RayLib-Masking/final_sp_TE_table"
TE_LIST = "/lustre/scratch/npaulat/RayLib-Masking/headers"
TYPE_LIST = "/lustre/scratch/npaulat/RayLib-Masking/header_type.txt"
SP_LIST = "/lustre/scratch/npaulat/RayLib-Masking/final_species_list.txt"
OUTPUT1 = os.path.join(OUTDIR, "final_sp_TE_heatmap.csv")
OUTPUT2 = os.path.join(OUTDIR, "final_sp_TE_heatmap_min90.csv")
OUTPUT3 = os.path.join(OUTDIR, "final_sp_TE_heatmap_min90_DNA_RC_only.csv")
OUTPUT4 = os.path.join(OUTDIR, "final_sp_TE_heatmap_min90_DNA_RC_max100sp.csv")

with open(SP_LIST) as f:
	HEADERS = list(line for line in (l.strip() for l in f) if line)

DF = pd.read_csv(INPUT, sep='\t', names=HEADERS)

with open(TE_LIST) as g:
	TES = list(line for line in (l.strip() for l in g) if line)

DF.insert(0, 'TE_names', np.array(TES))
with open(TYPE_LIST) as h:
	TYPES = list(line for line in (l.strip() for l in h) if line)

DF['TE_class'] = np.array(TYPES)
# starting table = 25676 rows (repetitive elements)
DF.to_csv(OUTPUT1, sep='\t', header=True, index=False)
DF.drop(columns=['TE_class'], inplace=True)
for column in DF.columns[1:]:
	DF.loc[DF[column] < 90, column] = 0
	DF.loc[DF[column] >= 90, column] = 1

DF['TE_class'] = np.array(TYPES)
DF.to_csv(OUTPUT2, sep='\t', header=True, index=False)
DF2 = DF[DF['TE_class'].isin(['DNA', 'RC'])].copy()
# result = 5616 rows
DF2.to_csv(OUTPUT3, sep='\t', header=True, index=False)
DF2 = DF2[(DF2 == 1).any(axis=1)]
DF2 = DF2[(DF2 == 0).any(axis=1)]
# result = 2656 rows
#col_list = list(DF2.columns[128:176])
col_list = list(DF2.columns[128:172])
DF2['bat_sum'] = DF2[col_list].sum(axis=1)
DF3 = DF2[(DF2.bat_sum > 0)].copy()
DF3.drop(columns=['bat_sum'], inplace=True)
# result = 1804 rows
col_list2 = list(DF3.columns[1:])
col_list3 = list(set(col_list2) - set(col_list))
DF3['sp_sum'] = DF3[col_list3].sum(axis=1)
# bat species = 48; other mammals = 214; conservative max sp for a more recent HT event = 100
DF4 = DF3[(DF3['sp_sum'] < 100)].copy()
# result = 1511 rows
DF4.drop(columns=['sp_sum'], inplace=True)
DF4.to_csv(OUTPUT4, sep='\t', header=True, index=False)

