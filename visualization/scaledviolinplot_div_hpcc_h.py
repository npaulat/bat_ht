# Violin Plots
# Jenna R. Grimshaw
# Oct 11, 2019
# Modified by Jenny Korstian 4/3/2020
# Modified by David Ray 5/21/2021

import pandas as pd
import os
import re
import glob
import numpy as np
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#plt.switch_backend('Agg')
from pylab import savefig

# 1. Loop through all rm bed files
# 2. Concatenate all into a single file
# 3. Create violin plots

#######
## This script assumes you want to generate a violin plot from the 'class' files generated
## by cat_data_props.py followed by filtering via filter_beds.py.
#######

##### Housekeeping #####
## Decide on your clade
CLADES = ['chiroptera']
#CLADES = ['chiroptera', 'sciuridae', 'cetartiodactyla', 'afrotheria', 'ruminantia', 'canoidea-pinnepedia', 'feloidea-plus', 'hystricomorpha', 'sciuromorpha', 'myomorpha-plus', 'primates']

## Decide if you want to do this on TE class or TE family
CATEGORY = 'class'
#CATEGORY = 'family'

## Update paths here if necessary
PROCESSEDBEDSPATH = '/lustre/scratch/npaulat/paulat_beds2/' + CATEGORY + '_files'
WORKPATH = '/lustre/scratch/npaulat/paulat_beds2/plots'

## Decide what class of TE you want to plot
#TETYPES = ['DNA'] #for testing
TETYPES = ['DNA', 'LINE', 'RC', 'SINE', 'LTR']
#TETYPES = ['Helitron', 'TcMariner', 'hAT', 'piggyBac']
## Decide whether you want 'all' rows of bed files or some subcategory filtered by filter_beds.py
RANGE = 'all'
#RANGE = '50my'

##Decide on a divergence or age value to use
#DIVERGENCE = 50
AGE = 50000000

##See sizefile input in first loop below

## Check file names. As written this limits to only DNAs and 50my files
for CLADE in CLADES:
	## Provide a file that has a list of the names of each species, one per line. This should be available as the size file provided to 'filter_beds.py'.
	SIZEFILE = WORKPATH + '/' + CLADE + '_sizefile.txt'
	#SIZEFILE = WORKPATH + '/' + 'genome_sizes_mrates.txt'
	# Read in the genomesizes to get a taxon list.
	print('Reading in genome sizes file, ' + SIZEFILE + '.')
	#GENOMESIZES = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu'], squeeze=True, index_col=0)
	GENOMESIZES = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu', 'taxon_abbrev'], squeeze=True, index_col=0)
	#GENOMESIZESFRAME = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu'], squeeze=True)
	GENOMESIZESFRAME = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu', 'taxon_abbrev'], squeeze=True)
	TAXA = GENOMESIZESFRAME['taxon'].tolist()
	TAXA.reverse()
	IDS = GENOMESIZESFRAME['taxon_abbrev'].tolist()
	IDS.reverse()
	TAXALENGTH = len(TAXA)
	print('Will plot ages and accumulation for ' + str(TAXALENGTH) + ' species.')
	NAMES = [x.replace("_", " ") for x in TAXA]

	for TETYPE in TETYPES:
		# Initialize a dataframe
		ALLTAXA = pd.DataFrame()
		
		#for TAXON in TAXA:
		for TAXON, ID in zip(TAXA, IDS):
			# Read in the filename for a given taxon
			print('Reading in ' + ID + '_' + TETYPE + '_' + RANGE + '_' + CATEGORY +'_processed_beds.txt')
			FILE = pd.read_csv(PROCESSEDBEDSPATH + '/' + ID + '_' + TETYPE + '_' + RANGE + '_' + CATEGORY +'_processed_beds.txt', sep="\t")
			# Rename the columns to make downstream processing a bit easier
			#FILE = FILE.rename(columns={TAXON + '_TE' : 'TE', TAXON + '_size' : 'Size', TAXON + '_class' : 'Class', TAXON + '_family' : 'Family', TAXON + '_div' : 'Div', TAXON + '_age' : 'Age'})
			FILE = FILE.rename(columns={ID + '_TE' : 'TE', ID + '_size' : 'Size', ID + '_class' : 'Class', ID + '_family' : 'Family', ID + '_div' : 'Div', ID + '_age' : 'Age'})
			# Add column with species name
			FILE['Taxon'] = TAXON
			# Filter to only keep elements that are more than 100 bp in length
			FILE=FILE[FILE.Size >= 100]
			# Filter to only keep TEs that occur at least 100 times
			FILE = (FILE.groupby("TE").filter(lambda x : len(x)>100))
			# Create new dataframe with only taxon id and age
			#FILE = FILE[['Taxon', 'Div']]
			FILE = FILE[['Taxon', 'Age']]
			print(FILE.head())
			# Append to growing dataframe
			ALLTAXA = pd.concat([ALLTAXA, FILE], ignore_index=True)
		#ALLTAXA.to_csv(WORKPATH + '/'  + CLADE + '_' + TETYPE + 's_div' + str(DIVERGENCE) + '_violinframe.txt', index=False)
		# Optional: filter to just young elements (ex. 10my)
		#ALLTAXA=ALLTAXA[ALLTAXA.Div <= DIVERGENCE]
		ALLTAXA=ALLTAXA[ALLTAXA.Age <= AGE]
		#DIMS = (20, 15)
		#DIMS = (18, 16)
		#DIMS = (9.5, 9)
		DIMS = (14.25, 13.25)
		FIG, ax = plt.subplots(figsize=DIMS)
		#Set font sizes
		LABEL_SIZE = 19
		AXES_SIZE = 19
		if TETYPE == 'other_DNA':
			TITLE = TETYPE.strip().replace('_', ' ')
			TITLE = re.sub("(^|\s)(\S)", convert_to_uppercase, TITLE)
		else:
			TITLE = TETYPE
		##Make horizontal violin plot
		#sns.violinplot('Age', 'Taxon', data=ALLTAXA, scale="count", ax=ax, order = IDS)
		#sns.violinplot('Age', 'Taxon', data=ALLTAXA, scale="count", ax=ax, order = TAXA, orient='h')
		sns.violinplot('Age', 'Taxon', data=ALLTAXA, scale="count", ax=ax, order = TAXA, orient='h', cut=0)
		plt.tick_params(width=1.5)
		ax.set_yticklabels(NAMES, size=AXES_SIZE, style='italic')
		#plt.xticks(rotation=60, ha="right")
		x_ticks = np.arange(0, 50000001, 10000000)
		plt.xticks(x_ticks)
		plt.xticks(fontsize=AXES_SIZE)
		plt.ticklabel_format(axis='x', style='sci', scilimits=(6,6), useMathText=True)
		#ax.xaxis.offsetText.set_fontsize(AXES_SIZE)
		plt.setp(ax.get_xaxis().get_offset_text(), visible=False)
		ax.set_title(TITLE + 's', size=LABEL_SIZE)
		#Set divergence limits
		#plt.ylim(0, 50)
		#Set age limits
		#plt.xlim(-9000000, 59000000)
		#plt.xlim(-1000000, 51000000)
		plt.xlim(0, 50000000)
		#ax.xaxis.grid(True)
		#ax.set_ylabel('Taxon', size=LABEL_SIZE)
		y_axis = ax.axes.get_yaxis()
		#y_axis.set_label_text('y')
		y_axis.label.set_visible(False)
		#ax.set_ylabel('Divergence')
		ax.set_xlabel('Age (million years)', size=LABEL_SIZE)
		plt.tight_layout()
		#FIG.savefig(WORKPATH + '/'  + CLADE + '_' + TETYPE + 's_div' + str(DIVERGENCE) + 'violin' + '.png')
		#FIG.savefig(WORKPATH + '/'  + CLADE + '_' + TETYPE + 's_div' + str(DIVERGENCE) + 'violin' + '.svg')
		FIG.savefig(WORKPATH + '/'  + CLADE + '_' + TETYPE + 's_age' + str(AGE) + 'violin_horizontal' + '.png', dpi=300)
		FIG.savefig(WORKPATH + '/'  + CLADE + '_' + TETYPE + 's_age' + str(AGE) + 'violin_horizontal' + '.svg')
		
