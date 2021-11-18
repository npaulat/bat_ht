import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
plt.switch_backend('Agg')
from pylab import savefig
import re
import argparse
import numpy as np


#Typical usage: python barplot_panel.py -i1 plots/all_all_taxa_classes_trimmed_species_updated3.txt -i2 plots/50my_all_taxa_DNA_families_trimmed_species_updated3.txt -l y -o h

####MAIN function
def main():

##Use the get_args function
	INPUT1, INPUT2, LEGEND, ORIENT, DATATYPE, BOXWIDTH, SUBPLOTHORIZONAL, SUBPLOTVERTICAL, SUBPLOTWIDTH, SUBPLOTHEIGHT, COLUMNS, SPACE = get_args()
	print('First input file = ' + INPUT1)
	print('Second input file = ' + INPUT2)

	PREFIX = re.split("[.]", INPUT1)[0]
	
	INPUTFRAME1 = pd.read_table(INPUT1, sep = '\t', index_col=0)
	with open(INPUT1, 'r') as f:
		line = f.readline().strip()
	NAMES1 = line.split("\t")
	NAMES1 = [x.replace("_", " ") for x in NAMES1]
	INPUTFRAME1 = INPUTFRAME1.transpose()
	
	INPUTFRAME2 = pd.read_table(INPUT2, sep = '\t', index_col=0)
	with open(INPUT2, 'r') as f:
		line = f.readline().strip()
	NAMES2 = line.split("\t")
	NAMES2 = [x.replace("_", " ") for x in NAMES2]
	INPUTFRAME2 = INPUTFRAME2.transpose()
	
	if NAMES1 != NAMES2:
		sys.exit("Names for the two inputs are not the same or are not in the same order. Please resolve.")
	
	#Set font sizes
	LABEL_SIZE = 20
	AXES_SIZE = 20
	#DIMS = (14.5, 15.5)
	DIMS = (20, 19)
	#Set custom color cycler (blue, orange, green, pink, purple, brown)
	#mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c', '#17becf', '#9467bd', '#A9561E'])
	mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#1f77b4', '#ff7f0e', '#9467bd', '#17becf', '#2ca02c', '#A9561E'])
	FIG, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=DIMS)
	FIG.subplots_adjust(hspace=0.1)
	
	if LEGEND == 'n' and ORIENT == 'v':
		FIG = INPUTFRAME1.plot.bar(stacked=True, width=SPACE, figsize=DIMS, legend=False)
		labels = range(len(NAMES1))
		plt.xticks(labels, NAMES1, size=AXES_SIZE, style='italic', rotation=60, ha="right")
		plt.yticks(size=AXES_SIZE)
		#plt.xlabel('Taxon', size=LABEL_SIZE)
		plt.ylabel('Genome Proportion', size=LABEL_SIZE)
		plt.tick_params(width=2)
		#Make sure the figure is adjusted to fit the axes labels into the frame
		plt.tight_layout()
		#FIG.figure.savefig(PREFIX + '_vertical_nolegend.png', dpi=300)
		#FIG.figure.savefig(PREFIX + '_vertical_nolegend.svg', format='svg')
	elif LEGEND == 'n' and ORIENT == 'h':
		FIG = INPUTFRAME1.plot.barh(stacked=True, width=SPACE, figsize=DIMS, legend=False)
		labels = range(len(NAMES1))
		plt.yticks(labels, NAMES1, size=AXES_SIZE, style='italic')
		plt.xticks(size=AXES_SIZE)
		#plt.ylabel('Taxon', size=LABEL_SIZE)
		plt.xlabel('Genome Proportion', size=LABEL_SIZE)
		plt.tick_params(width=2)
		#Make sure the figure is adjusted to fit the axes labels into the frame
		plt.tight_layout()
		#FIG.figure.savefig(PREFIX + '_horizontal_nolegend.png', dpi=300)
		#FIG.figure.savefig(PREFIX + '_horizontal_nolegend.svg', format='svg')
	elif LEGEND == 'y' and ORIENT == 'v':
		#plt.subplot(111)
		FIG = INPUTFRAME.plot.bar(stacked=True, width=SPACE, figsize=DIMS, legend=False)
		labels = range(len(NAMES1))
		plt.xticks(labels, NAMES1, size=AXES_SIZE, style='italic', rotation=60, ha="right")
		plt.yticks(size=AXES_SIZE)
		#plt.xlabel('Taxon', size=LABEL_SIZE)
		plt.ylabel('Genome Proportion', size=LABEL_SIZE)
		plt.tick_params(width=2)
		#lgd = plt.legend(loc=2, fontsize=AXES_SIZE, bbox_to_anchor=(1.01, 1), ncol=COLUMNS, borderaxespad=0.)
		lgd = plt.legend(loc=1, fontsize=AXES_SIZE, bbox_to_anchor=(0.99, 1), ncol=COLUMNS, markerscale=0.8, frameon=False, framealpha=1.0, edgecolor=None, borderaxespad=0.)
		#Make sure the figure is adjusted to fit the axes labels into the frame
		#plt.tight_layout()
		FIG.figure.savefig(PREFIX + '_standard_vertical_plot' + '.png', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
		FIG.figure.savefig(PREFIX + '_standard_vertical_plot' + '.svg', bbox_extra_artists=(lgd,), bbox_inches='tight', format='svg')
	else:
		#plt.subplot(111)
		#FIG = INPUTFRAME1.plot.barh(stacked=True, width=SPACE, figsize=DIMS, legend=False)
		#ax1 = INPUTFRAME1.plot.barh(stacked=True, width=SPACE, figsize=DIMS, legend=False)
		#INPUTFRAME1.plot.barh(ax=ax1, stacked=True, width=SPACE, figsize=DIMS, legend=False)
		ax1 = INPUTFRAME1.plot.barh(ax=ax1, stacked=True, width=SPACE, legend=False)
		labels = range(len(NAMES1))
		#plt.yticks(labels, NAMES1, size=AXES_SIZE, style='italic')
		ax1.yaxis.set_tick_params(labelsize=AXES_SIZE)
		#plt.xticks(size=AXES_SIZE)
		ax1.tick_params(width=2, labelsize=AXES_SIZE)
		#plt.ylabel('Taxon', size=LABEL_SIZE)
		plt.xlabel('Genome Proportion', size=LABEL_SIZE)
		#plt.tick_params(width=2)
		#legend positioned at upper right corner in plot; about 1/3 of the way down
		#ax1.legend(loc=1, fontsize=20, bbox_to_anchor=(0.985, 0.88), ncol=COLUMNS, markerscale=0.8, labelspacing = 0.3, handlelength = 1.8, handletextpad = 0.6, frameon=True, framealpha=1.0, edgecolor=None, borderaxespad=0.)
		ax1.legend(loc=1, fontsize=18, bbox_to_anchor=(0.985, 0.88), ncol=COLUMNS, markerscale=0.6, labelspacing = 0.3, handlelength = 1.7, handletextpad = 0.6, frameon=True, framealpha=1.0, edgecolor=None, borderaxespad=0.)
		ax1.set_title('All Transposable Elements', fontsize=LABEL_SIZE)
		
		#Generate subplot 2
		#ax2 = INPUTFRAME2.plot.barh(stacked=True, width=SPACE, figsize=DIMS, legend=False)
		#INPUTFRAME2.plot.barh(ax=ax2, stacked=True, width=SPACE, figsize=DIMS, legend=False)
		ax2 = INPUTFRAME2.plot.barh(ax=ax2, stacked=True, width=SPACE, legend=False)
		labels = range(len(NAMES2))
		plt.yticks(labels, NAMES1, size=AXES_SIZE, style='italic')
		plt.xticks(size=AXES_SIZE)
		plt.tick_params(width=2)
		#Make legend for subplot2
		lgd2 = plt.legend(loc=1, fontsize=18, bbox_to_anchor=(0.985, 0.88), ncol=COLUMNS, markerscale=0.7, labelspacing = 0.3, handlelength = 1.7, handletextpad = 0.6, frameon=True, framealpha=1.0, edgecolor=None, borderaxespad=0.)
		ax2.set_title(r'DNA Transposons $\leq$50 My', fontsize=LABEL_SIZE)
		
		#fking fix italics formatting
		for i in range(0,47):
			ax1.get_yticklabels()[i].set_fontstyle('italic')
		
		#Make sure the figure is adjusted to fit the axes labels into the frame
		FIG.tight_layout()
		FIG.subplots_adjust(hspace=0)
		FIG.tight_layout()
		#FIG.figure.savefig(PREFIX + '_h_plot_panel' + '.png', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
		#FIG.figure.savefig(PREFIX + '_h_plot_panel' + '.svg', bbox_extra_artists=(lgd,), bbox_inches='tight', format='svg')
		FIG.savefig(PREFIX + '_h_plot_panel' + '.png', dpi=300, bbox_inches='tight')
		FIG.savefig(PREFIX + '_h_plot_panel' + '.svg', bbox_inches='tight', format='svg')

##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Will process a saved pandas dataframe of data from RM2bed.py and processed through any of several catdata python scripts", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i1', '--input1', type=str, help='Name of first input file to be parsed', required=True)
	parser.add_argument('-i2', '--input2', type=str, help='Name of second input file to be parsed', required=True)
	parser.add_argument('-l', '--legend', type=str, help='n = do not include legend. Default = y.', default='y', required=True)
	parser.add_argument('-o', '--orientation', type=str, help='Vertical bars (v) or horizontal (h).', default='v')
	parser.add_argument('-d', '--datatype', type=str, help='CLASS, FAMILY, SINE, LINE, LTR, etc., depending on file')
	parser.add_argument('-b', '--boxwidth', type=float, help='Fraction of the standard size plot to allow for plotting the legend outside of that box. 0.8 (80%) will usually work')
	parser.add_argument('-ho', '--subplothorizontal', type=float, help='Number to indicate where to horizontally place the legend. Default is 1.01, just to the right of the main plot.', default = 1.01)
	parser.add_argument('-v', '--subplotvertical', type=float, help='Number to indicate where to vertically place the legend. Default is 1, just at the top of the main plot.', default = 1)
	parser.add_argument('-w', '--subplotwidth', type=float, help='Number to indicate how wide the legend should be. Default is .1, enough for one column.', default = .1)
	parser.add_argument('-sh', '--subplotheight', type=float, help='Number to indicate how tall the legend should be. Usually done automatically. Default = .1', default = 1)
	parser.add_argument('-c', '--numcol', type=int, help='Number of columns in the legend. Default = 1', default = 1)
	parser.add_argument('-s', '--spacing', type=float, help='Bar spacing. 1 = no space. 0 = no bar. Default = 0.9.', default = 0.9)

	args = parser.parse_args()
	INPUT1 = args.input1
	INPUT2 = args.input2
	LEGEND = args.legend
	ORIENT = args.orientation
	DATATYPE = args.datatype
	BOXWIDTH = args.boxwidth
	SUBPLOTHORIZONAL = args.subplothorizontal
	SUBPLOTVERTICAL = args.subplotvertical
	SUBPLOTWIDTH = args.subplotwidth
	SUBPLOTHEIGHT = args.subplotheight
	COLUMNS = args.numcol
	SPACE = args.spacing

	return INPUT1, INPUT2, LEGEND, ORIENT, DATATYPE, BOXWIDTH, SUBPLOTHORIZONAL, SUBPLOTVERTICAL, SUBPLOTWIDTH, SUBPLOTHEIGHT, COLUMNS, SPACE

if __name__ =="__main__":main()	
	
