## Visualizing TE content
  * This requires the reformatted RepeatMasker output from RM2Bed.py from the TE annotation workflow as the initial input
  * The TE_rm.bed files are then processed using the [**filter_beds.py**](https://github.com/davidaray/bioinfo_tools/blob/master/filter_beds.py) script, once with the age filter of 50my and once without
    * This requires an additional input, genome_sizes_mrates.txt, which is a list file with 1) taxon abbreviations, 2) size of the genome in bp, and 3) the species-specific neutral mutation rate
    * See bat_genome_sizes_mrates.txt
    * This assumes the \*\_rm.bed files are in the current directory, see catdata.sh for convenient setup
  ```
  python filter_beds.py -g genome_sizes_mrates.txt -p 50my -a 50000000
  python filter_beds.py -g genome_sizes_mrates.txt
  ```
  
  * The original TE_rm.bed files must also be processed using the [**catdata_props_age.py**](https://github.com/davidaray/bioinfo_tools/blob/master/catdata_props_age.py) script, for creating processed files split by TE class, TE family, and TE subfamily, as well as merged files with summaries of TE proportions
  ```
  python catdata_props_age.py -g genome_sizes_mrates.txt -p all
  ```

### Stacked barplots of TE content
  * This requires you to take the output file all_all_taxa_classes_merged_cats.txt from catdata_props_age.py and trim it down to the desired TE classes (can use Excel)
    * In this case, LINE, SINE, RC, DNA, and Unknown >> all_all_taxa_classes_trimmed.txt
  * For the barplot of DNA/RC TEs 50My or younger, will need to concatenate the desired families files (\<TAXON>_\<TE_Family>_all_processed_beds.txt) for each TE superfamily
    * In this case, Helitron, hAT, piggyBac, TcMariner, and other DNA are the superfamilies
  * Can do this by creating list files with the TE families for each superfamily, and using the taxon list file in some bash loops
  ```
  for TAXON in $(cat list.txt); do for TE in hAT Helitron piggyBac TcMariner other_DNA; do while read line; cat ${TAXON}_${line}_all_family_processed_beds.txt >> ${TAXON}_${TE}_all_family_processed_beds.txt; done < ${TE}_family_list; done

for TAXON in $(cat list.txt); do for TE in hAT Helitron piggyBac TcMariner other_DNA; do header=${TAXON}_TE; sed -i "/^${header}/d" ${TAXON}_${TE}_all_family_processed_beds.txt; done < ${TE}_family_list; done
  ```
  * Again, create a trimmed file with the total proportions for each category (can either sum in Excel, or in python)
  * Finally, use these two trimmed TE proportion summary files as input for **barplot_panel.py**
  * Example:
  ``` 
  python barplot_panel.py -i1 plots/all_all_taxa_classes_trimmed_species_updated4.txt -i2 plots/50my_all_taxa_DNA_families_trimmed_species_updated4.txt -l y -o h
  ```
 #### Boxplot inset for total DNA/RC TE content
  * Summarize TE content by age categories using **calc_te_proportions.py** (or manually in Excel or python)
  ```
  python calc_te_proportions.py
  ```
  * Open the resulting bat_te_proportions_summary file, and create a sheet named "Total DNA_RC" with the summed total DNA + RC TE proportions for each taxon
  * This sheet should have three columns: 1) Genus_species, 2) Group (Chiroptera or Others), and 3) Total_Proportion
  * Run **make_boxplot.R** in RStudio, tweak figure as desired
 
 
 ### Violin plots for temporal TE accumulation patterns
  * This requires the <TAXON>_<TECLASS>_50my_family_processed_beds.txt files made in the previous steps, which should all be in a single directory
  * Also requires a <CLADE>_sizefile.txt, which lists 1) Genus_species, 2) genome size in bp, 3) species-specific neutral mutation rate, and 4) taxon abbreviation (GenSpe)
    * See **scaledviolinplot_div_hpcc_family_h.py** annotations for details on paths and figure options
  ```
  python scaledviolinplot_div_hpcc_family_h.py
  ```

### Phylogenetic tree with HT events
 * Requires a phylogeny in Nexus format (.nex), and a FAD_LAD file (see bat_fad_lad.tsv)
  * Can make a basic phylogeny in Mesquite and export as a .nex file
  * FAD_LAD file contains the genus, FAD (estimated time of divergence in My), and LAD (in this case, 0.000000001 since the R package will not accept 0)
 * Run **bat_geoscale.R**, which adds the geologic timescale to your phylogeny
   * The main issue with this package is that the current version overwrites the base option to change the figure size, so it produces a very small figure with weird proportions; can fix by 1) making the plot viewer as large as possible, 2) save figure, 3) import into Inkscape and convert to an svg, and 4) manually change the lines and labels to an appropriate width/length and font size for readability
 * Also in Inkscape, choose a series of shapes and color shading to represent your TE types and heatmapping of HT events, and add them to the appropriate branches
 * Create a figure legend
 * Export as as appropriate file format (.svg, .png) with desired dpi (300, 600)
