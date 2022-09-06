## Visualizing TE content
  * This requires the reformatted RepeatMasker output from RM2Bed.py from the TE annotation workflow as the initial input
  * The TE_rm.bed files are then processed using the [filter_beds.py](https://github.com/davidaray/bioinfo_tools/blob/master/filter_beds.py) script, once with the age filter of 50my and once without
    * This requires an additional input, genome_sizes_mrates.txt, which is a list file with 1) taxon abbreviations, 2) size of the genome in bp, and 3) the species-specific neutral mutation rate
    * See bat_genome_sizes_mrates.txt
  ```
  python filter_beds.py -g genome_sizes_mrates.txt -p 50my -a 50000000
  python filter_beds.py -g genome_sizes_mrates.txt
  ```
    * This assumes the *_rm.bed files are in the current directory, see catdata.sh for convenient setup
  * The original TE_rm.bed files must also be processed using the [catdata_props_age.py](https://github.com/davidaray/bioinfo_tools/blob/master/catdata_props_age.py) script, for creating processed files split by TE class, TE family, and TE subfamily, as well as merged files with summaries of TE proportions
  ```
  python catdata_props_age.py -g genome_sizes_mrates.txt -p all
  ```

### Stacked barplots of TE content
  * This requires you to take the output file all_all_taxa_classes_merged_cats.txt from catdata_props_age.py and trim it down to the desired TE classes (can use Excel)
    * In this case, LINE, SINE, RC, DNA, and Unknown >> all_all_taxa_classes_trimmed.txt
  * For the barplot of DNA/RC TEs 50My or younger, repeat the catdata_props_age.py processing step with the filtered 50My output from filter_beds.py (or set an age limit in python to get filtered TEs)
  * This time, will need to concatenate the desired families files (<TAXON>_<TE_Family>_50my_processed_beds.txt) for each TE superfamily
    * In this case, Helitron, hAT, piggyBac, TcMariner, and other DNA are the superfamilies
  *
  ```
  python catdata_props_age.py -g genome_sizes_mrates.txt -p 50my
