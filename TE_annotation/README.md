### Annotating TEs for a) TE content summary and b) identifying TEs with patchy species distributions
1. Create species list in order of appearance on your phylogeny (i.e. in order of relatedness within each clade, in this case, mammalian order, family, and genus)
2. Annotate TEs in each species with RepeatMasker (v4.1.0)
3. Reformat and filter RepeatMasker output with [**RM2Bed.py**](https://github.com/davidaray/bioinfo_tools/blob/master/RM2bed.py) with -min 90 and -o higher_score options
    * TE annotations can be summarized visually using the **scaledviolinplot.py** and **barplot_panel.py** scripts in the visualization folder
4. Get TE counts with **te-counts.sh**
5. Get third column from these TE count files
    * ```for i in *_HT; do awk -F ' ' '{print $3}' ${i} > ${i}2; done```
6. Paste together these HT2 files in order of the species phylogeny using **paste_final_table.sh**; copy final_table to working directory as final_sp_TE_table
7. Run **final_generate_heatmap_tables.py** in order to reformat and filter TE data for DNA/RC elements with limited distributions involving bats (relevant output is final_sp_TE_heatmap_min100_DNA_RC_only.csv)
    * Use this table as input for BLAST search workflow
