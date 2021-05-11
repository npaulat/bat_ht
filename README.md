# bat_ht_project

## Workflow
1. Create species list in order of appearance on your phylogeny
2. Annotate TEs in each species with RepeatMasker
3. RM2Bed conversion with -min 80 and -o higher_score
4. Get TE counts with te-counts.sh
5. Get third column from these TE count files
      * awk command
6. Paste together these HT2 files in order of the species phylogeny
7. Run generate_heatmap_tables.py
8. Use final_sp_TE_heatmap_min100_DNA_RC_only.csv as input for 
9. Split job submission scripts by 1990 lines (max jobs in queue = 2000)
    for f in blast_all_TEs_all_mammals.sh; do split -d -a 2 -l 5000 --additional-suffix=.sh "$f" "${f%.sh}-"; done
10. Run BLAST submission scripts
    * Use jobid of blast_all_TEs_all_mammals_mkdb.sh as hold dependency for the BLAST searches
      * for i in blast_all_TEs_all_mammals-\*; do sbatch --dependency=afterany:<jobid#> ${i}; done
11. Run make_blast_summary.py script
12. Open potential_ht_3_80_sp_hits_summary.csv in Excel; use VLOOKUP of original heatmap to fill in new hit counts
    * =VLOOKUP($A3&B$2,potential_ht_hits!$A$1:$D$22855,4,0)
    * Copy worksheet to new (only values); replace "#N/A" with "_"
    * **(DB3):** List of DB2 TEs meeting presence/absence cutoff of 20 copies (20 80/80/80 hits) = ht_te_list
13. Use ht_te_list as input file in blast_array.sh (make sure to change -a 1-N; N=length of ht_te_list), submit blast_array.sh to run blastn locally on eukaryote genome assemblies
14. Run generate_ht_summary_90_90.py to determine species with 20+ 90/90/90 hits (ext_ht_te_list)
15. Use ext_ht_te_list as input for generateto run array jobs of ext_align scripts on non-mammals
16. Use potential_ht_3_80_sp_hits_summary.csv and potential_ht_tes_reduced as input for generate_bash_mamm_cons.py to generate array jobs of ext_align scripts on mammals
17. Copy all consensus sequences into species_consensus_seqs/ subdirectory, then move into subdirectories by TE name
     * **For mammals**
         * for SUBDIR in \*/; do cp ${SUBDIR::-1}/final_consensuses/${SUBDIR::-1}\_rep.fa /lustre/scratch/npaulat/yin_yang/ext_align/species_consensus_seqs/; done
    
    * **For non-mammals**
         * LIST=\["list of all ht TEs"]
    
         * cd /lustre/scratch/npaulat/yin_yang/ext_align/nonmamm_species_cons_seqs
         * for i in \*.fa; do cp ${i} /lustre/scratch/npaulat/yin_yang/ext_align/species_consensus_seqs/; done
         * cd ../species_consensus_seqs/
         * for ITEM in $LIST: do mkdir ${ITEM}; mv "${ITEM}\_"* ${ITEM}/ 2>/dev/null; done
