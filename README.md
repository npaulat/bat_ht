## Identification of putative HT DNA transposons involving bats

This workflow is to computationally identify putative horizontally transferred (HT) DNA transposons, done here by searching across 251 mammalian species, specifically focused on those involved in the 37 bat species.\
\
Given the wide variety of DNA transposons and species involved, this is a broad-scale search, with _a priori_ search thresholds. The main steps are:
1. Transposable element (TE) annotation using [RepeatMasker](http://repeatmasker.org/) in 253 mammalian genome assemblies
2. Local genome [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) searches (blastn) of all DNA transposons annotated in 1+ bat species with 100+ copies of at least 90 bp length in the species; search criteria of at least 90% sequence identity and at least 90% query alignment
3. Local BLAST searches (blastn) for those same DNA transposons in all [NCBI eukaryote genome assemblies available](https://ftp.ncbi.nlm.nih.gov/blast/db/) (accessed April 06 2021); search criteria of at least 90% sequence identity and at least 90% query alignment
4. Generate species-specific consensus sequences for all transposons; for non-mammalian species, must have 20+ hits
5. Use [CD-HIT-EST](http://weizhongli-lab.org/cd-hit/) to confirm all species-specific consensus sequences meet the 90/90/90 rule for the original library element (at least 90% sequence identity, 90% sequence length, 90 bp length); exclude any that do not meet this criteria
6. Identify autonomous elements: 
     * Use [EMBOSS getorf](https://www.bioinformatics.nl/cgi-bin/emboss/help/getorf) utility to identify open-reading frames (ORFs) for all species-specific consensus sequences greater than 800 bp 
     * Perform blastx searches of non-redundant proteins for the largest ORFs from each consensus sequence (can be performed remotely or in [web browser](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome))
     * Once autonomous elements are identified, exclude any with less than 20 copies meeting the 90/90/90 criteria; exclude any non-autonomous elements with less than 100 copies
7. Generate [RAxML](https://github.com/stamatak/standard-RAxML) trees of species-specific consensus sequences for each element
     * Only useful if present in 3 or more species; used to search for phylogenetic incongruence with species tree. Since this project focused on putative HT specifically involving bats, and many elements were not found outside a single clade, this step was not particularly informative.
8. Identify and exclude deletion products of elements within the set of putative HT elements by clustering sequences via a) the cross_match utility of [Phrap](http://www.phrap.org/phredphrapconsed.html) v0.990319 with default settings, and b) a modified CD-HIT search with the utility [ClusterPartialMatchingSubs.pl](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpz1.154) using default settings.
    * Further identification of deletion products by manually comparing aligned sequences.
9. Estimate average age of a TE subfamily in a given species using the average modified Kimura 2-parameter (K2P) distance from the library consensus and the neutral mutation rate of the species. Modified K2P values calculated in [RepeatMasker](http://repeatmasker.org/)'s utilities.
10. Infer putative HT placement on bat phylogeny (see bat_phylogeny) based on presence/absence and average TE age across species.
    * Bat phylogeny based on Foley et al. 2021 and [Amador et al. 2018](https://link.springer.com/article/10.1007/s10914-016-9363-8) and used a combination of non-conflicting average or median divergnence estimates from [TimeTree](http://www.timetree.org/) (accessed 3 September 2021).
    * To be conservative, elements assigned to oldest possible branch based on presence/absence data within a given clade (i.e. there are four _Myotis_ species representing to sister species pairs (for this tree), _M. brandtii + lucifugus_ and _M. davidii + myotis_; if only _M. brandtii_ and _myotis_ were searched, and the element was found in both, it was assumed to also be in the other two species, and so the HT event would be inferred to have occurred in the ancestral _Myotis_ lineage.)
11. Estimate association between young (>50 My) TE accumulation in bats and species richess; association between putative HT diversity and species richness in bats.


### Workflow
1. Create species list in order of appearance on your phylogeny
2. Annotate TEs in each species with RepeatMasker
3. RM2Bed conversion with -min 80 and -o higher_score
4. Get TE counts with te-counts.sh
5. Get third column from these TE count files
    * for i in \*\_HT; do awk -F ' ' '{print $3}' ${i} > ${i}2; done
6. Paste together these HT2 files in order of the species phylogeny using paste_final_table.sh; copy final_table to working directory as final_sp_TE_table
7. Run generate_heatmap_tables.py
8. Use final_sp_TE_heatmap_min100_DNA_RC_only.csv as input for generate_blast_80_sh.py
9. Split job submission scripts by 1990 lines (max jobs in queue = 2000), add bash submission header to each script
    for f in blast_all_TEs_all_mammals.sh; do split -d -a 2 -l 1990 --additional-suffix=.sh "$f" "${f%.sh}-"; done
    * This was done due to maximum user job limit in queue = 2000; adjust as needed for your system
10. Run BLAST submission scripts
    * Use jobid of blast_all_TEs_all_mammals_mkdb.sh as hold dependency for the BLAST searches
      * for i in blast_all_TEs_all_mammals-\*; do sbatch --dependency=afterany:<jobid#> ${i}; done
11. Run make_blast_summary.py script
12. Open potential_ht_3_80_sp_hits_summary.csv in Excel; use VLOOKUP of original heatmap to fill in new hit counts
    * =VLOOKUP($A3&B$2,potential_ht_hits!$A$1:$D$22855,4,0)
    * Copy worksheet to new (only values); replace "#N/A" with "_"
    * **(DB3):** List of DB2 TEs meeting presence/absence cutoff of 100 copies (20 90/90/90 hits) = ht_te_list (old version was TE_LIST2)
13. Use ht_te_list as input file in blast_array.sh (make sure to change -a 1-N; N=length of ht_te_list), submit blast_array.sh to run blastn locally on eukaryote genome assemblies (see blastdb_download.sh script for databases used)
14. Run generate_ht_summary_90_90.py to determine species with 20+ 90/90/90 hits (ht_te_list)
    * Run extend_align to get species-specific consensus sequences for any species with 20+ hits
15. Use ht_te_list as input for generate_bash_nonmamm_cons.sh to run array jobs of ext_align scripts on non-mammals
    * Use EMBOSS getorf to find ORFs, then do a blastx search of nonredundant proteins to search for 90% ID, 90% length match to a transposase or such
    * If autonomous (based on ORF + blastx), include all species with 20+ hits; if non-autonomous, include all species with 100+ hits
16. Use potential_ht_3_90_sp_hits_summary.csv and ht_te_list as input for generate_bash_mamm_cons.py to generate array jobs of ext_align scripts on mammals
17. Copy all consensus sequences into species_consensus_seqs/ subdirectory, then move into subdirectories by TE name
     * For mammals:
         * for SUBDIR in \*/; do cp ${SUBDIR::-1}/final_consensuses/${SUBDIR::-1}\_rep.fa /lustre/scratch/npaulat/yin_yang/ext_align/species_consensus_seqs/; done
         * OR for i in */; do cp ${i}/final_consensuses/* ../sp_consensus_seqs/; done
    
    * For non-mammals:
         * LIST=\["list of all ht TEs"]
    
         * cd /lustre/scratch/npaulat/yin_yang/ext_align/nonmamm_species_cons_seqs
         * for i in \*.fa; do cp ${i} /lustre/scratch/npaulat/yin_yang/ext_align/species_consensus_seqs/; done
         * cd ../species_consensus_seqs/
         * for ITEM in $LIST: do mkdir ${ITEM}; mv "${ITEM}\_"* ${ITEM}/ 2>/dev/null; done
18. Run concatenate_te_sp_rep_seqs.py to group species' consensus sequences by element to form a multi line FASTA for each element, and generate files for te_reps_muscle_aln.sh.
19. Align TE consensus sequences with MUSCLE by running te_reps_muscle_aln.sh
20. For TEs with 3+ species sequences (tes_multi_sp_aln_files), generate RAxML trees.
    * 
