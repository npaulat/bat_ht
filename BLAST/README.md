## BLAST searches for all DNA transposons with patchy species distributions
  * Local genome [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) searches (blastn) of all DNA transposons annotated in 1+ bat species with 90+ copies of at least 90 bp length in the species; BLAST search criteria of at least 90% sequence identity and at least 90% query alignment
  * Local BLAST searches (blastn) for those same DNA transposons in all [NCBI eukaryote genome assemblies available](https://ftp.ncbi.nlm.nih.gov/blast/db/) (accessed April 06 2021); search criteria of at least 90% sequence identity and at least 90% query alignment
  

Local BLAST searches of genome assemblies:
  * Use final_sp_TE_heatmap_min100_DNA_RC_only.csv as input for **final_generate_blast_90_sh.py**, which will generate job submission scripts
  * Open potential_ht_3_90_sp_hits_summary.csv in Excel; use VLOOKUP of original heatmap to fill in new hit counts
    * =VLOOKUP($A3&B$2,potential_ht_hits!$A$1:$D$22855,4,0)
    * Copy worksheet to new (only values); replace "#N/A" with "_"
    * **(DB3):** List of DB2 TEs meeting presence/absence cutoff of 90 copies (20 90/90/90 hits) = ht_te_list
  * Download NCBI eukaryote genome assembly database; see **blastdb_download.sh** for databases used 
  * Make database files for each genome
    * Split job submission scripts by 1990 lines (max jobs in queue = 2000), add bash submission header to each script
    ```for f in blast_all_TEs_all_mammals.sh; do split -d -a 2 -l 1990 --additional-suffix=.sh "$f" "${f%.sh}-"; done```
    * This was done due to maximum user job limit in queue = 2000; adjust as needed for your system
    * Use jobid of blast_all_TEs_all_mammals_mkdb.sh as hold dependency for the BLAST searches
    ```for i in blast_all_TEs_all_mammals-\*; do sbatch --dependency=afterany:<jobid#> ${i}; done```
    * Command line example:
    ```
    ncbi-blast-2.11.0+/bin/makeblastdb -dbtype nucl -in <TAXON>.fa -out <TAXON>
    ```
  
  * Run local BLAST+ on eukaryote genome assembly database; see **blast_array.sh**
    * Use ht_te_list as input file (make sure to change the option -a 1-N; N=length of ht_te_list)
    * Runs BLAST+ and generate .asn and .out files for all TEs per taxon (run as array job submission)
  ```
  TAXA_FILE=</path/to/taxon/list/file>
  TE_FILE=</path/to/TE/list/file>
  TE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $TE_FILE)
  
  cd </path/to/working/directory>/blast_90_mamm/
  mkdir -p ${TE}
  cd ${TE}/
  
  echo "BLAST runs with ${TE}"
  cat $TAXA_FILE | while read line; 
   do <path/to/BLAST/installation>/ncbi-blast-2.11.0+/bin/blastn -query </path/to/TE/fasta/file/directory>/te_fastas/${TE}.fa -db </path/to/working/directory>/blast_90_mamm/${line} -perc_identity 90 -qcov_hsp_perc 90 -out </path/to/working/directory>/blast_90_mamm/${TE}/${TE}_${line}_90.asn -outfmt 11; 
   <path/to/BLAST/installation>/ncbi-blast-2.11.0+/bin/blast_formatter -archive </path/to/working/directory>/blast_90_mamm/${TE}/${TE}_${line}_90.asn -out </path/to/working/directory>/blast_90_mamm/${TE}/${TE}_${line}_90.out -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore"; 
  done
  echo "${TE} BLAST custom output created"
  ```
  
  * Create raw summary (potential_ht_3_90_sp_hits_summary.csv) file of BLAST hits for all TEs and taxa with **make_blast_90_summary.py**
    * Run generate_ht_summary_90_90.py to determine species with any 90/90/90 hits in non-mammal species separately (ht_te_list)
  * Filter BLAST results to resolve duplicate calls with best match from candidate HTT library (417 TEs); see **check_overlaps.py**
  * Reformat filtered BLAST results for downstream use in the orthologous TE searches and summarization; see **reformat_filtered_blast.py**
  * Make summary file of filtered BLAST results for all HTTs and all mammals; see **mk_mamm_filter_blast_summary.py**
