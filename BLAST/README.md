Local BLAST searches of genome assemblies:

  * Make database files for each genome
  ```
  ncbi-blast-2.11.0+/bin/makeblastdb -dbtype nucl -in <TAXON>.fa -out <TAXON>
  ```
  * Run BLAST+ and generate .asn and .out files for all TEs per taxon (run as array job submission)
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
  
  * Create summary file of BLAST hits for all TEs and taxa with **make_blast_90_summary.py**
  
  * Download NCBI eukaryote genome assembly database; see **blastdb_download.sh** 
  * Run BLAST+ on eukaryote genome assembly database; see **blast_array.sh**
  * Filter BLAST results to resolve duplicate calls with best match from candidate HTT library (417 TEs); see **check_overlaps.py**
  * Reformat filtered BLAST results for downstream use in the orthologous TE searches and summarization; see **reformat_filtered_blast.py**
  * Make summary file of filtered BLAST results for all HTTs and all mammals; see **mk_mamm_filter_blast_summary.py**
