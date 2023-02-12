#!/bin/bash
#SBATCH --job-name=TEblast
#SBATCH --output=%x.%j.out
#SBATCH --partition=nocona
#SBATCH --mem-per-cpu=5150MB
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH -a 1-417

###This array is for BLAST searches within bat genomes, BLAST output format is modified from that of blast_array.sh (excludes "qcovs qcovhsp staxids sscinames scomnames")

NAMESFILE=/lustre/scratch/npaulat/zoonomia_RayLib-Masking/blast2_90/ht_te_list
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)
NAMESFILE2=/lustre/scratch/npaulat/zoonomia_RayLib-Masking/blast2_90/bat_list

cd /lustre/scratch/npaulat/zoonomia_RayLib-Masking/blast2_90/

mkdir -p ${CHRANGE}
cd ${CHRANGE}/

echo "BLAST runs with ${CHRANGE}"
cat $NAMESFILE2 | while read line; do /lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blastn -query /lustre/scratch/npaulat/zoonomia_RayLib-Masking/te_fastas/${CHRANGE}.fa -db /lustre/scratch/npaulat/zoonomia_RayLib-Masking/blast2_90/${line} -perc_identity 90 -qcov_hsp_perc 90 -out /lustre/scratch/npaulat/zoonomia_RayLib-Masking/blast2_90/${CHRANGE}/${CHRANGE}_${line}_90.asn -outfmt 11; /lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blast_formatter -archive /lustre/scratch/npaulat/zoonomia_RayLib-Masking/blast2_90/${CHRANGE}/${CHRANGE}_${line}_90.asn -out /lustre/scratch/npaulat/zoonomia_RayLib-Masking/blast2_90/${CHRANGE}/${CHRANGE}_${line}_90.out -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore"; done

echo "${CHRANGE} BLAST custom output created"
