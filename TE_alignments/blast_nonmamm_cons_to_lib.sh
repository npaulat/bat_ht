#!/bin/bash 
#SBATCH --job-name=blast2lib
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1

cd /lustre/scratch/npaulat/RayLib-Masking/ext_align/blast_to_lib

#/lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/makeblastdb -dbtype nucl -in /lustre/scratch/npaulat/RayLib-Masking/te_fastas/final_mammal_library_reduced.fa -out mammlib

for i in /lustre/scratch/npaulat/RayLib-Masking/ext_align/nonmamm_cons_seqs/*.fa; do 
/lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blastn -query ${i} -db mammlib -perc_identity 90 -qcov_hsp_perc 90 -out $(basename ${i} .fa)_2lib.asn -outfmt 11; /lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blast_formatter -archive $(basename ${i} .fa)_2lib.asn -out $(basename ${i} .fa)_2lib.out -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore"; done

