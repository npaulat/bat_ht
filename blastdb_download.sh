#!/bin/bash
#SBATCH --job-name=blastdb
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -a 1-69


NAMESFILE=/lustre/scratch/aosmansk/blastdb/DB_LIST
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

#Navigate to directory

cd /lustre/scratch/aosmansk/blastdb

#Download blastdb

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_euk_rep_genomes.$CHRANGE.tar.gz

#Decompress the downloaded directory

tar -zxvf ref_euk_rep_genomes.$CHRANGE.tar.gz
