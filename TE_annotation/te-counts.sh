#!/bin/bash 
#SBATCH --job-name=te-count
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -a 1-258

NAMESFILE=/lustre/scratch/npaulat/RayLib-Masking/species_abbrev.txt
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

cd /lustre/scratch/npaulat/RayLib-Masking/te-counts

for i in $(cat headers); do echo $CHRANGE $i $(grep $'\t'$i$'\t' /lustre/scratch/npaulat/RayLib-Masking/rm2bed_files/$CHRANGE"_rm.bed" | wc -l) >> /lustre/scratch/npaulat/RayLib-Masking/te-counts/$CHRANGE"_HT"; done
