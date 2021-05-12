#!/bin/bash 
#SBATCH --job-name=te-count
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -a 1-253



NAMESFILE=/lustre/scratch/npaulat/yin_yang/te-counts/LIST
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

cd /lustre/scratch/npaulat/yin_yang/te-counts

for i in $(cat headers); do echo $CHRANGE $i $(grep $'\t'$i$'\t' /lustre/scratch/npaulat/yin_yang/bed_higher_score/$CHRANGE"_rm.bed" | wc -l) >> /lustre/scratch/npaulat/yin_yang/te-counts/$CHRANGE"_HT"; done
