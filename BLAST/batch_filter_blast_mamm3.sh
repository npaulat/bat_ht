#!/bin/bash 
#SBATCH --job-name=check_blast
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --mem-per-cpu=5150MB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -a 1-214

NAMESFILE=/lustre/scratch/npaulat/RayLib-Masking/filter_blast/mamm_list
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

cd /lustre/scratch/npaulat/RayLib-Masking/filter_blast/

python /lustre/scratch/npaulat/RayLib-Masking/filter_blast/check_hit_overlaps3.py -t ${CHRANGE} -d . -od .
