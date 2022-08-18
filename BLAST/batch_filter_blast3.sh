#!/bin/bash 
#SBATCH --job-name=check_blast
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --mem-per-cpu=5150MB
#SBATCH --nodes=2
#SBATCH --ntasks=37
#SBATCH -a 1-37

NAMESFILE=/lustre/scratch/npaulat/RayLib-Masking/filter_blast/species_abbrev_list
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

cd /lustre/scratch/npaulat/RayLib-Masking/filter_blast/

python /lustre/scratch/npaulat/RayLib-Masking/filter_blast/check_hit_overlaps3_bats.py -t ${CHRANGE} -d . -od .

