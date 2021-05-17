#!/bin/bash
#SBATCH --job-name=TEmuscle
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -a 1-1013

NAMESFILE=/lustre/scratch/npaulat/yin_yang/ext_align/ht_te_list
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

cd /lustre/scratch/npaulat/yin_yang/ext_align/to_align

/lustre/work/daray/software/muscle/muscle -in ${CHRANGE}_sp_reps.fa -out ${CHRANGE}_sp_reps.aln.fa
