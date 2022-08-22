#!/bin/bash
#SBATCH --job-name=TEmuscle
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -a 1-183

NAMESFILE=/lustre/scratch/npaulat/RayLib-Masking/ext_align/to_align/tes_multi_sp_reps_list
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

cd /lustre/scratch/npaulat/RayLib-Masking/ext_align/to_align

/lustre/work/daray/software/muscle/muscle -in ${CHRANGE}_sp_reps.fa -out ${CHRANGE}_sp_reps.aln.fa
