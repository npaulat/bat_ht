#!/bin/bash
#SBATCH --job-name=TEorfs
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -a 1-7

NAMESFILE=/lustre/scratch/npaulat/RayLib-Masking/ext_align/nonmamm_cons_getorf_list
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

cd /lustre/scratch/npaulat/RayLib-Masking/ext_align/

/lustre/work/daray/software/EMBOSS-6.6.0/emboss/getorf -sequence /lustre/scratch/npaulat/RayLib-Masking/ext_align/nonmamm_cons_seqs/${CHRANGE} -outseq /lustre/scratch/npaulat/RayLib-Masking/ext_align/nonmamm_cons_orfs/$(basename $CHRANGE .fa)_orfs.fa -find 3
