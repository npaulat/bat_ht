#!/bin/bash
#SBATCH --job-name=TEraxml
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -a 1-183

NAMESFILE=/lustre/scratch/npaulat/RayLib-Masking/ext_align/to_align/te_multi_sp_aln_files
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

#ALN_FILE=$1
#NAME=$(basename $ALN_FILE -all_aln.fa)
#/lustre/work/aosmansk/apps/standard-RAxML/raxmlHPC-PTHREADS -T 2 -f a -m GTRGAMMA -p 12345 -x 98765 -# 1000 -s $ALN_FILE -n trees/$NAME"-rBoot.tre"

cd /lustre/scratch/npaulat/RayLib-Masking/ext_align/trees

/lustre/work/aosmansk/apps/standard-RAxML/raxmlHPC-PTHREADS -T 1 -f a -m GTRGAMMA -p 12345 -x 98765 -# 1000 -s $CHRANGE -n $(basename $CHRANGE _sp_reps.aln.fa)"-rBoot.tre"
