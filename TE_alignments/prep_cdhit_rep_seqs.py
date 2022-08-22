import os
import sys
from Bio import SeqIO

#concatenate seq files for each TE for CD-HIT_EST alignment
#removes blank lines

cons_dir = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/to_align"
te_list_file = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/ht_te_list"
lib_fastas = "/lustre/scratch/npaulat/RayLib-Masking/te_fastas"
out_dir = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/cdhit"
with open(te_list_file) as g:
	te_list = list(line for line in g.read().splitlines() if line)

file_list = []
for te in te_list:
	sp_cons_file = os.path.join(cons_dir, te) + "_sp_reps.fa"
	#sp_cons_file = sp_cons + "_sp_reps.fa"
	lib_cons_file = os.path.join(lib_fastas, te) + ".fa"
	#lib_cons_file = lib_cons + ".fa"
	files = [lib_cons_file, sp_cons_file]
	output = os.path.join(out_dir, te) + "_all_reps.fa"
	file_list.append(output)
	with open(output, "w+") as out:
		seq_records = []
		for file in files:
			with open(file, "r") as i_file:
				for record in SeqIO.parse(i_file, "fasta"):
					seq_records.append(record)
		SeqIO.write(seq_records, out, 'fasta')

#make bash script to run
output2 = os.path.join(out_dir, "te_all_rep_seqs_list")
with open(output2, 'w+') as g:
	for line in file_list:
		g.write(line + "\n")

output3 = os.path.join(out_dir, "te_reps_cdhit.sh")
with open(output3, 'w+') as f:
	f.write("#!/bin/bash\n")
	f.write("#SBATCH --job-name=TEcdhit\n")
	f.write("#SBATCH --output=%x.%j.out\n")
	f.write("#SBATCH --error=%x.%j.err\n")
	f.write("#SBATCH --partition=nocona\n")
	f.write("#SBATCH --nodes=1\n")
	f.write("#SBATCH --ntasks=1\n")
	f.write("#SBATCH -a 1-{}\n\n".format(str(len(te_list))))
	f.write("NAMESFILE={}\n".format(te_list_file))
	f.write('CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)\n\n')
	f.write("cd {}\n\n".format(out_dir))
	f.write("/lustre/work/daray/software/cdhit-4.8.1/cd-hit-est -i ${CHRANGE}_all_reps.fa -o ${CHRANGE}_all_reps_cdhit90_aL09 -c 0.9 -aL 0.9 -n 5 -M 2200\n")
