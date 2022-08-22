import os
import sys
from Bio import SeqIO

#add species suffix to sequence headers with Biopython
#concatenate seq files for each TE for MUSCLE alignment
#remove blank lines

out_dir = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/to_align"
cons_dir = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/species_consensus_seqs"
te_list_file = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/ht_te_list"
with open(te_list_file) as g:
	te_list = list(line for line in g.read().splitlines() if line)

te_list2 = []
te_list3 = []
output4 = os.path.join(out_dir, "tes_multi_sp_reps_list")
output5 = os.path.join(out_dir, "tes_single_sp_rep_list")
output6 = os.path.join(out_dir, "te_multi_sp_aln_files")

file_list = []
for te in te_list:
	subdirectory = os.path.join(cons_dir, te)
	files = [os.path.join(subdirectory, file) for file in os.listdir(subdirectory)]
	output = os.path.join(out_dir, te) + "_sp_reps.fa"
	file_list.append(output)
	if len(files) > 2:
		te_list2.append(te)
	else:
		te_list3.append(te)
	with open(output, "w+") as out:
		seq_records = []
		for file in files:
			species = os.path.basename(file).split("_")[-2]
			with open(file, "r") as i_file:
				seq_record = SeqIO.read(i_file, 'fasta')
				seq_id = seq_record.id + "_" + species
				seq_record.id = seq_id
				seq_record.description = ""
				seq_records.append(seq_record)
		SeqIO.write(seq_records, out, 'fasta')

with open(output4, "w+") as f:
	for line in te_list2:
		f.write(line + "\n")

with open(output5, "w+") as f:
	for line in te_list3:
		f.write(line + "\n")

with open(output6, "w+") as f:
	for line in te_list2:
		f.write(out_dir + "/" + line + "_sp_reps.aln.fa\n")

#make bash script to run
output2 = os.path.join(out_dir, "te_rep_seqs_list")
with open(output2, 'w+') as g:
	for line in file_list:
		g.write(line + "\n")

output3 = os.path.join(out_dir, "te_reps_muscle_aln.sh")
with open(output3, 'w+') as f:
	f.write("#!/bin/bash\n")
	f.write("#SBATCH --job-name=TEmuscle\n")
	f.write("#SBATCH --output=%x.%j.out\n")
	f.write("#SBATCH --error=%x.%j.err\n")
	f.write("#SBATCH --partition=nocona\n")
	f.write("#SBATCH --nodes=1\n")
	f.write("#SBATCH --ntasks=1\n")
	f.write("#SBATCH -a 1-{}\n\n".format(str(len(te_list2))))
	f.write("NAMESFILE={}\n".format(output4))
	f.write('CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)\n\n')
	f.write("cd {}\n\n".format(out_dir))
	f.write("/lustre/work/daray/software/muscle/muscle -in ${CHRANGE}_sp_reps.fa -out ${CHRANGE}_sp_reps.aln.fa\n")
