from Bio import SeqIO
import os

out_dir = "/lustre/scratch/npaulat/yin_yang/ext_align/to_align"
cons_dir = "/lustre/scratch/npaulat/yin_yang/ext_align/species_consensus_seqs"
te_list_file = "/lustre/scratch/npaulat/yin_yang/ext_align/ht_te_list"
with open(te_list_file) as g:
	te_list = list(line for line in g.read().splitlines() if line)

for te in te_list:
	subdirectory = os.path.join(cons_dir, te)
	files = [os.path.join(subdirectory, file) for file in os.listdir(subdirectory)]
	output = os.path.join(out_dir, te) + "_sp_reps.fa"
	with open(output, "w") as out:
		for file in files:
			species = os.path.basename(file).split("_")[1]
			with open(file, "r") as i_file:
				seq_record = SeqIO.read(i_file, 'fasta')
				seq_id = seq_record.id + "_" + species
				seq_record.id = seq_id
				seq_record.description = ""
				SeqIO.write(seq_record, out, 'fasta')
