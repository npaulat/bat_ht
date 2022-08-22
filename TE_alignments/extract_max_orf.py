import os
import sys
from Bio import SeqIO

te_list = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/nonmamm_cons_getorf_list"
orf_dir = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/nonmamm_cons_orfs"
out_dir = "/lustre/scratch/npaulat/RayLib-Masking/ext_align/nonmamm_cons_orfs/max_orf"

with open(te_list) as g:
	tes = list(line for line in (l.strip() for l in g) if line)

blast_list = []
for te in tes:
	basename = os.path.basename(te).split(".fa")[0]
	orf_file = os.path.join(orf_dir, basename) + "_orfs.fa"
	if os.path.isfile(orf_file) and os.path.getsize(orf_file) > 0:
		basename = os.path.basename(orf_file).split("_orfs")[0]
		out_file = os.path.join(out_dir, basename) + "_max_orf.fa"
		blast_list.append(te)
		records = list(SeqIO.parse(orf_file, "fasta"))
		records.sort(key=lambda r: -len(r))
		max_orf = records[0]
		SeqIO.write(max_orf, out_file, "fasta")

file_list = os.path.join(out_dir, "te_max_orfs_fasta_list")
with open(file_list, 'w+') as g:
	for line in blast_list:
		g.write(line + "\n")

