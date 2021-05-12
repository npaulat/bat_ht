import sys
import os

blast_dir = "/lustre/scratch/npaulat/yin_yang/blast_80"
summary_file = os.path.join(blast_dir, "potential_ht_3_80_sp_hits_summary.csv")
#te_list_file = "/lustre/scratch/npaulat/yin_yang/TE_LIST2"
ht_te_file = "/lustre/scratch/npaulat/yin_yang/ext_align/ht_te_list"
output_file = "/lustre/scratch/npaulat/yin_yang/ext_align/mammal_blast_ht_te_list"
output_file2 = "/lustre/scratch/npaulat/yin_yang/ext_align/mamm/batch_mamm_generate_cons.sh"

te_dict = {}
#with open(te_list_file) as g:
#	te_list = list(line for line in g.read().splitlines() if line)
#
#for te in te_list:
#	te_dict[te] = "None"

with open(summary_file) as g:
	file_lines = list(line for line in g.read().splitlines() if line)

bat_sp_list = ['MacSob', 'EidHel', 'RouAeg', 'PteVam', 'PteAle', 'RhiFer', 'rCfl', 'rTri', 'rSed', 'rAff', 'RhiSin', 'hCyc', 'aSto', 'HipGal', 'hLar', 'HipArm', 'MegLyr', 'CraTho', 'MurFea', 'MyoLuc', 'MyoBra', 'MyoDav', 'MyoMyo', 'PipPip', 'PipKuh', 'aPal', 'LasBor', 'eFus', 'MinSch', 'MinNat', 'MolMol', 'tBra', 'nLep', 'PtePar', 'MorBla', 'aJam', 'sHon', 'CarPer', 'AnoCau', 'pHas', 'PhyDis', 'TonSau', 'DesRot', 'MicHir']

for line in file_lines:
	te = line.split(",")[0].strip()
	species_id = line.split(",")[1].strip()
	hit_count = int(line.split(",")[2].strip())
	if species_id in bat_sp_list:
		if hit_count >= 100:
				te_dict[te] = "None"

with open(ht_te_file, "w+") as f:
	for key in te_dict:
		f.write(key + "\n")

for line in file_lines:
	te = line.split(",")[0].strip()
	species_id = line.split(",")[1].strip()
	hit_count = int(line.split(",")[2].strip())
	if te in te_dict:
		if hit_count >= 100:
			try:
				te_dict[te].append((species_id, str(hit_count)))
			except AttributeError:
				te_dict[te] = [(species_id, str(hit_count))]

with open(output_file, 'a+') as g:
	for key, values in te_dict.items():
		for tuple_value in values:
			file_row = [key] + list(tuple_value)
			g.write("\t".join(str(x) for x in file_row) + "\n")

with open(output_file, "r") as f:
	te_lines = list(line for line in f.read().splitlines() if line)

with open(output_file2, 'a+') as d:
	for line in te_lines:
		te = line.split("\t")[0].strip()
		species_id = line.split("\t")[1].strip()
		hit_count = int(line.split("\t")[2].strip())
		field1 = "sbatch /lustre/scratch/npaulat/yin_yang/template_extend_align_npaulat.sh /lustre/scratch/npaulat/yin_yang/masked_assemblies/"  + species_id + ".fa"
		field2 = "/lustre/scratch/npaulat/yin_yang/ext_align/mamm"
		field3 = "/lustre/scratch/npaulat/yin_yang/te_fastas/" + te + ".fa"
		out_row = "{} {} {}".format(field1, field2, field3)
		d.write(out_row + "\n")
