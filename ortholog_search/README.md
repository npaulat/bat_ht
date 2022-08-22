## Orthologous TE insertion site searches across genomes
  * Required programs: [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/), [BEDtools](https://bedtools.readthedocs.io/en/latest/), halLiftover, hal2fasta, [Cross_Match](https://www.repeatmasker.org/RepeatMasker/)
  * This requires a whole-genome HAL alignment of the taxa of interest
  * In this case, the [241-way mammalian alignment](https://cglgenomics.ucsc.edu/data/cactus/) (v2) from the Zoonomia project was used

### Generate BED files of putative HTT insertions with flanking sequence (in this case 50 bp) for each taxon
  * This is output from **reformat_filtered_blast.py**

### Generate Liftover BED files of orthologous sites in pairwise species comparisons from the HAL alignment; see **halLiftover_array.sh**
  * Example:
  ```
  ./halLiftover --inBedVersion 3 --keepExtra /lustre/scratch/npaulat/RayLib-Masking/241-mammalian-2020v2.hal Myotis_myotis /lustre/scratch/RayLib-Masking/hal_comp/filtered_beds/MyoMyo_blast90_flank50.bed Microgale_talazaci MyoMyo_to_MicTal.bed
  ```
  
### Group lifts by coordinates and merge overlapping/near-overlapping (e.g. 2 bp distant) lifts
  * This requires BEDtools
  * Example:
  ```
  cd /lustre/scratch/npaulat/RayLib-Masking/hal_comp/lifts
  
  for i in *.bed; 
    do NAME=$(basename ${i} .bed); 
    sort -k1,1 -k2,2n ${i} | awk '{print$1"\t"$2"\t"$3"\t"$4"*"$5}' | sort -k1,1 -k2,2n | groupBy -grp 4 -c 2,3 -o min,max -full | awk '{print$1"\t"$5"\t"$6"\t"$4}' | awk '{print$1"\t"$2"\t"$3"\t"$3-$2"\t"$4}' | sed 's/:/\t/g' | sed 's/*/\t/g' | awk '{gsub("-","\t",$6)}1' | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$7-$6"\t"$8}' | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$8-100"\t"$9"\t"$8-$4"\t"$9-$4}' > /lustre/scratch/npaulat/RayLib-Masking/hal_comp/merged_lifts/${NAME}_flank50.merged.bed; 
  done
  ```
  
### Convert HAL alignment to aligned taxon genome FASTA
  * Example using *Bos taurus*:
  ```
  ./hal2fasta /u1/local/genomes/200Mammals/Version2/241-mammalian-2020v2.hal Bos_taurus > BosTau_genome.fa
  ```
  * Example using array job syntax (requires two lists, one of genus_species, another of GenSpe):
  ```
  ./hal2fasta /lustre/scratch/npaulat/RayLib-Masking/241-mammalian-2020v2.hal ${CHRANGE1} > /lustre/scratch/npaulat/RayLib-Masking/hal_genomes/${CHRANGE2}_genome.fa
  ```

### Extract the FASTA sequences for the lifted orthologous sites from the HAL alignment
  * Note: For some reason, the following bash command does not work as copy/paste, but will work if typed manually
  * Example with *Tonatia saurophila* HTT sites (source taxon) in *Bos taurus* (target taxon):
  ```
  bedtools getfasta -fi BosTau_genome.fa -bed TonSau_to_BosTau_unsure.bed -fo TonSau_to_BosTau_unsure.fa
  ```
  * Example using array syntax (uses all files ending in \*\_flank50.merged.bed in the current directory):
  ```
  cd merged_lifts

  for i in *_flank50.merged.bed; 
    do id=$(awk -F_ '{print $3}' <<< ${i}); id2=$(basename ${i} .merged.bed); 
    bedtools getfasta –fi /lustre/scratch/npaulat/RayLib-Masking/hal_genomes/${id}_genome.fa –bed /lustre/scratch/npaulat/RayLib-Masking/hal_comp/merged_lifts/${i} –fo /lustre/scratch/npaulat/RayLib-Masking/hal_comp/unsure_fa/${id2}.merged.fa; 
  done
  ```
  
### Run RepeatMasker on the lifted orthologous site FASTA sequences using custom TE library FASTA file (-lib)
  * See library prep steps further down if combining/creating a library from Dfam database RepeatMasker's famdb.py utility
  * Generic example:
  ```
  RepeatMasker-4.1.2-p1/RepeatMasker -s -lib Dfam3.6_curated_wNames.fa TonSau_to_BosTau_unsure.fa -pa 14 –nolow
  ```
  * Job submission example:
  ```
  cd unsure_fa/
  /lustre/work/daray/software/RepeatMasker-4.1.2-p1/RepeatMasker -s -lib /lustre/scratch/npaulat/RayLib-Masking/combined_library.fa ${CHRANGE} -pa 14 –nolow
  ```
### Make list file of putative HTT names from HTT library
```
grep ">" ht_te_221_library.fa | sed 's/>//g' > 221titles
```

### Pull out any target HTT IDs from the RepeatMasker output
  * This uses diff and fgrep to first find non-matching IDs, then making a subset that have the HTT IDs
  ```
  cp *flank50.merged.fa.out ../merged_rm/merged50_out/
  cd ../merged_rm/merged50_out
  cd /lustre/scratch/npaulat/RayLib-Masking/hal_comp/merged_rm/merged50_out
  for i in *BosTau_flank50.merged.fa.out; do id=$(basename ${i} .fa.out); fgrep -vf ../../ht_te_221_list ${i} > ../merged50_check/${id}_no221.out ; done
  
  cd ../merged50_check
  for i in *BosTau_flank50.merged_no221.out; do id=$(basename ${i} _no221.out); diff ../merged50_out/${id}.fa.out ${i} > ${id}_221match.out; done
  ```

### For any matching IDs, pull the coordinates and then the relevant FASTA sequences
```
for i in *_221match.out; do id=$(basename ${i} .out); grep "<" ${i} | awk '{print $6"\t"$7"\t"$8}' | sed 's/:/\t/g' | sed 's/-/\t/g' | awk '{print$1"\t"$2+$4"\t"$2+$5}' > ${id}_lines.bed; done
for i in *_221match_lines.bed; do id=$(basename ${i} .bed); id2=$(awk -F_ '{print $3}' <<< ${i}); bedtools getfasta -fi ../../../hal_genomes/${id2}_genome.fa -bed ${i} -fo ${id}.fa; done
```

### For possible matches, run CrossMatch to confirm TE identity
  * In cases where you wish to see the alignment (i.e. results show an HTT match), rerun CrossMatch with the -alignments flag
```
for i in *_221match_lines.fa; do if [ -s ${i} ]; then id=$(basename ${i} _lines.fa); /lustre/work/daray/software/cross_match/cross_match ${i} ../../fastalib > ${id}_crossmatch.out; fi; done
```

### Reformat CrossMatch output and check for HTT IDs
```
for i in *221match_crossmatch.out; do id=$(basename ${i} .out); /lustre/scratch/npaulat/RayLib-Masking/GrepCrossmatch ${i} > ${id}_hits.out; done

for i in *221match_crossmatch_hits.out; do if [ -s ${i} ]; id=$(basename ${i} _hits.out); then fgrep -f ../../ht_te_221_list ${i} > ${id}_match.out; fi; done
```

### For any hits with HTT IDs, check results to source annotations
  * Check the source_to_target.merged.bed to see if TE matches that in the source taxon
    * If it doesn't match the TE ID or size, congratulations, you found a parallel insertion, not an orthologous insertion
  * Can also check for large gaps in the lifted site
    * If large gap (e.g. 10 or 100 kb), there is a gap in the alignment, and the "match" is actually sequence bits from two very different locations OR at the very least from a likely non-orthologous site
  * Can also rerun CrossMatch to include the alignment of the TE consensus sequence to the target sequence to check the quality of the ID
    
&nbsp;  
&nbsp;
&nbsp;
### If creating/combining libraries with Dfam database elements, will need to get IDs and edit headers "name=#DFXXXXX.XX" with actual TE name
```
/lustre/work/daray/software/RepeatMasker-4.1.2-p1/famdb.py –i Dfam_curatedonly.h5 families –ad root –f fasta acc > Dfam3.6_curated.fa
sed 's/ name=/_/g' Dfam3.6_curated.fa > Dfam3.6_curated_names.fa
```

If needed, combine Dfam elements with another library (in this case, the full mammal library) using BioPython interactively:
```  
from Bio import SeqIO

dfam="/lustre/scratch/npaulat/RayLib-Masking/hal_comp/Dfam3.6_curated_names.fa"
zoo="/lustre/scratch/npaulat/RayLib-Masking/te_fastas/final_mammal_library_reduced.fa"
zoo_ids=[]
with open(zoo) as handle:
	for record in SeqIO.parse(handle, “fasta”):
		zoo_ids.append(record.id)

dfam_ids=[]
with open(dfam) as handle:
	for record in SeqIO.parse(handle, “fasta”):
		dfam_ids.append(record.id)

dfam_names=[]
for x in dfam_ids:
	if ‘_’ in x:
		name=x.split(‘_’,1)[1]
		dfam_names.append(name)
	else:
		dfam_names.append(x)

out1=”dfam_ids_list”
with open(out1, “w”) as file:
	for x in dfam_ids:
		file.write(x+”\n”)
out1=”dfam_names_list”
with open(out1, “w”) as file:
	for x in dfam_names:
		file.write(x+”\n”)
out1=”zoo_ids_list”
with open(out1, “w”) as file:
	for x in zoo_ids:
		file.write(x+”\n”)

new_ids=[]
for id in zoo_ids:
	if id not in dfam_names:
		new_ids.append(id)

out1=”zoo_only_ids_list”
with open(out1, “w”) as file:
	for x in new_ids:
		file.write(x+”\n”)

combined=[]
with open(dfam) as handle:
	for record in SeqIO.parse(handle, “fasta”):
	combined.append(record)

rec_dict = SeqIO.index(zoo, “fasta”)
for item in new_ids:
	record = rec_dict[item]
	combined.append(record)

with open(combined_library.fa”, “w”) as handle:
	SeqIO.write(combined, handle, “fasta”)
```
  
  
