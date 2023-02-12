#!/bin/bash
#SBATCH --job-name=TEblast
#SBATCH --output=%x.%j.out
#SBATCH --partition=nocona
#SBATCH --mem-per-cpu=5150MB
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH -a 1-417

###Find non-mammalian hits from blast output files.

NAMESFILE=/lustre/scratch/npaulat/RayLib-Masking/ht_te_list
CHRANGE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $NAMESFILE)

cd /lustre/scratch/npaulat/RayLib-Masking/external_ht/

echo "BLAST run with ${CHRANGE}"

#Run a remote blast with the te in question
/lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blastn -query /lustre/scratch/npaulat/RayLib-Masking/te_fastas/$CHRANGE.fa -db /lustre/scratch/aosmansk/blastdb/ref_euk_rep_genomes -perc_identity 90 -qcov_hsp_perc 90 -out /lustre/scratch/npaulat/RayLib-Masking/external_ht/blast_asn/${CHRANGE}_90.asn -outfmt 11
#
echo "${CHRANGE} BLAST dictionary created"
#
/lustre/work/aosmansk/apps/ncbi-blast-2.11.0+/bin/blast_formatter -archive /lustre/scratch/npaulat/RayLib-Masking/external_ht/blast_asn/${CHRANGE}_90.asn -out /lustre/scratch/npaulat/RayLib-Masking/external_ht/blast_out/${CHRANGE}_90_90.out -outfmt "6 qseqid sseqid pident qcovs qcovhsp length mismatch gapopen gaps qstart qend sstart send evalue bitscore staxids sscinames scomnames"

echo "${CHRANGE} BLAST custom output created"

#Use the blast output file to find the taxonID & record the taxonomic information to a new file
for i in $(awk -F '\t' '{print $16}' /lustre/scratch/npaulat/RayLib-Masking/external_ht/blast_out/${CHRANGE}_90_90.out | sort -V | uniq); do awk -v taxID=$i '$1==taxID {print $0}' /lustre/scratch/aosmansk/external_ht/rankedlineage.dmp >> /lustre/scratch/npaulat/RayLib-Masking/external_ht/taxa_files/${CHRANGE}_90.taxa; done

echo "${CHRANGE} taxa file created"

#Delete all lines that have "Mammalia" in the .taxa file
grep -v "Mammalia" /lustre/scratch/npaulat/RayLib-Masking/external_ht/taxa_files/${CHRANGE}_90.taxa > /lustre/scratch/npaulat/RayLib-Masking/external_ht/non_mammal/${CHRANGE}_90.nonMammal

echo "${CHRANGE} nonMammal taxa list created"

#Search for Prototheria/Metatheria hits (have to do this at the Order level, as they are nested under the "Mammalia" clade
grep -e "Monotremata" -e "Dasyuromorphia" -e "Didelphimorphia" -e "Diprotodontia" -e "Microbiotheria" -e "Notoryctemorphia" -e "Paucituberculata" -e "Peramelemorphia" /lustre/scratch/npaulat/RayLib-Masking/external_ht/taxa_files/${CHRANGE}_90.taxa > /lustre/scratch/npaulat/RayLib-Masking/external_ht/metatheria/${CHRANGE}_90_nonEutherian.out

echo "${CHRANGE} nonEutherian taxa list created"
