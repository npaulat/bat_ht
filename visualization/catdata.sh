#!/bin/bash 
#SBATCH --job-name=catdata
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err 
#SBATCH --partition=nocona 
#SBATCH --nodes=1 
#SBATCH --mem-per-cpu=60G
#SBATCH --ntasks=1

. ~/conda/etc/profile.d/conda.sh
conda activate extend_env

cd /lustre/scratch/daray/paulat_beds2

#cp rm_beds/*_rm.bed .

#python ~/gitrepositories/bioinfo_tools/filter_beds.py -g genome_sizes_mrates.txt -p 50my -a 50000000

#python ~/gitrepositories/bioinfo_tools/filter_beds.py -g genome_sizes_mrates.txt 

#rm *rm.bed

#python ~/gitrepositories/bioinfo_tools/catdata_props_age.py -g genome_sizes_mrates.txt -p 50my

#mkdir filtered_beds
#mv *50my_filtered.bed filtered_beds/
#mv *50my_processed.bed processed_beds/
#mv *50my_processed_beds.txt processed_beds/
#mkdir processed_cats
#mv *processed_cats.txt processed_cats/

python ~/gitrepositories/bioinfo_tools/catdata_props_age.py -g genome_sizes_mrates.txt -p all

mv *all_filtered.bed filtered_beds/
mv *all_processed.bed processed_beds/
mv *all_processed_beds.txt processed_beds/
mv *processed_cats.txt processed_cats/
