## Extracting TE insertions, generating alignments and species-specific consensus sequences
  * In this case, we can build off of our original BLAST searches, but focusing on hits that met our filtering standards
  * We want to extract the TE sequences to a) make a consensus sequence, and b) use for orthologous TE site searches
 
### Extracting TE sequences from BLAST hits
  * Using the same directories we made during the original BLAST searches, make subdirectories for extracting TE insertion FASTAs
  * This uses a wrapper script, **template_extend_align_mod.sh**, for various python scripts (as before), the new one being **extract_all.py**
  * Example for job submission:
  ```
  sbatch </path/to/wrapper/script>/template_extend_align_mod5.sh </path/to/genomes>/assemblies/aJam.fa </path/to/BLAST/search/working/dir>/extend_align/mamm </path/to/TE/consensus/FASTA>/te_fastas/ArtJam-1.129.fa
  ```
  
### Copy all species-specific consensus sequences to a new directory
  * Example for mammals:
  ```
  cd </path/to/BLAST/working/dir>/extend_align/
  mkdir species_consensus_seqs
  for SUBDIR in \*/; do cp ${SUBDIR::-1}/final_consensuses/${SUBDIR::-1}\_rep.fa </path/to/BLAST/working/dir>/extend_align/species_consensus_seqs/; done
  ```
    OR
  ```
  for i in \*/; do cp ${i}/final_consensuses/* ../sp_consensus_seqs/; done
  ```
  
  * For non-mammal sequences, repeat process into a different directory (nonmamm_species_cons_seqs)
  * Copy to species_consensus_seqs directory:
  ```
  LIST=["list of all ht TEs"]
  cd </path/to/BLAST/working/dir>/extend_align/nonmamm_species_cons_seqs
  for i in \*.fa; do cp ${i} </path/to/BLAST/working/dir>/extend_align/species_consensus_seqs/; done
  cd ../species_consensus_seqs/
  for ITEM in $LIST: do mkdir ${ITEM}; mv "${ITEM}\_"* ${ITEM}/ 2>/dev/null; done
  ```
  
### Create multi-FASTA files for each element with all species-specific consensus sequences
  * Run **concatenate_te_sp_rep_seqs.py**
  * This will also output lists of elements that are a) present in only 1 species, b) present in 2+ species, which can be used for batch submissions for alignments
  
### Align consensus sequences for TEs present in 2+ species
  * This example with **te_reps_muscle_aln.sh** uses MUSCLE, but CLUSTAL works just as well
  ```
  sbatch te_reps_muscl_aln.sh
  ```
  
### Build sequence tree with RaxML (optional)
  * Generally only helpful for TEs present in 3+ species across multiple clades (i.e. bats, lemurs, and tenrecs)
  * See **raxml.sh**

### Checking for TE autonomy (optional)
  * If you wish to determine if a TE (either library consensus or species-specific consensus sequence) has open reading frames (ORFs), and do a protein BLAST search 
  * Make a list of TE_GenSpe_rep.fa files (TE by species consensus sequence file names)
  * Run EMBOSS getorf; see **nonmamm_getorfs_array.sh** for example of batch job submission
  * Example:
  ```
  </path/to/EMBOSS/installation>/EMBOSS-6.6.0/emboss/getorf -sequence </path/to/consensus/sequence/files>/extend_align/nonmamm_cons_seqs/Mariner2_pKuh_ErpCal_rep.fa -outseq </path/to/output/dir>/extend_align/nonmamm/Mariner2_pKuh_ErpCal_orfs.fa -find 3
  ```
  
### Checking for sequence matches in full TE library (optional)
  * If you wish to check that a non-mammal (or other) TE consensus sequence is either the same as the original element, or for related elements
  * See **blast_nonmamm_cons_to_lib.sh**
  
  
  
  
  
  
  
