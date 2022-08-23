## Identification of putative HT DNA transposons involving bats

This workflow is to computationally identify putative horizontally transferred (HT) DNA transposons, done here by searching across 251 mammalian species, specifically focused on those present in any of 37 bat species.\
\
Given the wide variety of DNA transposons and species involved, this is a broad-scale search, with _a priori_ search thresholds. The main steps are:
1. Transposable element (TE) annotation using [RepeatMasker](http://repeatmasker.org/) in 253 mammalian genome assemblies
    * Detailed workflow in TE_annotation folder
2. Local [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) searches (blastn) of all DNA transposons with patchy species distributions in all mammal and other eukaryote genome assemblies; BLAST search criteria of at least 90% sequence identity and at least 90% query alignment
    * Detailed workflow in BLAST folder
3. Generate species-specific consensus sequences for all candidate DNA transposons
    * Detailed workflow in TE_alignments folder
4. Identify autonomous elements: 
     * Use [EMBOSS getorf](https://www.bioinformatics.nl/cgi-bin/emboss/help/getorf) utility to identify open-reading frames (ORFs) for all species-specific consensus sequences greater than 800 bp 
     * Perform blastx searches of non-redundant proteins for the largest ORFs from each consensus sequence (can be performed remotely or in [web browser](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome))
     * Once autonomous elements are identified, exclude any with less than 20 copies meeting the 90/90/90 criteria; exclude any non-autonomous elements with less than 100 copies
5. Generate [RAxML](https://github.com/stamatak/standard-RAxML) trees of species-specific consensus sequences for each element (optional)
     * Only useful if present in 3 or more species; used to search for phylogenetic incongruence with species tree. Since this project focused on putative HT specifically involving bats, and many elements were not found outside a single clade, this step was not particularly informative.
     * Detailed workflow in TE_alignments folder
6. Identify and exclude deletion products of elements within the set of putative HT elements by clustering sequences via a) the cross_match utility of [Phrap](http://www.phrap.org/phredphrapconsed.html) v0.990319 with default settings, and b) a modified CD-HIT search with the utility [ClusterPartialMatchingSubs.pl](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpz1.154) using default settings.
    * Detailed workflow in TE_alignments folder
7. Estimate average age of a TE subfamily in a given species using the average modified Kimura 2-parameter (K2P) distance from the library consensus and the neutral mutation rate of the species. Modified K2P values calculated in [RepeatModeler](http://repeatmasker.org/)'s utilities.
    * Detailed workflow in TE_alignments folder
8. Infer putative HT placement on bat phylogeny (see bat_phylogeny) based on presence/absence and average TE age across species.
    * Bat phylogeny based on Foley et al. 2021 and [Amador et al. 2018](https://link.springer.com/article/10.1007/s10914-016-9363-8) and used a combination of non-conflicting average or median divergnence estimates from [TimeTree](http://www.timetree.org/) (accessed 3 September 2021).
    * To be conservative, elements assigned to oldest possible branch based on presence/absence data within a given clade (i.e. there are four _Myotis_ species representing to sister species pairs (for this tree), _M. brandtii + lucifugus_ and _M. davidii + myotis_; if only _M. brandtii_ and _myotis_ were searched, and the element was found in both, it was assumed to also be in the other two species, and so the HT event would be inferred to have occurred in the ancestral _Myotis_ lineage.)
9. Estimate association between young (>50 My) TE accumulation in bats and species richess; association between putative HT diversity and species richness in bats.
    * Detailed workflow in association_modeling folder
