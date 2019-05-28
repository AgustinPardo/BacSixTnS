# Bacterial(Bac)16S(SixTNS)

Pipeline develop in R to work with bacterial 16S ribosomal RNA gene amplicons from Next-generation Sequencing.
The pipeline starts from Illumina-sequenced paired-end fastq files that have been split (or “demultiplexed”) by sample and from which the barcodes/adapters have already been removed. The end product is an amplicon sequence variant (ASV) table, a higher-resolution analogue of the traditional OTU table, .rds format (RData format). 

Input

Steps

Filter

Taxonomic asignation

Data bases

**Dependecies**
+ DADA-2
