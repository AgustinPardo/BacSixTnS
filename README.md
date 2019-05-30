# Bacterial(Bac)16S(SixTNS)
Pipeline develop in R to work with bacterial 16S ribosomal RNA gene amplicons from Next-generation Sequencing.
The pipeline starts from Illumina-sequenced paired-end fastq files that have been split (or “demultiplexed”) by sample and from which the barcodes/adapters have already been removed. The end product is an amplicon sequence variant (ASV) table, a higher-resolution analogue of the traditional OTU table, .rds format (RData format). 

## Data requirements

This workflow assumes that your sequencing data meets certain criteria:

* Samples have been demultiplexed.
* Non-biological nucleotides have been removed, e.g. primers, adapters, linkers, etc.
* If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order.

## Starting point

#### Set the work directory
```R
setwd("C:/Users/apardo/Documents/metagenomica/scripts-github")
```
#### Set the samples forlders. Each folder cntain th demultiplexed files of the sample.
```R
directory_list=c("Invierno2017","Primavera2014")
```

## Filter and trim
Trim the reads using truncLen. First number is for the forwards read, and the second for the reverse read. 
Filter the reads using maxEE and maxN. Reads with higherthan maxEE "expected errors" will be discarded. Expected errors are calculatedfrom the nominal definition of the quality score: EE = sum(10^(-Q/10)). Sequences with more than maxN Ns will be discarded.
```R
filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,240), maxEE=c(10,10), maxN=0, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
```
Considerations: Your reads must still overlap after truncation in order to merge them later. That depende on the size of your     
primer set used and, and how much the forward and reverse reads overlap after trimming. Your truncLen must be large enough to maintain at least 20 nucleotides of overlap between them.

## Merge paired reads
The maximum mismatches allowed in the overlap region could be seted by maxMismatch function.
```R
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, maxMismatch = 15)
```
Considerations: Most of your reads should successfully merge. If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon. 

## Track reads through the pipeline
  The lost of reads of each step of the pipeline ("input", "filtered", "denoisedF", "denoisedR", "merged") is tracked in a file for each sample(track_sample)


### Assign taxonomy
```R
taxa <- assignTaxonomy(seqtab.nochim, paste(getwd(),"silva_nr_v132_train_set.fa.gz" ,sep="/"), multithread=TRUE, tryRC=TRUE)
```
```R
taxa <- addSpecies(taxa, paste(getwd(),"silva_species_assignment_v132.fa.gz", sep="/"), tryRC=TRUE)
```
## Dependecies
DADA-2

## SILVA Data Sets
* [silva_nr_v132_train_set.fa.gz](https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1)

* [silva_species_assignment_v132.fa.gz](https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1)

