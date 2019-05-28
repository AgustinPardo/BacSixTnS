# Bacterial(Bac)16S(SixTNS)
Pipeline develop in R to work with bacterial 16S ribosomal RNA gene amplicons from Next-generation Sequencing.
The pipeline starts from Illumina-sequenced paired-end fastq files that have been split (or “demultiplexed”) by sample and from which the barcodes/adapters have already been removed. The end product is an amplicon sequence variant (ASV) table, a higher-resolution analogue of the traditional OTU table, .rds format (RData format). 

### Starting point

This workflow assumes that your sequencing data meets certain criteria:

* Samples have been demultiplexed.
* Non-biological nucleotides have been removed, e.g. primers, adapters, linkers, etc.
* If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order.


### Filter and trim
```
filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,240),maxEE=c(10,10), maxN=0, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
```
  

### Merge paired reads
```
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, maxMismatch = 15)
```

### Track reads through the pipeline

### Remove chimeras

### Assign taxonomy

## Dependecies
DADA-2

## SILVA Data Sets
* [silva_nr_v132_train_set.fa.gz](https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1)

* [silva_species_assignment_v132.fa.gz](https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1)

