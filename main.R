library(dada2); packageVersion("dada2")

# Put the name of the folders of the samples demultiplexed
directory_list=c("Invierno2017")

setwd("C:/Users/apardo/Documents/metagenomica/scripts-github")
dir.create("DADA2_out")
out_path=paste(getwd(),"DADA2_out",sep = "/")

for (temporada in directory_list) {
  
  path=paste(getwd(),temporada,sep = "/")
  
  fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
  # Extract sample names
  sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
  
  # Filter and trim
  # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  ###
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,240),maxEE=c(10,10),
                       maxN=0, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE) 
  
  # Learn error rates
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  plotErrors(errF, nominalQ=TRUE)
  plotErrors(errR, nominalQ=TRUE)
  
  # Dereplication
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  # Sample Inference
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  # Merge paired reads
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, maxMismatch = 15)
  
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  
  seqtab_out=paste(out_path,"/","seqtab_",temporada,".rds", sep = "")
  saveRDS(seqtab, seqtab_out)
  
  # Track reads through the pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
  rownames(track) <- sample.names
  
  track_out=paste(out_path,"/track_",temporada, sep = "" )
  write.table(track, track_out , sep="\t")
}

# Remuevo variables
rm(list=ls())
# y la RAM: Garbage colection
gc()

library(dada2); packageVersion("dada2")
setwd("C:/Users/apardo/Documents/metagenomica/scripts-github")
out_path=paste(getwd(),"DADA2_out",sep = "/")

# Read each seqtable separately
st1 <- readRDS(paste(out_path,"seqtab_Invierno2017.rds",sep="/"))
st2 <- readRDS(paste(out_path,"seqtab_Primavera2014.rds",sep="/"))

# Merge multiple runs (if necessary)
st.all <- mergeSequenceTables(st1,st2)

#Mi agregado. Saco las quimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)

# Assign taxonomy. Correr en el cecar
taxa <- assignTaxonomy(seqtab.nochim, paste(getwd(),"silva_nr_v132_train_set.fa.gz" ,sep="/"), multithread=TRUE, tryRC=TRUE)
taxa <- addSpecies(taxa, paste(getwd(),"silva_species_assignment_v132.fa.gz", sep="/"), tryRC=TRUE)

# Write to disk
saveRDS(seqtab.nochim, paste(out_path,"seqtab_nochim.rds",sep="/")) 
saveRDS(taxa, paste(out_path, "tax_final.rds",sep="/"))
