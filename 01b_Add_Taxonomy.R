# Load packages ####
library(dada2); packageVersion("dada2")
library(tidyverse)
library(ggplot2)
library(readxl)
library(phyloseq)
library(decontam)

# load sequenc table from previous step
seqtab.nochim <- readRDS("./output/seqtab.nochim.RDS")

# keep only fire sites and negatives - remove newly empty ASVs
seqtab.nochim <- seqtab.nochim[grep("FS-",rownames(seqtab.nochim)),]
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]


# Assign taxonomy - Custom database including LOTS of outgroups ####
taxa <- assignTaxonomy(seqs = seqtab.nochim,refFasta = "./taxonomy/UNITE_Euk_2020-02-04_non-dev.fasta.gz", multithread = TRUE,verbose = TRUE)

# save it
saveRDS(taxa,"./output/taxa.RDS")


# inspect taxonomy assignments
unique(taxa[,1])
unique(taxa[,2])
unique(taxa[,3])

# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
dim(seqtab.nochim)
#++++++++++++++++++++++++++++++++++

# Create PhyloSeq object ####

# read in metadata
meta = read_csv("./Fire_metadata_plus_GDC.csv")
# meta$controls <- meta$SampleSource
# meta$controls[meta$controls == "Extraction Negative"] <- TRUE
# meta$controls[meta$controls != TRUE] <- FALSE
# meta$controls = as.logical(meta$controls)
seqtab.nochim <- seqtab.nochim[grep("NEG",names(seqtab.nochim[,1]),invert = TRUE),]
# subset meta to just remaining ITS samples
in.meta = which(names(seqtab.nochim[,1]) %in% meta$SampleID == TRUE)
seqtab.fwd.nochim = seqtab.nochim[in.meta,]
dim(seqtab.fwd.nochim)
in.seqtab = which(meta$SampleID %in% names(seqtab.fwd.nochim[,1]))
meta = meta[in.seqtab,]

#re-order
meta = meta[order(meta$SampleID),]
row.names(meta) <- meta$SampleID

# df = as.data.frame(seqtab.fwd.nochim, row.names = 1)
# row.names(df) <- meta$SampleID
identical(row.names(seqtab.nochim), meta$SampleID)
dim(taxa)
dim(seqtab.fwd.nochim)
dim(meta)

row.names(taxa)
colnames(seqtab.nochim)


# look at ps components, make sure they match names etc.
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)

row.names(met) <- row.names(meta)
identical(sample_names(met),sample_names(otu))
sample_names(met)



# make phyloseq object 
ps <- phyloseq(otu,met,tax)

# save it
saveRDS(ps, "./output/phyloseq_object_ITS.RDS")


