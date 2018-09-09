library(ape)
packageVersion("ape")

cytb <- read.table("C:/Users/derek/Documents/Metabolic-rate-analyses/Vikram blogpost GenBank Max Likelihood phylo tree/cytb_sculpins_pcrit_temp_somo.csv", quote="\"", stringsAsFactors=FALSE)

#convert to character lists
as.list(cytb)$V1->cytbL

#use read.GenBank to acquire the sequences
cytbgen<-read.GenBank(cytbL,species.names=T) # object class = "DNAbin"

#use species names
names_cytb<-data.frame(species=attr(cytbgen,"species"),accs=names(cytbgen))
names(cytbgen)<-attr(cytbgen,"species")

## A dataframe with species names as one column and accession number as the other!
names_cytb

#set WD
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Vikram blogpost GenBank Max Likelihood phylo tree")

#export list as a FASTA file for further use in other programs
write.dna(cytbgen,"renamed_sculpin_practice_CYTB.fasta", format="fasta")
