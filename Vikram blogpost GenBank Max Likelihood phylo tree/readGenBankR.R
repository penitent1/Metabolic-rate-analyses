library(ape)
packageVersion("ape")

#import each sequence list from its CSV file
#NO "strings as factor"
#coi <- read.table("C:/MyWorkingDirectory/csv_lists/coi.csv", quote="\"", stringsAsFactors=FALSE)

#name them: 
  #coi      =coi
  #cytb     =cyt


#convert to character lists
as.list(coi)$V1->coiL
as.list(cyt)$V1->cytL


#use read.GenBank to acquire the sequences

coigen<-read.GenBank(coiL,species.names=T)
cytgen<-read.GenBank(cytL,species.names=T)

#use species names
names_coi<-data.frame(species=attr(coigen,"species"),accs=names(coigen))
names(coigen)<-attr(coigen,"species")
names_cyt<-data.frame(species=attr(cytgen,"species"),accs=names(cytgen))
names(cytgen)<-attr(cytgen,"species")



#set WD

#export each

write.dna(coigen,"renamed_COI.fasta", format="fasta")
write.dna(cytgen,"renamed_CYTB.fasta", format="fasta")


