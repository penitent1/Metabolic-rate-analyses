library(tidyverse)
library(ape)
library(rotl)
library(geiger)
library(ggthemes)
library(MCMCglmm)
library(kinship2)
library(phytools)
library(phangorn)
library(nlme)

#####################################################################
## Read sculpin phyologeny into R and prune for species in analysis
#####################################################################

ramon_tree <- read.nexus("T5000_FigTree_TimeScaledRoot1.nex")
keepers_ramon <- c("Oligocottus_maculosus",
                   "Clinocottus_globiceps",
                   "Artedius_harringtoni",
                   "Artedius_lateralis",
                   "Artedius_fenestaolis",
                   "Scorpaenichthys_marmoratus",
                   "Enophrys_bison",
                   "Hemilepidotus_hemilepidotus")
ramon_phy <- drop.tip(ramon_tree, setdiff(ramon_tree$tip.label, keepers_ramon))
ramon_phy <- compute.brlen(ramon_phy, method = "Grafen")
ramon_phy$tip.label[ramon_phy$tip.label == "Artedius_fenestaolis"] <- "Artedius_fenestralis"

sculpin_data <- read.csv(file.choose(), stringsAsFactors = FALSE,
                         strip.white = TRUE, na.strings = c("NA","."))
sculpin_data <- as_data_frame(sculpin_data)
sculpin_data <- sculpin_data[sculpin_data$species != "Blepsias_cirrhosus",]

#sculpin_data$spps_mean

inv.phylo<-inverseA(ramon_phy,nodes="TIPS",scale=TRUE)

prior2<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),
             R=list(V=1,nu=0.02))
prior2.1<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)),
             R=list(V=1,nu=0.002))
prior2.2<-list(G=list(G1=list(V=1,nu=0.2),G2=list(V=1,nu=0.2)),
               R=list(V=1,nu=0.2))
prior1<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02), G3=list(V=1,nu=0.02)),
             R=list(V=1,nu=0.02))
prior1.1<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002), G3=list(V=1,nu=0.002)),
             R=list(V=1,nu=0.002))
prior1.2<-list(G=list(G1=list(V=1,nu=0.2),G2=list(V=1,nu=0.2), G3=list(V=1,nu=0.2)),
               R=list(V=1,nu=0.02))

model_repeat1<-MCMCglmm(pcrit.r~temp,random=~species+spps+fish.id,
                        family="gaussian",
                        ginverse=list(species=inv.phylo$Ainv),
                        prior=prior3,
                        data=sculpin_data,
                        nitt=5000000,burnin=1000,thin=500)

plot(model_repeat1)
summary(model_repeat1)

model_repeat1.1<-MCMCglmm(pcrit.r~temp,random=~species+spps+fish.id,
                        family="gaussian",
                        ginverse=list(species=inv.phylo$Ainv),
                        prior=prior3.1,
                        data=sculpin_data,
                        nitt=5000000,burnin=1000,thin=500)

plot(model_repeat1.1)
summary(model_repeat1.1)

model_repeat1.2<-MCMCglmm(pcrit.r~temp,random=~species+spps+fish.id,
                          family="gaussian",
                          ginverse=list(species=inv.phylo$Ainv),
                          prior=prior1.2,
                          data=sculpin_data,
                          nitt=5000000,burnin=1000,thin=500)

plot(model_repeat1.2)
summary(model_repeat1.2)

model_repeat2<-MCMCglmm(pcrit.r~temp,random=~species+fish.id,
                        family="gaussian",
                        ginverse=list(species=inv.phylo$Ainv),
                        prior=prior2,
                        data=sculpin_data,
                        nitt=5000000,burnin=1000,thin=500)

plot(model_repeat2)
summary(model_repeat2)

model_repeat2.1<-MCMCglmm(pcrit.r~temp,random=~species+fish.id,
                        family="gaussian",
                        ginverse=list(species=inv.phylo$Ainv),
                        prior=prior2.1,
                        data=sculpin_data,
                        nitt=5000000,burnin=1000,thin=500)

plot(model_repeat2.1)
summary(model_repeat2.1)

model_3 <- MCMCglmm(pcrit.r~temp+spps+temp*spps,random=~species+fish.id,
                    family="gaussian",
                    ginverse=list(species=inv.phylo$Ainv),
                    prior=prior2,
                    data=sculpin_data,
                    nitt=5000000,burnin=1000,thin=500)

plot(model_3)
summary(model_3)

model_3.1 <- MCMCglmm(pcrit.r~temp+temp*spps,random=~species+fish.id,
                    family="gaussian",
                    ginverse=list(species=inv.phylo$Ainv),
                    prior=prior2,
                    data=sculpin_data,
                    nitt=5000000,burnin=1000,thin=500)

plot(model_3)
summary(model_3)
#lambda_m1 <- (model_repeat1$VCV[,'species']/
#  (model_repeat1$VCV[,'species']+
#     model_repeat1$VCV[,'spps']+
#     model_repeat1$VCV[,'units']
#     )
#)





#model_repeat1<-MCMCglmm(pcrit.r~temp,random=~species+spps,
#                        family="gaussian",ginverse=list(species=inv.phylo$Ainv),
#                        prior=prior2,
#                        data=sculpin_data,
#                        nitt=5000000,burnin=1000,thin=500)


#####
#####
#####
#####
phylo<-read.nexus("phylo.nex")
data<-read.table("data_repeat.txt",header=TRUE)
head(data)

data$spec_mean_cf<-sapply(split(data$cofactor,data$phylo),mean)[data$phylo]

inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)
prior2<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),
             R=list(V=1,nu=0.02))
model_repeat1<-MCMCglmm(phen~spec_mean_cf,random=~phylo+species,
                        family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                        prior=prior2,data=data,nitt=5000000,burnin=1000,thin=500)
