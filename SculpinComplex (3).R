library(ape)
library(geiger)
library(nlme)

setwd("") #set working directory to folder

complexes <- read.table("PcritvsBrainLiverComplexes.txt", row.names=1)

str(complexes) #structure of the file 

#OLS Pcrit vs complexes ordinary least squares 

ols.PcritBCI <- lm(Pcrit~BrainCI, data=complexes) #labelled the file, tested the two, Pcrit is the y, brain CI is x, and reading the data from "complexes"
summary(ols.PcritBCI)
AIC(ols.PcritBCI) #testing which model fits the best, smallest numerical value, does go into the negatives 

ols.PcritBCII <- lm(Pcrit~BrainCII, data=complexes)
summary(ols.PcritBCII)
AIC(ols.PcritBCII)

ols.PcritBCIII <- lm(Pcrit~BrainCIII, data=complexes)
summary(ols.PcritBCIII)
AIC(ols.PcritBCIII)

ols.PcritBCIV <- lm(Pcrit~BrainCIV, data=complexes)
summary(ols.PcritBCIV)
AIC(ols.PcritBCIV)

ols.PcritBCV <- lm(Pcrit~BrainCV, data=complexes)
summary(ols.PcritBCV)
AIC(ols.PcritBCV)

ols.PcritLCI <- lm(Pcrit~LiverCI, data=complexes)
summary(ols.PcritLCI)
AIC(ols.PcritLCI)

ols.PcritLCII <- lm(Pcrit~LiverCII, data=complexes)
summary(ols.PcritLCII)
AIC(ols.PcritLCII)

ols.PcritLCIII <- lm(Pcrit~LiverCIII, data=complexes)
summary(ols.PcritLCIII)
AIC(ols.PcritLCIII)

ols.PcritLCIV <- lm(Pcrit~LiverCIV, data=complexes)
summary(ols.PcritLCIV)
AIC(ols.PcritLCIV)

ols.PcritLCV <- lm(Pcrit~LiverCV, data=complexes)
summary(ols.PcritLCV)
AIC(ols.PcritLCV)

Knope <- read.nexus("bayes_tree_final.nex") #keep it as a nexus file inside R 

plot.phylo(Knope) #really crowded tree
plot.phylo(Knope, cex=0.5) #makes the font smaller, easier to see

#To verify that the taxa listed in my data table match those at the tips of the phylogeny
name.check(Knope, complexes) #same number of species and names in both, list of species NOT matched in the data

#To prune the tree; also have the opposite way to do this re: Georgie's shark data, to select for the tips that you do want in your analyses and get the ones you do not want pruned 
delete <- c("Alcichthys_elongatus", "Archaulus_biseriatus", "Artediellus_fuscimentus", "Artediellus_pacificus", "Artediellus_scaber", "Artedius_corallinus", "Artedius_harringtoni", "Artedius_notospilotus", "Atopocottus_tribranchius", "Bero_elegans", "Blepsias_bilobus", "Chitonotus_pugetensis", "Clinocottus_acuticeps", "Clinocottus_analis", "Clinocottus_embryum", "Clinocottus_recalvus", "Cottiusculus_gonez", "Cottus_amblystomopsis", "Cottus_cognatus",  "Cottus_kazika", "Cottus_pollux", "Dasycottus_setiger", "Enophrys_diceraus", "Enophrys_lucasi", "Enophrys_taurina", "Furcina_osimae", "Gymnocanthus_galeatus", "Gymnocanthus_pistilliger", "Gymnocanthus_tricuspis", "Hemilepidotus_gilberti", "Hemilepidotus_jordani", "Hemilepidotus_papilio", "Hemilepidotus_spinosus", "Hemilepidotus_zapus", "Hemitripterus_bolini", "Icelinus_borealis", "Icelinus_burchami", "Icelinus_filamentosus", "Icelinus_fimbriatus", "Icelinus_tenuis", "Icelus_canaliculatus", "Icelus_euryops", "Icelus_ochotensis", "Icelus_spatula", "Icelus_spiniger", "Icelus_toyamensis", "Icelus_uncinalis", "Jordania_zonope", "Leiocottus_hirundo", "Megalocottus_platycephalus", "Microcottus_sellaris", "Myoxocephalus_brandtii", "Myoxocephalus_cf.scorpioides", "Myoxocephalus_jaok", "Myoxocephalus_scorpius", "Myoxocephalus_stelleri", "Myoxocephalus_verrucosus", "Nautichthys_oculofasciatus", "Nautichthys_pribilovius", "Oligocottus_rimensis", "Oligocottus_rubellio", "Oligocottus_snyderi", "Orthonopias_triacis", "Porocottus_allisi", "Porocottus_camtschaticus", "Pseudoblennius_percoides", "Pseudoblennius_sp3", "Psychrolutes_phrictus", "Radulinus_asprellus", "Rastrinus_scutiger", "Rhamphocottus_richardsonii", "Ricuzenius_pinetorum", "Ruscarius_creaseri", "Ruscarius_meanyi", "Scorpaenichthys_marmoratus", "Stellerina_xyosterna", "Stlengis_misakia", "Synchirus_gilli", "Thyriscus_anoplus", "Trichocottus_brashnikovi", "Triglops_forficatus", "Triglops_macellus", "Triglops_metopias", "Triglops_murrayi", "Triglops_nybelini", "Triglops_pingelli", "Triglops_quadricornis", "Triglops_scepticus", "Triglops_xenostethus", "Zesticelus_profundorum")

KnopeShort <- drop.tip(Knope, delete) 

#check again

name.check(KnopeShort, complexes)

#Trait measurement must be in the same row order as the names of the species in the tree object, rows need to match 

#Check again
rownames(complexessort)

____________________________________________________________________________
####Pcrit and Brain CI

#Brownian Motion model, PIC
bm.PcritBCI <- gls(Pcrit~BrainCI, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritBCI)

#Ornstein-Uhlenbeck Motion model
ou.PcritBCI <- gls(Pcrit~BrainCI, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritBCI)

#alter branch lengths, this is to get around models that gives you error message 
tempTree <- KnopeShort
tempTree$edge.length <- tempTree$edge.length*100
ou.PcritBCI <- gls(Pcrit~BrainCI, data=complexessort, correlation=corMartins(1,tempTree))
summary(ou.PcritBCI)

#Pagel's model
pagel.PcritBCI <- gls(Pcrit~BrainCI, data=complexessort, correlation=corPagel(1, KnopeShort, fixed=FALSE))
summary(pagel.PcritBCI)

pagel.PcritBCI <- gls(Pcrit~BrainCI, data=complexessort, correlation=corPagel(1, tempTree, fixed=FALSE))
summary(pagel.PcritBCI)

#Grafen's 
grafen.PcritBCI <- gls(Pcrit~BrainCI, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritBCI)

#ACDC (Blomberg)
acdc.PcritBCI <- gls(Pcrit~BrainCI, data=complexessort, correlation=corBlomberg(1, KnopeShort))
summary(acdc.PcritBCI)

____________________________________________________________________________
####Pcrit and Brain CII

#Brownian Motion model
bm.PcritBCII <- gls(Pcrit~BrainCII, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritBCII)

#Ornstein-Uhlenbeck Motion model
ou.PcritBCII <- gls(Pcrit~BrainCII, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritBCII)

#Pagel's model
pagel.PcritBCII <- gls(Pcrit~BrainCII, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritBCII)

#Grafen's 
grafen.PcritBCII <- gls(Pcrit~BrainCII, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritBCII)

#ACDC (Blomberg)
acdc.PcritBCII <- gls(Pcrit~BrainCII, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=TRUE))
summary(acdc.PcritBCII)

____________________________________________________________________________
####Pcrit and Brain CIII

#Brownian Motion model
bm.PcritBCIII <- gls(Pcrit~BrainCIII, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritBCIII)

#Ornstein-Uhlenbeck Motion model
ou.PcritBCIII <- gls(Pcrit~BrainCIII, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritBCIII)

#alter branch lengths
tempTree <- KnopeShort
tempTree$edge.length <- tempTree$edge.length*100
ou.PcritBCIII <- gls(Pcrit~BrainCIII, data=complexessort, correlation=corMartins(1,tempTree))
summary(ou.PcritBCIII)

#Pagel's model
pagel.PcritBCIII <- gls(Pcrit~BrainCIII, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritBCIII)

#Grafen's 
grafen.PcritBCIII <- gls(Pcrit~BrainCIII, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritBCIII)

#ACDC (Blomberg)
acdc.PcritBCIII <- gls(Pcrit~BrainCIII, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=T))
summary(acdc.PcritBCIII)

____________________________________________________________________________
####Pcrit and Brain CIV

#Brownian Motion model
bm.PcritBCIV <- gls(Pcrit~BrainCIV, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritBCIV)

#Ornstein-Uhlenbeck Motion model
ou.PcritBCIV <- gls(Pcrit~BrainCIV, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritBCIV)

#Pagel's model
pagel.PcritBCIV <- gls(Pcrit~BrainCIV, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritBCIV)

#Grafen's 
grafen.PcritBCIV <- gls(Pcrit~BrainCIV, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritBCIV)

#ACDC (Blomberg)
acdc.PcritBCIV <- gls(Pcrit~BrainCIV, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=T))
summary(acdc.PcritBCIV)

____________________________________________________________________________
####Pcrit and Brain CV

#Brownian Motion model
bm.PcritBCV <- gls(Pcrit~BrainCV, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritBCV)

#Ornstein-Uhlenbeck Motion model
ou.PcritBCV <- gls(Pcrit~BrainCV, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritBCV)

#Pagel's model
pagel.PcritBCV <- gls(Pcrit~BrainCV, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritBCV)

#Grafen's 
grafen.PcritBCV <- gls(Pcrit~BrainCV, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritBCV)

#ACDC (Blomberg)
acdc.PcritBCV <- gls(Pcrit~BrainCV, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=T))
summary(acdc.PcritBCV)

____________________________________________________________________________
####Pcrit and Liver CI

#Brownian Motion model
bm.PcritLCI <- gls(Pcrit~LiverCI, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritLCI)

#Ornstein-Uhlenbeck Motion model
ou.PcritLCI <- gls(Pcrit~LiverCI, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritLCI)

#Pagel's model
pagel.PcritLCI <- gls(Pcrit~LiverCI, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritLCI)

#Grafen's 
grafen.PcritLCI <- gls(Pcrit~LiverCI, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritLCI)

#ACDC (Blomberg)
acdc.PcritLCI <- gls(Pcrit~LiverCI, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=T))
summary(acdc.PcritLCI)

____________________________________________________________________________
####Pcrit and Liver CII

#Brownian Motion model
bm.PcritLCII <- gls(Pcrit~LiverCII, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritLCII)

#Ornstein-Uhlenbeck Motion model
ou.PcritLCII <- gls(Pcrit~LiverCII, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritLCII)

#Pagel's model
pagel.PcritLCII <- gls(Pcrit~LiverCII, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritLCII)

#Grafen's 
grafen.PcritLCII <- gls(Pcrit~LiverCII, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritLCII)

#ACDC (Blomberg)
acdc.PcritLCII <- gls(Pcrit~LiverCII, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=TRUE))
summary(acdc.PcritLCII)

____________________________________________________________________________
####Pcrit and Liver CIII

#Brownian Motion model
bm.PcritLCIII <- gls(Pcrit~LiverCIII, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritLCIII)

#Ornstein-Uhlenbeck Motion model
ou.PcritLCIII <- gls(Pcrit~LiverCIII, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritLCIII)

#alter branch lengths
tempTree <- KnopeShort
tempTree$edge.length <- tempTree$edge.length*100
ou.PcritLCIII <- gls(Pcrit~LiverCIII, data=complexessort, correlation=corMartins(1,tempTree))
summary(ou.PcritLCIII)

#Pagel's model
pagel.PcritLCIII <- gls(Pcrit~LiverCIII, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritLCIII)

#Grafen's 
grafen.PcritLCIII <- gls(Pcrit~LiverCIII, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritLCIII)

#ACDC (Blomberg)
acdc.PcritLCIII <- gls(Pcrit~LiverCIII, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=T))
summary(acdc.PcritLCIII)

____________________________________________________________________________
####Pcrit and Liver CIV

#Brownian Motion model
bm.PcritLCIV <- gls(Pcrit~LiverCIV, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritLCIV)

#Ornstein-Uhlenbeck Motion model
ou.PcritLCIV <- gls(Pcrit~LiverCIV, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritLCIV)

#Pagel's model
pagel.PcritLCIV <- gls(Pcrit~LiverCIV, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritLCIV)

#Grafen's 
grafen.PcritLCIV <- gls(Pcrit~LiverCIV, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritLCIV)

#ACDC (Blomberg)
acdc.PcritLCIV <- gls(Pcrit~LiverCIV, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=T))
summary(acdc.PcritLCIV)

____________________________________________________________________________
####Pcrit and Liver CV

#Brownian Motion model
bm.PcritLCV <- gls(Pcrit~LiverCV, data=complexessort, correlation=corBrownian(1, KnopeShort))
summary(bm.PcritLCV )

#Ornstein-Uhlenbeck Motion model
ou.PcritLCV  <- gls(Pcrit~LiverCV, data=complexessort, correlation=corMartins(1,KnopeShort))
summary(ou.PcritLCV)

#alter branch lengths
tempTree <- KnopeShort
tempTree$edge.length <- tempTree$edge.length*100
ou.PcritLCV <- gls(Pcrit~LiverCV, data=complexessort, correlation=corMartins(1,tempTree))
summary(ou.PcritLCV)

#Pagel's model
pagel.PcritLCV <- gls(Pcrit~LiverCV, data=complexessort, correlation=corPagel(1, KnopeShort))
summary(pagel.PcritLCV)

#Grafen's 
grafen.PcritLCV <- gls(Pcrit~LiverCV, data=complexessort, correlation=corGrafen(1, KnopeShort))
summary(grafen.PcritLCV)

#ACDC (Blomberg)
acdc.PcritLCV <- gls(Pcrit~LiverCV, data=complexessort, correlation=corBlomberg(1, KnopeShort, fixed=T))
summary(acdc.PcritLCV)
