library(ape)
library(MCMCglmm)
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/MCMCglmm practice")
phylo <- read.nexus("phylo.nex")
data <- read.table("data_simple.txt", header = TRUE)
head(data)

inv.phylo <- inverseA(phylo, nodes = "TIPS", scale = TRUE)
prior <- list(G=list(G1=list(V=1, nu=0.02)), R=list(V=1, nu=0.02))

model_simple <- MCMCglmm(phen ~ cofactor, random = ~phylo,
                         family = "gaussian", ginverse=list(phylo=inv.phylo$Ainv), prior=prior,
                         data = data, nitt=5000000, burnin = 1000, thin = 500)

## ***********************
##
## Using MCMCglmm package
##
## ***********************

inv.phylo <- inverseA(phylo, nodes = "TIPS", scale = TRUE)
prior <- list(G=list(G1=list(V=1, nu=0.02)), R=list(V=1, nu=0.02))

model_simple <- MCMCglmm(phen ~ cofactor, random = ~phylo,
                         family = "gaussian", ginverse=list(phylo=inv.phylo$Ainv), prior=prior,
                         data = data, nitt=5000000, burnin = 1000, thin = 500)

mcmcglmm_data <- lm_pcrit_smr_temp %>%
  dplyr::select(species, slope_pcrit_low_temps, slope_pcrit_high_temps, mean_pcrit_12) %>% 
  as.data.frame()

inv_phylo_mandic <- inverseA(mandic_phy, nodes = "TIPS") #, scale = TRUE
inv_phylo_mandic_all <- inverseA(mandic_phy, nodes = "ALL")

inv_phylo_simtree <- inverseA(tree, nodes = "TIPS")

prior <- list(G=list(G1=list(V=1, nu=0.02)), R=list(V=1, nu=0.02))
mcmcglmm_beta_low_temps_pcrit12 <- MCMCglmm(slope_pcrit_low_temps ~ mean_pcrit_12,
                                            random = ~ species,
                                            family = "gaussian",
                                            ginverse = list(species=inv_phylo_mandic$Ainv),
                                            prior = prior,
                                            data = mcmcglmm_data,
                                            nitt = 5000000, burnin = 1000, thin = 500)

summary(mcmcglmm_beta_low_temps_pcrit12)

mcmcglmm_beta_high_temps_pcrit12 <- MCMCglmm(slope_pcrit_high_temps ~ mean_pcrit_12,
                                            random = ~ species,
                                            family = "gaussian",
                                            ginverse = list(species=inv_phylo_mandic$Ainv),
                                            prior = prior,
                                            data = mcmcglmm_data,
                                            nitt = 5000000, burnin = 1000, thin = 500)

summary(mcmcglmm_beta_high_temps_pcrit12)

## CTmax ~ high temp Beta Pcrit

mcmcglmm_ctmax_data <- beta_pcrit_ct_max_data %>%
  dplyr::select(species, slope_pcrit_low_temps, slope_pcrit_high_temps, mean_pcrit_12, avg_ct_max) %>% 
  as.data.frame()

inv_phylo_mandic_ctmax <- inverseA(mandic_ctmax_phy, nodes = "TIPS")

mcmcglmm_beta_high_temps_ctmax <- MCMCglmm(slope_pcrit_high_temps ~ avg_ct_max,
                                             random = ~ species,
                                             family = "gaussian",
                                             ginverse = list(species=inv_phylo_mandic_ctmax$Ainv),
                                             prior = prior,
                                             data = mcmcglmm_ctmax_data,
                                             nitt = 5000000, burnin = 1000, thin = 500)

summary(mcmcglmm_beta_high_temps_ctmax)

obj_lambda<-(obj$VCV[,1])/(obj$VCV[,1]+obj$VCV[,2])
plot(obj_lambda);mean(obj_lambda)

mcmcglmm_beta_high_temps_ctmax_lambda <- (mcmcglmm_beta_high_temps_ctmax$VCV[,1])/(mcmcglmm_beta_high_temps_ctmax$VCV[,1]+
                                                          mcmcglmm_beta_high_temps_ctmax$VCV[,2])
plot(mcmcglmm_beta_high_temps_ctmax_lambda);mean(mcmcglmm_beta_high_temps_ctmax_lambda)

## Pcrit ~ MO2 + TPO + SPPS i.e. phylo signal = spps

mcmcglmm_pcrit_mo2_tpo_data <- mass_corr_smr_pcrit_data %>%
  dplyr::select(species, pcrit.r, smr.mass.corr.ms) %>%
  mutate(phylo = species,
         species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus"),
         tpo = case_when(species == "Oligocottus_maculosus" ~ "yes",
                         species == "Clinocottus_globiceps" ~ "yes",
                         species == "Artedius_harringtoni" ~ "yes",
                         species == "Artedius_lateralis" ~ "yes",
                         species == "Artedius_fenestralis" ~ "yes",
                         species == "Blepsias_cirrhosus" ~ "no",
                         species == "Enophrys_bison" ~ "no",
                         species == "Hemilepidotus_hemilepidotus" ~ "no",
                         species == "Scorpaenichthys_marmoratus" ~ "no")) %>%
  as.data.frame()

## Eq 6 from Ch 11 text: MPCM Evolution book
data$spec_mean_cf<-sapply(split(data$cofactor,data$phylo),mean)[data$phylo]

mcmcglmm_pcrit_mo2_tpo_data$spec_mean_cf <- sapply(split(mcmcglmm_pcrit_mo2_tpo_data$smr.mass.corr.ms,
                                                         mcmcglmm_pcrit_mo2_tpo_data$phylo),
                                                   mean)[mcmcglmm_pcrit_mo2_tpo_data$phylo]

inv_phylo_mandic <- inverseA(mandic_phy, nodes = "TIPS")

prior2<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

mcmcglmm_pcrit_mo2_tpo <- MCMCglmm(pcrit.r ~ smr.mass.corr.ms + tpo,
                                           random = ~ phylo+species,
                                           family = "gaussian",
                                           ginverse = list(species=inv_phylo_mandic$Ainv),
                                           prior = prior2,
                                           data = mcmcglmm_pcrit_mo2_tpo_data,
                                           nitt = 5000000, burnin = 1000, thin = 500)

summary(mcmcglmm_pcrit_mo2_tpo)

mcmcglmm_pcrit_mo2_tpo_lambda <- (mcmcglmm_pcrit_mo2_tpo$VCV[,1])/(mcmcglmm_pcrit_mo2_tpo$VCV[,1]+
                                                                     mcmcglmm_pcrit_mo2_tpo$VCV[,2])
plot(mcmcglmm_pcrit_mo2_tpo_lambda);mean(mcmcglmm_pcrit_mo2_tpo_lambda)

