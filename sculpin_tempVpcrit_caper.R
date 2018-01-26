library(ctv)
library(tidyverse)
library(ape)
library(rotl)
library(geiger)
library(ggthemes)
library(MCMCglmm)
library(kinship2)
library(phytools)

#####################################################################
## Read sculpin phyologeny into R and prune for species in analysis
#####################################################################

knope <- read.tree("knope_rooted_ultrametric")
#knope.tree.new <- read.nexus("bayes_tree_final_NEW.nex")
#knope2 <- read.nexus("bayes_tree_final_knope.nex")

is.ultrametric(knope)
is.rooted(knope)

plot.phylo(knope, cex = 0.5)

## Species I have Pcrit data for, plus one outgroup to keep tree rooted
keepers <- c("Oligocottus_maculosus",
             "Clinocottus_globiceps",
             "Artedius_harringtoni",
             "Artedius_lateralis",
             "Artedius_fenestralis",
             "Blepsias_cirrhosus",
             "Scorpaenichthys_marmoratus",
             "Enophrys_bison",
             "Hemilepidotus_hemilepidotus")

## Prune phylogeny for species in analysis
knope_pruned <- drop.tip(knope, setdiff(knope$tip.label, keepers))
plot.phylo(knope_pruned)

############################################
## Read sculpin pcrit and meta-data into R
############################################

md <- read.csv("SculpinPcritData_ComparativeAnalysisFormat.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
md <- as_tibble(md)

## Pcrit and temp data, NA's filtered out, only "trial.no == 1" inlcuded
md_pcrit_smr <- md %>%
  select(species, spps, date, temp, fish.id, mass.g, smr.best, pcrit.r, trial.no) %>%
  filter(!is.na(smr.best)) %>%
  filter(!is.na(pcrit.r)) %>%
  filter(trial.no == 1)

## Averages for Pcrit and SMR per species per temperature, all 3 temps included
md_all_temps <- md_pcrit_smr %>%
  group_by(temp, species) %>%
  summarise(avg_pcrit = mean(pcrit.r), sd_pcrit = sd(pcrit.r), n_pcrit = length(pcrit.r),
            avg_smr = mean(smr.best), sd_smr = sd(smr.best), n_smr = length(smr.best)) %>%
  group_by(temp, species) %>%
  mutate(sem_pcrit = sd_pcrit/n_pcrit,
         sem_smr = sd_smr/n_smr)

## Averages Pcrit, smr per species per temperature, no data at 16C (prep for q10)
md_temp_groups_q10 <- md_pcrit_smr %>%
  select(species, spps, date, temp, fish.id, mass.g, smr.best, pcrit.r, trial.no) %>%
  filter(temp != 16) %>%
  group_by(temp, species) %>%
  summarise(avg_pcrit = mean(pcrit.r), sd_pcrit = sd(pcrit.r),
            avg_smr = mean(smr.best), sd_smr = sd(smr.best))

## Q10s for pcrit and smr
md_q10 <- md_temp_groups_q10 %>%
  group_by(species) %>%
  summarise(q10_pcrit = (avg_pcrit[temp == 20]/avg_pcrit[temp == 12])^
              (10/(temp[temp == 20]-temp[temp == 12])),
            q10_smr = (avg_smr[temp == 20]/avg_smr[temp == 12])^
              (10/(temp[temp == 20]-temp[temp == 12])))

## Open ct_max data and add to q10 data
ct_max_depth_df <- read.csv("ct_max_data.csv", stringsAsFactors = FALSE,
                      strip.white = TRUE, na.strings = c("NA","."))
md_q10_ctmax <- full_join(md_q10, ct_max_df, by = "species")
md_q10_pcrit12 <- md_temp_groups_q10[md_temp_groups_q10$temp==12,]

## Use this df to test Pcrit at 12 against Q10 for Pcrit
md_q10_ctmax_pcrit12 <- left_join(md_q10_ctmax, md_pcrit12, by = "species")

############################################
############################################



## Statistical models of trait correlations/regression



############################################
############################################

##########
##########
##########

## Pcrit at 12C vs max depth: All 9 species!

md_q10_ctmax_pcrit12 <- md_q10_ctmax_pcrit12[md_q10_ctmax_pcrit12$species!="Rhamphocottus_richardsoni",]
md_pcrit_depth <- md_q10_ctmax_pcrit12[,c("species","max_depth_mandic","avg_pcrit")]
md_pcrit_depth$max_depth_mandic[md_pcrit_depth$species=="Blepsias_cirrhosus"] <- 150
md_pcrit_depth$max_depth_mandic[md_pcrit_depth$species=="Hemilepidotus_hemilepidotus"] <- 450
md_pcrit_depth <- as.data.frame(md_pcrit_depth)
md_pcrit_depth_rn <- column_to_rownames(md_pcrit_depth, var = "species")

model_pcrit_depth_data <- comparative.data(knope_pruned,
                                         md_pcrit_depth,
                                         names.col = "species",
                                         force.root = TRUE,
                                         vcv = TRUE,
                                         warn.dropped = TRUE)
model_pgls_pcrit_depth <- pgls(avg_pcrit~max_depth_mandic, model_pcrit_depth_data, lambda = "ML")
model_profile_pgls <- pgls.profile(model_pgls_pcrit_depth, "lambda", N=100, param.CI = 0.95)
plot(model_pgls_pcrit_depth)
summary(model_pgls_pcrit_depth)
lm.lk_pcritVdepth<-pgls.profile(model_pgls_pcrit_depth, which="lambda")
plot(lm.lk_pcritVdepth)# Put these in a slide after end of SICB talk
res_pgls <- residuals.pgls(model_pgls_pcrit_depth)
plot(model_pgls_pcrit_depth)
qqnorm(res_pgls)
qqline(res_pgls)
plot(res_pgls)
shapiro.test(res_pgls)
#> shapiro.test(res_pgls)
#Shapiro-Wilk normality test
#data:  res_pgls
#W = 0.96051, p-value = 0.8037

model_ols_pcrit_depth <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit_depth_rn, method="ML")
summary(model_ols_pcrit_depth)
res_ols <- residuals(model_ols_pcrit_depth)
plot(model_ols_pcrit_depth)
qqnorm(res_ols)
qqline(res_ols)
plot(res_ols)
shapiro.test(res_ols)
#> shapiro.test(res_ols)
#Shapiro-Wilk normality test
#data:  res_ols
#W = 0.96051, p-value = 0.8037

AIC(model_pgls_pcrit_depth)
AIC(model_ols_pcrit_depth)

anova.pgls(model_pgls_pcrit_depth)
#> anova.pgls(model_pgls_pcrit_depth)
#Analysis of Variance Table
#Sequential SS for pgls: lambda = 0.00, delta = 1.00, kappa = 1.00
#
#Response: avg_pcrit
#Df Sum Sq Mean Sq F value Pr(>F)
#max_depth_mandic  1  7.532  7.5320    1.37 0.2801
#Residuals         7 38.485  5.4979

anova(model_ols_pcrit_depth)

pcrit_depth_plot <- ggplot(md_pcrit_depth, aes(x=max_depth_mandic, y=avg_pcrit)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd, size=1.25)  +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Max depth (m)"),
       y = expression("P"["crit"]~("Torr"))) +
  #geom_abline(slope = 0.03985455, intercept = 29.61927908, 
  #            color="black", size=1) +
  #geom_abline(slope = 0.05706547, intercept = 28.81887854,
  #            color="blue", size=1) +
  coord_cartesian(ylim = c(10,60))



  ## Q10_Pcrit vs max depth: All 9 species!
  
md_q10_pcrit_depth <- md_q10_ctmax_pcrit12[,c("species","max_depth_mandic","q10_pcrit")]
md_q10_pcrit_depth$max_depth_mandic[md_pcrit_depth$species=="Blepsias_cirrhosus"] <- 150
md_q10_pcrit_depth$max_depth_mandic[md_pcrit_depth$species=="Hemilepidotus_hemilepidotus"] <- 450
md_q10_pcrit_depth <- as.data.frame(md_q10_pcrit_depth)
md_q10_pcrit_depth_rn <- column_to_rownames(md_q10_pcrit_depth, var = "species")

model_q10_pcrit_depth_data <- comparative.data(knope_pruned,
                                           md_q10_pcrit_depth,
                                           names.col = "species",
                                           force.root = TRUE,
                                           vcv = TRUE,
                                           warn.dropped = TRUE)
model_pgls_q10_pcrit_depth <- pgls(q10_pcrit~max_depth_mandic, model_q10_pcrit_depth_data, lambda = "ML")
model_profile_pgls_q10_pcrit_depth <- pgls.profile(model_pgls_q10_pcrit_depth, "lambda", N=100, param.CI = 0.95)
plot(model_pgls_q10_pcrit_depth)
summary(model_pgls_q10_pcrit_depth)
lm.lk_q10_pcritVdepth<-pgls.profile(model_pgls_q10_pcrit_depth, which="lambda")
plot(lm.lk_q10_pcritVdepth)# Put these in a slide after end of SICB talk

res_pgls_q10 <- residuals.pgls(model_pgls_q10_pcrit_depth)
plot(model_pgls_q10_pcrit_depth)
qqnorm(res_pgls_q10)
qqline(res_pgls_q10)
#plot(res_pgls_q10)
shapiro.test(res_pgls_q10)

anova.pgls(model_pgls_q10_pcrit_depth)

model_ols_q10_pcrit_depth <- gls(q10_pcrit ~ max_depth_mandic, data=md_q10_pcrit_depth_rn, method="ML")
summary(model_ols_q10_pcrit_depth)
res_ols_q10 <- residuals(model_ols_q10_pcrit_depth)
plot(model_ols_q10_pcrit_depth)
qqnorm(res_ols_q10)
qqline(res_ols_q10)
plot(res_ols_q10)
shapiro.test(res_ols_q10)

anova(model_ols_q10_pcrit_depth)

q10_pcrit_depth_plot <- ggplot(md_q10_pcrit_depth, aes(x=max_depth_mandic, y=q10_pcrit)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd, size=1.25)  +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Max depth (m)"),
       y = expression("Q"[10]~"of"~"P"["crit"])) #+
#geom_abline(slope = 0.03985455, intercept = 29.61927908, 
#            color="black", size=1) +
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(10,60))

###   ###   ###   ###   ###   ###   ###   ###
   ###   ###   ###   ###   ###   ###   ###

## Q10_Pcrit vs Pcrit at 12C: All 9 species!

md_q10pcrit_pcrit12 <- md_q10_ctmax_pcrit12[,c("species","q10_pcrit","avg_pcrit")]
md_q10pcrit_pcrit12 <- as.data.frame(md_q10pcrit_pcrit12)
md_q10pcrit_pcrit12_rn <- column_to_rownames(md_q10pcrit_pcrit12, var = "species")

model_q10pcrit_pcrit12_data <- comparative.data(knope_pruned,
                                                md_q10pcrit_pcrit12,
                                                names.col = "species",
                                                force.root = TRUE,
                                                vcv = TRUE,
                                                warn.dropped = TRUE)
model_pgls_q10pcrit_pcrit12 <- pgls(q10_pcrit~avg_pcrit, model_q10pcrit_pcrit12_data, lambda = "ML")
model_profile_pgls_q10pcrit_pcrit12 <- pgls.profile(model_pgls_q10pcrit_pcrit12, "lambda", N=100, param.CI = 0.95)
plot(model_pgls_q10pcrit_pcrit12)
summary(model_pgls_q10pcrit_pcrit12)
lm.lk_q10pcrit_pcrit12<-pgls.profile(model_pgls_q10pcrit_pcrit12, which="lambda")
plot(lm.lk_q10pcrit_pcrit12)# Put these in a slide after end of SICB talk

res_pgls_q10pcrit_pcrit12 <- residuals.pgls(model_pgls_q10pcrit_pcrit12)
plot(model_pgls_q10_pcrit_depth)
qqnorm(res_pgls_q10pcrit_pcrit12)
qqline(res_pgls_q10pcrit_pcrit12)
#plot(res_pgls_q10)
shapiro.test(res_pgls_q10pcrit_pcrit12) # Passed test of normality of model residuals
#anova.pgls(model_pgls_q10_pcrit_depth)

model_ols_q10pcrit_pcrit12 <- gls(q10_pcrit ~ avg_pcrit, data=md_q10pcrit_pcrit12_rn, method="ML")
summary(model_ols_q10pcrit_pcrit12)
res_ols_q10pcrit_pcrit12 <- residuals(model_ols_q10pcrit_pcrit12)
plot(model_ols_q10pcrit_pcrit12)
qqnorm(res_ols_q10pcrit_pcrit12)
qqline(res_ols_q10pcrit_pcrit12)
plot(res_ols_q10pcrit_pcrit12)
shapiro.test(res_ols_q10pcrit_pcrit12) # Passed test of normality of model residuals
anova(model_ols_q10pcrit_pcrit12)

x_min <- min(md_q10pcrit_pcrit12$avg_pcrit)
x_max <- max(md_q10pcrit_pcrit12$avg_pcrit)
y_min <- 3.1143425 + -0.0327242*(min(md_q10pcrit_pcrit12$avg_pcrit))
y_max <- 3.1143425 + -0.0327242*(max(md_q10pcrit_pcrit12$avg_pcrit))

md_pcrit12sem <- md_all_temps[md_all_temps$temp==12,
                              c("species","sem_pcrit")]
md_q10pcrit_pcrit12sem <- full_join(md_q10pcrit_pcrit12, md_pcrit12sem, by="species")

q10pcrit_pcrit12_plot <- ggplot(md_q10pcrit_pcrit12sem, aes(x=avg_pcrit, y=q10_pcrit)) +
  geom_errorbarh(mapping = aes(xmax = avg_pcrit+sem_pcrit, xmin = avg_pcrit-sem_pcrit, height = 0.1), size=1) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.2, size=1) +
  #geom_line(position=pd, size=1.25)  +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("P"["crit"]~("Torr")),
       y = expression("Q"[10]~"of"~"P"["crit"])) +
  geom_segment(aes(x=x_min, y=y_min, xend=x_max, yend=y_max),
               color = "red", size=1)
  #geom_abline(slope = -0.0327242, intercept = 3.1143425, 
  #            color="red", size=1) #+
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(10,60))

###   ###   ###   ###   ###   ###   ###   ###
###   ###   ###   ###   ###   ###   ###

## CT_max vs Pcrit at 12C: 7 species

bye_tips <- c("Hemilepidotus_hemilepidotus", "Blepsias_cirrhosus")
knope_ctmax <- drop.tip(knope_pruned, bye_tips)
plot.phylo(knope_ctmax)

md_pcrit_ct_max <- md_q10_ctmax_pcrit12[md_q10_ctmax_pcrit12$species!="Hemilepidotus_hemilepidotus"&
                                             md_q10_ctmax_pcrit12$species!="Blepsias_cirrhosus",
                                           c("species","avg_pcrit","ct_max_avg")]
md_pcrit_ct_max <- as.data.frame(md_pcrit_ct_max)
md_pcrit_ct_max_rn <- column_to_rownames(md_pcrit_ct_max, var = "species")

model_pcrit_ctmax_data <- comparative.data(knope_ctmax,
                                           md_pcrit_ct_max,
                                           names.col = "species",
                                           force.root = TRUE,
                                           vcv = TRUE,
                                           warn.dropped = TRUE)
model_pgls_pcrit_ctmax <- pgls(ct_max_avg~avg_pcrit, model_pcrit_ctmax_data, lambda = "ML")
model_profile_pgls_pcrit_ctmax <- pgls.profile(model_pgls_pcrit_ctmax, "lambda", N=100, param.CI = 0.95)
plot(model_pgls_pcrit_ctmax)
summary(model_pgls_pcrit_ctmax)
lm.lk_pcrit_ctmax<-pgls.profile(model_pgls_pcrit_ctmax, which="lambda")
plot(lm.lk_pcrit_ctmax)# Put these in a slide after end of SICB talk

res_pgls_pcrit_ctmax <- residuals.pgls(model_pgls_pcrit_ctmax)
plot(model_pgls_pcrit_ctmax)
qqnorm(res_pgls_pcrit_ctmax)
qqline(res_pgls_pcrit_ctmax)
#plot(res_pgls_q10)
shapiro.test(res_pgls_pcrit_ctmax) # Passed test of normality of model residuals
#anova.pgls(model_pgls_q10_pcrit_depth)

model_ols_pcrit_ctmax <- gls(ct_max_avg ~ avg_pcrit, data=md_pcrit_ct_max_rn, method="ML")
summary(model_ols_pcrit_ctmax)
res_ols_pcrit_ctmax <- residuals(model_ols_pcrit_ctmax)
plot(model_ols_pcrit_ctmax)
qqnorm(res_ols_pcrit_ctmax)
qqline(res_ols_pcrit_ctmax)
plot(res_ols_pcrit_ctmax)
shapiro.test(res_ols_pcrit_ctmax) # Passed test of normality of model residuals
anova(model_ols_pcrit_ctmax)

x_min_ctmax_pcrit <- min(md_pcrit_ct_max$avg_pcrit)
x_max_ctmax_pcrit <- max(md_pcrit_ct_max$avg_pcrit)
y_min_ctmax_pcrit <- 28.784377 + -0.078107*(min(md_pcrit_ct_max$avg_pcrit))
y_max_ctmax_pcrit <- 28.784377 + -0.078107*(max(md_pcrit_ct_max$avg_pcrit))

md_pcrit12sem <- md_all_temps[md_all_temps$temp==12,
                              c("species","sem_pcrit")]
md_pcrit12sem_Vctmax <- md_pcrit12sem[md_pcrit12sem$species!="Blepsias_cirrhosus"&
                                        md_pcrit12sem$species!="Hemilepidotus_hemilepidotus",]
md_q10_ctmax_pcrit12_7spps <- md_q10_ctmax_pcrit12[md_q10_ctmax_pcrit12$species!="Blepsias_cirrhosus"&
                                                     md_q10_ctmax_pcrit12$species!="Hemilepidotus_hemilepidotus",]
md_ctmaxsem_pcrit12sem <- full_join(md_q10_ctmax_pcrit12_7spps, md_pcrit12sem_Vctmax, by="species")


# #  
  # PLOT : CTmax vs Pcrit at 12C

ctmax_pcrit_plot <- ggplot(md_ctmaxsem_pcrit12sem, aes(x=avg_pcrit, y=ct_max_avg)) +
  geom_errorbarh(mapping = aes(xmax = avg_pcrit+sem_pcrit, xmin = avg_pcrit-sem_pcrit, height = 0.1), size=1) +
  geom_errorbar(mapping = aes(x=avg_pcrit, ymin=ct_max_avg-ct_max_sem, ymax=ct_max_avg+ct_max_sem), width=0.2, size=1) +
  #geom_line(position=pd, size=1.25)  +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("P"["crit"]~("Torr")),
       y = expression(paste("CT"["max"], " (",degree~C,")"))) #+
  #geom_segment(aes(x=x_min_ctmax_pcrit, y=y_min_ctmax_pcrit, xend=x_max_ctmax_pcrit, yend=y_max_ctmax_pcrit),
  #             color = "red", size=1)
#geom_abline(slope = -0.0327242, intercept = 3.1143425, 
#            color="red", size=1) #+
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(10,60))
#expression(paste("Temperature (",degree~C,")"))

###   ###   ###   ###   ###   ###   ###   ###
###   ###   ###   ###   ###   ###   ###

## CT_max vs Q10-Pcrit: 7 species

#bye_tips <- c("Hemilepidotus_hemilepidotus", "Blepsias_cirrhosus")
#knope_ctmax <- drop.tip(knope_pruned, bye_tips)
#plot.phylo(knope_ctmax)

md_q10pcrit_ct_max <- md_q10_ctmax_pcrit12[md_q10_ctmax_pcrit12$species!="Hemilepidotus_hemilepidotus"&
                                          md_q10_ctmax_pcrit12$species!="Blepsias_cirrhosus",
                                        c("species","q10_pcrit","ct_max_avg")]
md_q10pcrit_ct_max <- as.data.frame(md_q10pcrit_ct_max)
md_q10pcrit_ct_max_rn <- column_to_rownames(md_q10pcrit_ct_max, var = "species")

model_q10pcrit_ctmax_data <- comparative.data(knope_ctmax,
                                              md_q10pcrit_ct_max,
                                              names.col = "species",
                                              force.root = TRUE,
                                              vcv = TRUE,
                                              warn.dropped = TRUE)
model_pgls_q10pcrit_ctmax <- pgls(ct_max_avg~q10_pcrit, model_q10pcrit_ctmax_data, lambda = "ML")
model_profile_pgls_pcrit_ctmax <- pgls.profile(model_pgls_q10pcrit_ctmax, "lambda", N=100, param.CI = 0.95)
plot(model_pgls_q10pcrit_ctmax)
summary(model_pgls_q10pcrit_ctmax)
lm.lk_q10pcrit_ctmax<-pgls.profile(model_pgls_q10pcrit_ctmax, which="lambda")
plot(lm.lk_q10pcrit_ctmax)# Put these in a slide after end of SICB talk

res_pgls_q10pcrit_ctmax <- residuals.pgls(model_pgls_q10pcrit_ctmax)
plot(model_pgls_q10pcrit_ctmax)
qqnorm(res_pgls_q10pcrit_ctmax)
qqline(res_pgls_q10pcrit_ctmax)
#plot(res_pgls_q10)
shapiro.test(res_pgls_q10pcrit_ctmax) # Passed test of normality of model residuals
#anova.pgls(model_pgls_q10_pcrit_depth)

model_ols_q10pcrit_ctmax <- gls(ct_max_avg ~ q10_pcrit, data=md_q10pcrit_ct_max_rn, method="ML")
summary(model_ols_q10pcrit_ctmax)
res_ols_q10pcrit_ctmax <- residuals(model_ols_q10pcrit_ctmax)
plot(model_ols_q10pcrit_ctmax)
qqnorm(res_ols_q10pcrit_ctmax)
qqline(res_ols_q10pcrit_ctmax)
plot(res_ols_q10pcrit_ctmax)
shapiro.test(res_ols_q10pcrit_ctmax) # Passed test of normality of model residuals
anova(model_ols_q10pcrit_ctmax)

x_min_ctmax_q10pcrit <- min(md_q10pcrit_ct_max$q10_pcrit)
x_max_ctmax_q10pcrit <- max(md_q10pcrit_ct_max$q10_pcrit)
y_min_ctmax_q10pcrit <- 26.712788 + -0.143574*(min(md_q10pcrit_ct_max$q10_pcrit))
y_max_ctmax_q10pcrit <- 26.712788 + -0.143574*(max(md_q10pcrit_ct_max$q10_pcrit))

md_pcrit12sem <- md_all_temps[md_all_temps$temp==12,
                              c("species","sem_pcrit")]
md_pcrit12sem_Vctmax <- md_pcrit12sem[md_pcrit12sem$species!="Blepsias_cirrhosus"&
                                        md_pcrit12sem$species!="Hemilepidotus_hemilepidotus",]
md_q10_ctmax_pcrit12_7spps <- md_q10_ctmax_pcrit12[md_q10_ctmax_pcrit12$species!="Blepsias_cirrhosus"&
                                                     md_q10_ctmax_pcrit12$species!="Hemilepidotus_hemilepidotus",]
md_ctmaxsem_pcrit12sem <- full_join(md_q10_ctmax_pcrit12_7spps, md_pcrit12sem_Vctmax, by="species")


# #  
# PLOT : CTmax vs Q10_Pcrit

ctmax_q10pcrit_plot <- ggplot(md_ctmaxsem_pcrit12sem, aes(x=q10_pcrit, y=ct_max_avg)) +
  #geom_errorbarh(mapping = aes(xmax = avg_pcrit+sem_pcrit, xmin = avg_pcrit-sem_pcrit, height = 0.1), size=1) +
  geom_errorbar(mapping = aes(x=q10_pcrit, ymin=ct_max_avg-ct_max_sem, ymax=ct_max_avg+ct_max_sem), width=0.05, size=1) +
  #geom_line(position=pd, size=1.25)  +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Q"[10]~"of"~"P"["crit"]), #"Q"[10]~"of"~"P"["crit"]
       y = expression(paste("CT"["max"], " (",degree~C,")"))) #+
  #geom_segment(aes(x=x_min_ctmax_pcrit, y=y_min_ctmax_pcrit, xend=x_max_ctmax_pcrit, yend=y_max_ctmax_pcrit),
  #             color = "red", size=1)
#geom_abline(slope = -0.0327242, intercept = 3.1143425, 
#            color="red", size=1) #+
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(10,60))
#expression(paste("Temperature (",degree~C,")"))