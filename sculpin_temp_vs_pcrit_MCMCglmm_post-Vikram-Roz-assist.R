library(ctv)
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
library(visreg)
library(RColorBrewer)
library(caper)

#####################################################################
## Read sculpin phyologeny into R and prune for species in analysis
#####################################################################
ramon_tree <- read.nexus("T5000_FigTree_TimeScaledRoot1.nex")

  ## NOTE! This does not include Blepsias cirrhosus (not in Ramon tree)
keepers_ramon <- c("Oligocottus_maculosus",
             "Clinocottus_globiceps",
             "Artedius_harringtoni",
             "Artedius_lateralis",
             "Artedius_fenestaolis", ## Ramon tree has typo here, should be "fenestralis"
             "Scorpaenichthys_marmoratus",
             "Enophrys_bison",
             "Hemilepidotus_hemilepidotus")

## Prune phylogeny for species in analysis
ramon_phy <- drop.tip(ramon_tree, setdiff(ramon_tree$tip.label, keepers_ramon))
## Fix "fenestralis" typo
ramon_phy$tip.label[ramon_phy$tip.label == "Artedius_fenestaolis"] <- "Artedius_fenestralis"
## Compute branch lengths, tree didn't have any
ramon_phy <- compute.brlen(ramon_phy, method = "Grafen")

plot.phylo(ramon_phy)
nodelabels()
tiplabels()
edgelabels()

## variance-covariance matrix for ramon tree pruned for species in study
vcv.phylo(ramon_phy, cor = TRUE)

############################################
## Read sculpin pcrit and meta-data into R
############################################

md <- read.csv("SculpinPcritData_ComparativeAnalysisFormat.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
md <- as_tibble(md)

## Pcrit and temp data, NA's filtered out, only "trial.no == 1" inlcuded
md_pcrit_smr <- md %>%
  #select(species, spps, date, temp, fish.id, mass.g, smr.best, pcrit.r, trial.no) %>%
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

## MO2:Pcrit ratios for comparing Pcrits against MO2
  ## Think carefully about what this means!
md_ratios <- md_all_temps %>%
  dplyr::select(temp, species, avg_pcrit, avg_smr) %>% ## use dplyr::select to avoid clash between MASS and dplyr packages
  mutate(ratios = avg_smr/avg_pcrit)

md_ratios$spps_names <- species_names

  ## Plot of MO2:pcrit ratios against temperatures
ratios_plot <-  ggplot(md_ratios, aes(x=temp, y=ratios, colour = species)) +
  geom_point() +
  geom_line()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # Plot of Pcrit vs temp for each species mean +/- sem

species_names <- c(
  "Artedius fenestralis",
  "Artedius harringtoni",
  "Artedius lateralis",
  "Blepsias cirrhosus",
  "Clinocottus globiceps",
  "Enophrys bison",
  "Hemilepidotus hemilepidotus",
  "Oligocottus maculosus",
  "Scorpaenichthys marmoratus"
)

md_all_temps$spps_names <- species_names
#md_lognat$spps_names <- species_names

pd <- position_dodge(0.2)

### Max depth data frame
max_depth_data <- read.csv("ct_max_maxdepth_data.csv", stringsAsFactors = FALSE,
                           strip.white = TRUE, na.strings = c("NA","."))

max_depth_df <- as.tibble(max_depth_data)
max_depth_df <- max_depth_df[max_depth_df$species!="Rhamphocottus_richardsoni",
                             c("species","spps","max_depth_mandic","center_habitat")]
max_depth_df <- max_depth_df[order(max_depth_df$max_depth_mandic),]
max_depth_df$rank_depth <- c(1,2,3,4,5,6,7,8,9)

md_depth_data <- left_join(md_all_temps, max_depth_df, key = species)

# # # 
# Colored by species ID: Pcrit
pcrit_vs_temp_plot <- ggplot(md_all_temps, aes(x=temp, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  geom_line(position=pd, size=1.25)  +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Species"))
pcrit_vs_temp_plot

# Colored by max depth: Pcrit
pcrit_vs_temp_plot_maxdepth <- ggplot(md_depth_data, aes(x=temp, y=avg_pcrit, group=species, colour = rank_depth)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  geom_line(position=pd, size=1.25) +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.title = element_text(size = 20),#legend.position = "none",
        legend.text = element_text(size = 20),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Habitat depth rank")) +
  scale_colour_gradient(low = "#56B1F7", high = "#132B43")# +
  #guides(fill = guide_legend(title = "Max depth (m)"))

md_depth_data_olma <- md_depth_data %>%
  filter(species == "Oligocottus_maculosus")
md_depth_data_enbi <- md_depth_data %>%
  filter(species == "Enophrys_bison")
md_depth_data_enbi_olma <- bind_rows(md_depth_data_olma,md_depth_data_enbi)

# Colored by max depth: Pcrit Tidepool vs Buffalo
pcrit_vs_temp_plot_maxdepth_enbi_olma <- ggplot(md_depth_data_enbi_olma, aes(x=temp, y=avg_pcrit, group=species, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  geom_line(position=pd, size=1.25) +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.title = element_text(size = 20),#legend.position = "none",
        legend.text = element_text(size = 20),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression("P"["crit"]~("Torr")),
       legend = expression("Species")) +
  coord_cartesian(ylim = c(0,160))#+#,
       #colour = expression("Habitat depth rank")) +
  #scale_colour_gradient(low = "#56B1F7", high = "#132B43")# +
#guides(fill = guide_legend(title = "Max depth (m)"))

# Colored by species ID: SMR
smr_vs_temp_plot_noleg <- ggplot(md_all_temps, aes(x=temp, y=avg_smr, colour = species)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_smr-sem_smr), ymax=(avg_smr+sem_smr)), position=pd, width=1, size=1.5) +
  geom_line(position=pd, size=1.25) +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")))

# Colored by species max depth: SMR
smr_vs_temp_plot_maxdepth <- ggplot(md_depth_data, aes(x=temp, y=avg_smr, group=species, colour = rank_depth)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_smr-sem_smr), ymax=(avg_smr+sem_smr)), position=pd, width=1, size=1.5) +
  geom_line(position=pd, size=1.25) +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.title = element_text(size = 20),#legend.position = "none",
        legend.text = element_text(size = 20),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Habitat depth rank")) +
  scale_colour_gradient(low = "#56B1F7", high = "#132B43")

# Colored by species ID: SMR vs Pcrit
smr_vs_pcrit_plot <- ggplot(md_all_temps, aes(x=avg_smr, y=avg_pcrit, colour = species)) +
  geom_errorbar(mapping = aes(ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 2.5), size=1.5) + 
  geom_line(size=1.25) +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(y = expression("P"["crit"]~("Torr")),
       x = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")))

# Colored by species ID with LEGEND: SMR vs Pcrit
smr_vs_pcrit_plot <- ggplot(md_all_temps, aes(x=avg_smr, y=avg_pcrit, colour = species)) +
  geom_errorbar(mapping = aes(ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 2.5), size=1.5) + 
  geom_line(size=1.25) +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(y = expression("P"["crit"]~("Torr")),
       x = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")))
       
# Colored by species max depth: SMR vs Pcrit
smr_vs_pcrit_plot <- ggplot(md_depth_data, aes(x=avg_smr, y=avg_pcrit, group = species, colour = rank_depth)) +
  geom_errorbar(mapping = aes(ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 2.5), size=1.5) + 
  geom_line(size=1.25) +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(y = expression("P"["crit"]~("Torr")),
       x = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Max depth (m)")) +
  scale_colour_gradient(low = "#56B1F7", high = "#132B43")

# Colored by species habitat location: SMR vs Pcrit
smr_vs_pcrit_plot <- ggplot(md_depth_data, aes(x=avg_smr, y=avg_pcrit, group = species, colour = center_habitat)) +
  geom_errorbar(mapping = aes(ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 2.5), size=1.5) + 
  geom_line(size=1.25) +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(y = expression("P"["crit"]~("Torr")),
       x = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Habitat location")) #+
  #scale_colour_gradient(low = "#56B1F7", high = "#132B43")


# Colored by species ID: Q10_SMR vs Q10_Pcrit
q10smr_vs_q10pcrit_plot <- ggplot(md_q10, aes(x=q10_smr, y=q10_pcrit, colour = species)) +
  #geom_errorbar(mapping = aes(ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.5) +
  #geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 2.5), size=1.5) + 
  #geom_line(size=1.25) +
  geom_point(size = 5) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(y = expression("Q"[10]~"of"~"P"["crit"]),
       x = expression("Q"[10]~"of"~"Mo"[2~"min"]))



pcrit_vs_temp_plot_noleg_lognat <- ggplot(md_lognat, aes(x=temp, y=lognat_pcrit, colour = species)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(lognat_pcrit-sem_ln_pcrit), ymax=(lognat_pcrit+sem_ln_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd) +
  geom_smooth(method = lm, se=FALSE) +
  geom_point(position=pd, size = 3.5) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16)) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = "ln (Pcrit)")

pcrit_vs_temp_plot_leg <- ggplot(md_all_temps, aes(x=temp, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 3.5) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16)) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = "Pcrit (Torr)")

pcrit_vs_temp_plot_leg_lognat <- ggplot(md_lognat, aes(x=temp, y=lognat_pcrit, colour = spps_names)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(lognat_pcrit-sem_ln_pcrit), ymax=(lognat_pcrit+sem_ln_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd) +
  geom_smooth(method = lm, se=FALSE) +
  geom_point(position=pd, size = 3.5) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16)) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = "ln (Pcrit)")

pcrit_vs_temp_plot_leg_lines <- ggplot(md_lognat, aes(x=temp, y=avg_pcrit, colour = spps_names)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(lognat_pcrit-sem_ln_pcrit), ymax=(lognat_pcrit+sem_ln_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd) +
  geom_smooth(method = lm, se=FALSE) +
  geom_point(position=pd, size = 3.5) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16)) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = "Pcrit (Torr))")

## playing with idea of plotting with color determined by CTmax, or 'MaxDepth'
pcrit_vs_temp_plot_noleg_sppsDepthColor <- ggplot(md_all_temps, aes(x=temp, y=avg_pcrit, group = species, colour = )) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(lognat_pcrit-sem_ln_pcrit), ymax=(lognat_pcrit+sem_ln_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd) +
  geom_smooth(method = lm, se=FALSE) +
  geom_point(position=pd, size = 3.5) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16)) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = "ln (Pcrit)")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Averages Pcrit, smr per species per temperature, no data at 16C (prep for q10)
md_temp_groups_q10 <- md_pcrit_smr %>%
  #select(species, spps, date, temp, fish.id, mass.g, smr.best, pcrit.r, trial.no) %>%
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
ct_max_df <- read.csv("ct_max_data.csv", stringsAsFactors = FALSE,
                      strip.white = TRUE, na.strings = c("NA","."))
 #md_q10_ctmax <- full_join(md_q10, ct_max_df, by = "species")
 #md_pcrit12 <- md_temp_groups_q10[md_temp_groups_q10$temp==12,]

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

## Pcrit at 12C vs max depth

md_pcrit12_maxd <- md_q10_ctmax_pcrit12[md_q10_ctmax_pcrit12$species!="Rhamphocottus_richardsoni"&
                                         md_q10_ctmax_pcrit12$species!="Blepsias_cirrhosus",
                                       c("species","q10_pcrit","q10_smr","avg_pcrit","max_depth_mandic")]

md_pcrit12_maxd$max_depth_mandic[md_pcrit12_maxd$species=="Hemilepidotus_hemilepidotus"] <- 450
md_pcrit12_maxd_rn <- column_to_rownames(md_pcrit12_maxd, var = "species")

## Regular OLS
ols_pcrit_v_maxd <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_rn, method="ML")

## PGLS assuming a Brownian correlation
pgls_bm_pcrit_v_maxd <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_rn, correlation = corBrownian(phy=ramon_phy), method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_lambda_pcrit_v_maxd <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_rn, 
                                    correlation = corPagel(value=1, phy=ramon_phy, fixed=FALSE), method="ML")
pgls_lambda_pcrit_v_maxd

## Compare models using AIC (lower is beter)
AIC(ols_pcrit_v_maxd)
AIC(pgls_bm_pcrit_v_maxd)
AIC(pgls_lambda_pcrit_v_maxd)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_lambda_pcrit_v_maxd, ols_pcrit_v_maxd)
anova(pgls_bm_pcrit_v_maxd, ols_pcrit_v_maxd)

## Is lambda significantly better than Brownian motion
anova(pgls_bm_pcrit_v_q10pcrit, pgls_lambda_pcrit_v_q10pcrit)

res <- residuals.lm(pgls_lambda_pcrit_v_q10pcrit)
plot(pgls_lambda_pcrit_v_maxd)
hist(res)
qqnorm(pgls_lambda_pcrit_v_q10pcrit)
qqline(pgls_lambda_pcrit_v_q10pcrit)

#plot(md_q10$q10_pcrit~md_q10$q10_smr)
#abline(pgls_lambda)
anova(pgls_lambda_pcrit_v_q10pcrit)

## Pcrit at 12C vs max depth --- NO HEHE!

md_pcrit12_maxd_nohehe <- md_pcrit12_maxd[md_pcrit12_maxd$species!="Hemilepidotus_hemilepidotus",]

md_pcrit12_maxd_nohehe_rn <- column_to_rownames(md_pcrit12_maxd_nohehe, var = "species")

## Regular OLS
ols_pcrit_v_maxd_nohehe <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_nohehe_rn, method="ML")

## PGLS assuming a Brownian correlation
pgls_bm_pcrit_v_maxd_nohehe <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_nohehe_rn, correlation = corBrownian(phy=ramon_phy_ctmax), method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_lambda_pcrit_v_maxd_nohehe <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_nohehe_rn, 
                                correlation = corPagel(value=1, phy=ramon_phy_ctmax, fixed=FALSE), method="ML")
pgls_lambda_pcrit_v_maxd_nohehe
summary(pgls_lambda_pcrit_v_maxd_nohehe)
## Compare models using AIC (lower is beter)
AIC(ols_pcrit_v_maxd_nohehe)
AIC(pgls_bm_pcrit_v_maxd_nohehe)
AIC(pgls_lambda_pcrit_v_maxd_nohehe)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_lambda_pcrit_v_maxd_nohehe, ols_pcrit_v_maxd_nohehe)
anova(pgls_bm_pcrit_v_maxd_nohehe, ols_pcrit_v_maxd_nohehe)

## Is lambda significantly better than Brownian motion
anova(pgls_bm_pcrit_v_maxd_nohehe, pgls_lambda_pcrit_v_maxd_nohehe)

res <- residuals.lm(pgls_lambda_pcrit_v_q10pcrit)
plot(pgls_lambda_pcrit_v_maxd)
hist(res)
qqnorm(pgls_lambda_pcrit_v_q10pcrit)
qqline(pgls_lambda_pcrit_v_q10pcrit)

#plot(md_q10$q10_pcrit~md_q10$q10_smr)
#abline(pgls_lambda)
anova(pgls_lambda_pcrit_v_q10pcrit)

########  ########  ########

# Checking lambda likelihood landscape

#define the range for lambda and at which intervals we are sampling
lambda_1 <- seq(0, 1, length.out = 100)

#Calculate likelihood over the entire range[modified from Symonds&Blombergs tutorial!]
lik <- sapply(lambda_1, function(lambda_1) logLik(gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_nohehe_rn, 
                                                      correlation = corPagel(value=lambda_1, phy=ramon_phy_ctmax, fixed=TRUE))))
plot(lik ~ lambda_1, type = "l", main = expression(paste("Likelihood Plot for ",lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#adding our ML estimate as a red line
abline(v=pgls_lambda_pcrit_v_maxd_nohehe$modelStruct, col = "red")

sculpins_pcrit12Vdepth <- comparative.data()

## PGLS model that returned a result

comparative.data(XTBGTree, PhyloDataUN, names.col = "Species", vcv =TRUE,
                 vcv.dim=3, warn.dropped=TRUE) -> ComparedData

pgls(TLP_DS ~ AverageOP_DS, ComparedData, lambda = "ML") -> model1

pgls.profile(model1, "lambda", N =100, param.CI = 0.95) -> modelprofile


pic_pcrit_v_q10pcrit <- pic(md_pcrit12_maxd_rn$avg_pcrit, ramon_phy)
pic_pcrit_v_q10pcrit

# # # >>>>> PLOT data with both slopes

pcrit12_v_maxd_plot <- ggplot(md_pcrit12_maxd, aes(x=max_depth_mandic, y=avg_pcrit)) +
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
  geom_abline(slope = 0.03985455, intercept = 29.61927908, 
              color="black", size=1) +
  #geom_abline(slope = 0.05706547, intercept = 28.81887854,
  #            color="blue", size=1) +
  coord_cartesian(ylim = c(10,60))



##################################

## END PCRIT_12 VS MAX_DEPTH ANALYSIS

##################################
##################################
##################################

##########
##########
##########

## Q10-Pcrit at 12C vs max depth

md_pcrit12_maxd <- md_q10_ctmax_pcrit12[md_q10_ctmax_pcrit12$species!="Rhamphocottus_richardsoni"&
                                          md_q10_ctmax_pcrit12$species!="Blepsias_cirrhosus",
                                        c("species","q10_pcrit","q10_smr","avg_pcrit","max_depth_mandic")]

md_pcrit12_maxd$max_depth_mandic[md_pcrit12_maxd$species=="Hemilepidotus_hemilepidotus"] <- 450
md_pcrit12_maxd_rn <- column_to_rownames(md_pcrit12_maxd, var = "species")

## Regular OLS
ols_q10pcrit_v_maxd <- gls(q10_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_rn, method="ML")

## PGLS assuming a Brownian correlation
pgls_bm_q10pcrit_v_maxd <- gls(q10_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_rn, correlation = corBrownian(phy=ramon_phy), method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_lambda_q10pcrit_v_maxd <- gls(q10_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_rn, 
                                correlation = corPagel(value=1, phy=ramon_phy, fixed=FALSE), method="ML")
pgls_lambda_q10pcrit_v_maxd

## Compare models using AIC (lower is beter)
AIC(ols_q10pcrit_v_maxd)
AIC(pgls_bm_q10pcrit_v_maxd)
AIC(pgls_lambda_q10pcrit_v_maxd)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_lambda_q10pcrit_v_maxd, ols_q10pcrit_v_maxd)
anova(pgls_bm_q10pcrit_v_maxd, ols_q10pcrit_v_maxd)

## Is lambda significantly better than Brownian motion
anova(pgls_bm_q10pcrit_v_maxd, pgls_lambda_q10pcrit_v_maxd)

res <- residuals.lm(pgls_lambda_q10pcrit_v_maxd)
plot(pgls_lambda_q10pcrit_v_maxd)
hist(res)
qqnorm(pgls_lambda_q10pcrit_v_maxd)
qqline(pgls_lambda_q10pcrit_v_maxd)

#plot(md_q10$q10_pcrit~md_q10$q10_smr)
#abline(pgls_lambda)
anova(pgls_lambda_q10pcrit_v_maxd)

## Pcrit at 12C vs max depth --- NO HEHE!

md_pcrit12_maxd_nohehe <- md_pcrit12_maxd[md_pcrit12_maxd$species!="Hemilepidotus_hemilepidotus",]

md_pcrit12_maxd_nohehe_rn <- column_to_rownames(md_pcrit12_maxd_nohehe, var = "species")

## Regular OLS
ols_pcrit_v_maxd_nohehe <- gls(q10_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_nohehe_rn, method="ML")

## PGLS assuming a Brownian correlation
pgls_bm_pcrit_v_maxd_nohehe <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_nohehe_rn, correlation = corBrownian(phy=ramon_phy_ctmax), method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_lambda_pcrit_v_maxd_nohehe <- gls(avg_pcrit ~ max_depth_mandic, data=md_pcrit12_maxd_nohehe_rn, 
                                       correlation = corPagel(value=1, phy=ramon_phy_ctmax, fixed=FALSE), method="ML")
pgls_lambda_pcrit_v_maxd_nohehe
summary(pgls_lambda_pcrit_v_maxd_nohehe)
## Compare models using AIC (lower is beter)
AIC(ols_pcrit_v_maxd_nohehe)
AIC(pgls_bm_pcrit_v_maxd_nohehe)
AIC(pgls_lambda_pcrit_v_maxd_nohehe)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_lambda_pcrit_v_maxd_nohehe, ols_pcrit_v_maxd_nohehe)
anova(pgls_bm_pcrit_v_maxd_nohehe, ols_pcrit_v_maxd_nohehe)

## Is lambda significantly better than Brownian motion
anova(pgls_bm_pcrit_v_maxd_nohehe, pgls_lambda_pcrit_v_maxd_nohehe)

res <- residuals.lm(pgls_lambda_pcrit_v_q10pcrit)
plot(pgls_lambda_pcrit_v_maxd)
hist(res)
qqnorm(pgls_lambda_pcrit_v_q10pcrit)
qqline(pgls_lambda_pcrit_v_q10pcrit)

#plot(md_q10$q10_pcrit~md_q10$q10_smr)
#abline(pgls_lambda)
anova(pgls_lambda_pcrit_v_q10pcrit)


##################################

## END Q10_PCRIT_12 VS MAX_DEPTH ANALYSIS

##################################
##################################
##################################



########################################
########################################
########################################

##########
##########
##########

## Pcrit at 12C vs q10 for Pcrit

## Take out Grunt sculpin data
md_q10_pcrit12 <- md_q10_ctmax_pcrit12[md_q10_ctmax_pcrit12$species!="Rhamphocottus_richardsoni"&
                                         md_q10_ctmax_pcrit12$species!="Blepsias_cirrhosus",
                                     c("species","q10_pcrit","q10_smr","avg_pcrit")]

md_q10_pcrit12_rn <- column_to_rownames(md_q10_pcrit12, var = "species")

## Regular OLS
ols_pcrit_v_q10pcrit <- gls(q10_pcrit ~ avg_pcrit, data=md_q10_pcrit12_rn, method="ML")

## PGLS assuming a Brownian correlation
pgls_bm_pcrit_v_q10pcrit <- gls(q10_pcrit ~ avg_pcrit, data=md_q10_pcrit12_rn, correlation = corBrownian(phy=ramon_phy), method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_lambda_pcrit_v_q10pcrit <- gls(q10_pcrit ~ avg_pcrit, data=md_q10_pcrit12_rn, 
                   correlation = corPagel(value=1, phy=ramon_phy, fixed=FALSE), method="ML")
pgls_lambda_pcrit_v_q10pcrit

## Compare models using AIC (lower is beter)
AIC(ols_pcrit_v_q10pcrit)
AIC(pgls_bm_pcrit_v_q10pcrit)
AIC(pgls_lambda_pcrit_v_q10pcrit)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_lambda_pcrit_v_q10pcrit, ols_pcrit_v_q10pcrit)
anova(pgls_bm_pcrit_v_q10pcrit, ols_pcrit_v_q10pcrit)

## Is lambda significantly better than Brownian motion
anova(pgls_bm_pcrit_v_q10pcrit, pgls_lambda_pcrit_v_q10pcrit)

#check for departures from normal distribution of residuals
res <- resid(pgls_lambda_pcrit_v_q10pcrit, type="n")
qqnorm(res, col="blue")
qqline(res, col="blue")

res <- residuals.lm(pgls_lambda_pcrit_v_q10pcrit)
plot(pgls_lambda_pcrit_v_q10pcrit)
summary(pgls_lambda_pcrit_v_q10pcrit)
#hist(res)
#qqnorm(pgls_lambda_pcrit_v_q10pcrit)
#qqline(pgls_lambda_pcrit_v_q10pcrit)

#plot(md_q10$q10_pcrit~md_q10$q10_smr)
#abline(pgls_lambda)
anova(pgls_lambda_pcrit_v_q10pcrit)

# # # >>>>> PLOT for Q10 vs Pcrit at 12C
md_q10_pcrit12_plotting <- md_all_temps[md_all_temps$temp==12,
                                        c("species","avg_pcrit","sem_pcrit")]
md_q10_pcrit12_plotting <- left_join(md_q10_pcrit12_plotting, md_q10_ctmax, by="species")

summary(pgls_lambda_pcrit_v_q10pcrit)
vcv.phylo(ramon_phy)
max(nodeHeights(ramon_phy))

pcrit12_v_q10pcrit_plot <- ggplot(md_q10_pcrit12_plotting, aes(x=avg_pcrit, y=q10_pcrit)) +
  geom_errorbarh(mapping = aes(xmax = avg_pcrit+sem_pcrit, xmin = avg_pcrit-sem_pcrit, height = 0.2), size=1.5) +
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
  geom_abline(slope = -0.032497, intercept = 3.177915, 
              color="black", size=1) +
  #geom_abline(slope = 0.05706547, intercept = 28.81887854,
  #            color="blue", size=1) +
  coord_cartesian(ylim = c(10,60))



##################################

## END PCRIT VS Q10-PCRIT ANALYSIS

##################################
##################################
##################################


## Data pruned for species in Q10_pcrit/Q10_smr analysis
md_pcrit <- md_q10_ctmax[md_q10_ctmax$species!="Blepsias_cirrhosus"&md_q10_ctmax$species!="Rhamphocottus_richardsoni",1:3]
  
  ## Using gls to test for an association of body size with hind limb (dumb, i know but bear with me)
  
## Models of correlations between Q10-pcrit and Q10-SMR
md_pcrit_rn <- column_to_rownames(md_pcrit, var = "species")

  ## Regular OLS
  ols <- gls(q10_pcrit ~ q10_smr, data=md_pcrit_rn, method="ML")
  
  ## PGLS assuming a Brownian correlation
  pgls_bm <- gls(q10_pcrit ~ q10_smr, data=md_pcrit_rn, correlation = corBrownian(phy=ramon_phy), method="ML")
  
  ## PGLS assuming a Pagel's lambda correlation
  pgls_lambda <- gls(q10_pcrit ~ q10_smr, data=md_pcrit_rn, 
                     correlation = corPagel(value=1, phy=ramon_phy, fixed=FALSE), method="ML")
  
  pgls_lambda
  
  ## Compare models using AIC (lower is beter)
  AIC(ols)
  AIC(pgls_bm)
  AIC(pgls_lambda)
  
  ## Using a likelihood ratio test to evaluate whether pgls is justified
  anova(pgls_lambda, ols)
  
  ## Is lambda significantly better than Brownian motion
  anova(pgls_bm, pgls_lambda)

res <- residuals.lm(pgls_lambda)
plot(pgls_lambda)
hist(res)
qqnorm(pgls_lambda)
qqline(pgls_lambda)

plot(md_q10$q10_pcrit~md_q10$q10_smr)
abline(pgls_lambda)
anova(pgls_lambda)

###########################################
###########################################
###########################################
###########################################

###### CTmax versus Pcrit analyses

###########################################
###########################################
###########################################
###########################################

ramon_phy_ctmax <- drop.tip(ramon_phy, "Hemilepidotus_hemilepidotus")

md_ctmax <- md_q10_ctmax[md_q10_ctmax$species!="Hemilepidotus_hemilepidotus"&
                           md_q10_ctmax$species!="Blepsias_cirrhosus"&
                           md_q10_ctmax$species!="Rhamphocottus_richardsoni",]
md_ctmax_rn <- column_to_rownames(md_ctmax, var = "species")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # Some plots of phylogeny colored by trait values, or phenogram of trait values


md_phylo_ctmax <- md_ctmax[,c(1,5)]
md_phylo_ctmax_rn <- column_to_rownames(md_phylo_ctmax, var = "species")
md_phylo_ctmax_matrix <- as.matrix(md_phylo_ctmax_rn)[,1]

#############
# # # Data for ct max phenogram
#############
# Using ct_max_df, no q10 values or anything besides ct_max data
md_phylo_ctmax_only <- ct_max_df[ct_max_df$spps!="rhri",c(1,3)]
md_phylo_ctmax_only_rn <- column_to_rownames(md_phylo_ctmax_only, var = "species")
md_phylo_ctmax_only_matrix <- as.matrix(md_phylo_ctmax_only_rn)[,1]
#############
# # # END Data for ct max phenogram
#############

md_pcrit
md_phylo_pcrit <- md_pcrit[,c(1:2)]
md_phylo_q10_pcrit_rn <- column_to_rownames(md_phylo_pcrit, var = "species")
md_phylo_q10_pcrit_matrix <- as.matrix(md_phylo_q10_pcrit_rn)[,1]

md_phylo_pcrit_noHehe <- md_pcrit[md_pcrit$species!="Hemilepidotus_hemilepidotus",]
md_phylo_pcrit_noHehe_rn <- column_to_rownames(md_phylo_pcrit_noHehe, var = "species")
md_phylo_pcrit_noHehe_matrix <- as.matrix(md_phylo_pcrit_noHehe_rn)[,1]

## Not great visually
ct_max_phylo_plot <- contMap(ramon_phy_ctmax, md_phylo_ctmax_matrix)

##
plot(ramon_phy_ctmax, adj = 0, label.offset = 0.75)
tiplabels(pch = 21, col="blue", adj = 1, bg=ctmax_labels, cex=2)

plot(mini.phy, adj=0, label.offset=0.75) #adjust if necessary/desired!
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=2) #ditto! try different values!

## This seems more promising!
ctmax_phylo_plot <- plotBranchbyTrait(ramon_phy_ctmax,
                                      md_phylo_ctmax_matrix,
                                      mode = "tips",
                                      palette = "rainbow")

## Also potentially good
pheno_ctmax <- phenogram(ramon_phy_ctmax, md_phylo_ctmax_only_matrix,
                         ylab=expression(paste("CTmax (",degree~C,")")),
                         xlab=NULL,
                         xaxt='n')

pheno_q10_pcrit <- phenogram(ramon_phy_ctmax, md_phylo_pcrit_noHehe_matrix,
                             ylab=expression('Q'[10]:~'P'[crit]),
                             xlab=NULL)


pcrit_q10_phylo_plot <- plotBranchbyTrait(ramon_phy_ctmax, 
                                          md_phylo_ctmax_matrix,
                                          mode = "tips", 
                                          palette = "heat.colors")


plotBranchbyTrait(tree, x, mode=c("edges","tips","nodes"), palette="rainbow", 
                  legend=TRUE, xlims=NULL, ...)

############################
## Q10-Pcrit against CTmax
## Regular OLS
ols_ctmax_pcrit <- gls(q10_pcrit ~ ct_max_avg, 
                    data=md_ctmax_rn, 
                    method="ML")

## PGLS assuming a Brownian correlation
pgls_bm_ctmax_pcrit <- gls(q10_pcrit ~ ct_max_avg, 
                     data=md_ctmax_rn, 
                     correlation = corBrownian(phy=ramon_phy_ctmax), 
                     method="ML")

## PGLS assuming a corPagel vcv structure (multiplies phylo vcv by lambda, value of lambda estimated by ML)
pgls_lambda_ctmax_pcrit <- gls(q10_pcrit ~ ct_max_avg, 
                              data=md_ctmax_rn, 
                              correlation = 
                                corPagel(value=1, phy=ramon_phy_ctmax, fixed=FALSE), 
                              method="ML")

## Compare models using AIC (lower is beter)
AIC(ols_ctmax_pcrit)
AIC(pgls_bm_ctmax_pcrit)
AIC(pgls_lambda_ctmax_pcrit)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_lambda_ctmax_pcrit, ols_ctmax_pcrit)
anova(pgls_bm_ctmax_pcrit, ols_ctmax_pcrit)
## Is lambda significantly better than Brownian motion
anova(pgls_bm_ctmax_pcrit, pgls_lambda_ctmax_pcrit)

res <- residuals.lm(pgls_lambda_ctmax_pcrit)
plot(pgls_lambda_ctmax_pcrit)
hist(res)
qqnorm(pgls_lambda)
qqline(pgls_lambda)

plot(md_ctmax_rn$q10_pcrit~md_ctmax_rn$ct_max_avg)
abline(pgls_lambda_ctmax_pcrit)
anova(pgls_lambda)























##########
##########
##
## PGLS assuming a Pagel's lambda correlation
##
##########
##########

# The lambda, corPagel model with fixed = FALSE and value = 1

## Lambda value close to 0

pgls_lambda_ctmax_12.1 <- gls(pcrit_12 ~ ct_max, 
                         data=md_ctmax, 
                         correlation = 
                           corPagel(value=0.1, phy=ramon_phy_ctmax, fixed=TRUE), 
                         method="ML")
pgls_lambda_ctmax_12.2 <- gls(pcrit_12 ~ ct_max, 
                              data=md_ctmax, 
                              correlation = 
                                corPagel(value=0.5, phy=ramon_phy_ctmax, fixed=TRUE), 
                              method="ML")
pgls_lambda_ctmax_12.3 <- gls(pcrit_12 ~ ct_max, 
                              data=md_ctmax, 
                              correlation = 
                                corPagel(value=0.9, phy=ramon_phy_ctmax, fixed=TRUE), 
                              method="ML")

### Compare AIC and logLik to determine whether there actually is phylo signal
AIC(pgls_lambda_ctmax_12.1) ## AIC = 55.25429 THIS is lowest AIC and lowest logLik
AIC(pgls_lambda_ctmax_12.2) ## AIC = 56.42812
AIC(pgls_lambda_ctmax_12.3) ## AIC = 58.4175

summary(pgls_lambda_ctmax_12.1) ## logLik = -24.62714
summary(pgls_lambda_ctmax_12.2) ## logLik = -25.21406
summary(pgls_lambda_ctmax_12.3) ## logLik = -26.20875

## Compare models using AIC (lower is beter)
AIC(ols_ctmax_12)
AIC(pgls_bm_ctmax_12)
AIC(pgls_lambda_ctmax_12.1)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_bm_ctmax_12, ols_ctmax_12)
anova(pgls_lambda_ctmax_12.1, ols_ctmax_12)

## Is lambda significantly better than Brownian motion
anova(pgls_bm_ctmax_12, pgls_lambda_ctmax_12.1)

## ct_max appears not to be a significant model term
anova(ols_ctmax_12)

#################################################
#################################################
#################################################
### BLEEP BOOP BLEEP BOOP BLEEP BOOP BLEEP BOOP

## Pcrit at 16 against CTmax
## Regular OLS
ols_ctmax_16 <- gls(pcrit_16 ~ ct_max, 
                    data=md_ctmax, 
                    method="ML")

## PGLS assuming a Brownian correlation
pgls_bm_ctmax_16 <- gls(pcrit_16 ~ ct_max, 
                        data=md_ctmax, 
                        correlation = corBrownian(phy=ramon_phy_ctmax), 
                        method="ML")
##########
##########
##
## PGLS assuming a Pagel's lambda correlation
##
##########
##########

pgls_lambda_ctmax_16 <- gls(pcrit_16 ~ ct_max, 
                              data=md_ctmax, 
                              correlation = 
                                corPagel(value=1, phy=ramon_phy_ctmax, fixed=FALSE), 
                              method="ML")

## Lambda value close to 0

pgls_lambda_ctmax_16.1 <- gls(pcrit_16 ~ ct_max, 
                              data=md_ctmax, 
                              correlation = 
                                corPagel(value=0.1, phy=ramon_phy_ctmax, fixed=TRUE), 
                              method="ML")
pgls_lambda_ctmax_16.2 <- gls(pcrit_16 ~ ct_max, 
                              data=md_ctmax, 
                              correlation = 
                                corPagel(value=0.5, phy=ramon_phy_ctmax, fixed=TRUE), 
                              method="ML")
pgls_lambda_ctmax_16.3 <- gls(pcrit_16 ~ ct_max, 
                              data=md_ctmax, 
                              correlation = 
                                corPagel(value=0.9, phy=ramon_phy_ctmax, fixed=TRUE), 
                              method="ML")

### Compare AIC and logLik to determine whether there actually is phylo signal
AIC(pgls_lambda_ctmax_16.1) ## AIC = 55.25429 THIS is lowest AIC and lowest logLik
AIC(pgls_lambda_ctmax_16.2) ## AIC = 56.42812
AIC(pgls_lambda_ctmax_16.3) ## AIC = 58.4175

summary(pgls_lambda_ctmax_16.1) ## logLik = -24.62714
summary(pgls_lambda_ctmax_16.2) ## logLik = -25.21406
summary(pgls_lambda_ctmax_16.3) ## logLik = -26.20875

## Compare models using AIC (lower is beter)
AIC(ols_ctmax_16)
AIC(pgls_bm_ctmax_16)
AIC(pgls_lambda_ctmax_16.1)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_bm_ctmax_16, ols_ctmax_16)
anova(pgls_lambda_ctmax_16.1, ols_ctmax_16)

## Is lambda significantly better than Brownian motion
anova(pgls_bm_ctmax_16, pgls_lambda_ctmax_16.1)

## ct_max appears not to be a significant model term
anova(ols_ctmax_16)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

## RANDOM EXTRA CODE

#############################################################################
#############################################################################




# Use mutate to calculate a q10 for EACH INDIVIDUAL(?) or each species?
# Try both, see if different.
#q10_pcrit <- (md_ids_q10s$pcrit.r[md_ids_q10s$temp == 20]/
#                md_ids_q10s$pcrit.r[md_ids_q10s$temp == 12])^
#  (10/(md_ids$temp[md_ids$temp == 20]-md_ids$temp[md_ids$temp == 12]))

## Mean values data
#md_means <- md_pcrit_smr %>%
#  group_by(spps, temp) %>%
#  summarise(mean_pcrit = mean(pcrit.r), mean_smr = mean(smr.best))

#write.csv(md_means, file = "sculpin_pcrit-smr_temp_means.csv")

## Q10s
#md_q10 <- md_means %>%
#  group_by(spps) %>%
#  summarise(q10_pcrit = (mean_pcrit[temp == 20]/mean_pcrit[temp == 12])/(10/(temp[temp == 20]-temp[temp == 12])),
#            q10_smr = (mean_smr[temp == 20]/mean_smr[temp == 12])/(10/(temp[temp == 20]-temp[temp == 12])))

#species <- c("Artedius_fenestralis",
#             "Artedius_harringtoni",
#             "Artedius_lateralis",
#             "Blepsias_cirrhosus",
#             "Clinocottus_globiceps",
#             "Enophrys_bison",
#             "Hemilepidotus_hemilepidotus",
#             "Oligocottus_maculosus",
#             "Scorpaenichthys_marmoratus")
#md_q10$species <- species
#md_q10 <- as.data.frame(md_q10)



#write.csv(md_q10, file="C:/Users/derek/Documents/Metabolic-rate-analyses//md_q10.csv")
#write.csv(md_q10, file = "sculpin_pcrit-smr_q10")
#md_q10 <- read.csv("md_q10.csv", row.names=1)
sculpin_td <- treedata(sculpin_phy, md_q10)
str(sculpin_td)  

pcrit_q10 <- sculpin_td$data[,"q10_pcrit"]
sculpin_tree <- sculpin_td$phy

pcrit_q10 <- pcrit_q10[sculpin_tree$tip.label]
##############################################
#states <- anolis_td$data[,"AVG_SVL"]
#tree <- anolis_td$phy

#states <- states[tree$tip.label]

## ## Estimating phylogenetic signal (Pagel's lambda)
fit_lambda <- fitContinuous(sculpin_tree, pcrit_q10, model="lambda")

## The estimate will be between 0 and 1
fit_lambda$opt$lambda


##############################################
## From Dolph's website

#pagel_model <- gls(mean_pcrit ~ temp, data=md_means, correlation = corPagel(1,sculpin_phy,fixed=FALSE))

#gls(y ~ x, data=mydata, correlation=corPagel(1,mytree, fixed = FALSE))

#write.nexus(sculpin_phy, file = "sculpin_phy_nex")
#write.tree(sculpin_phy, file = "sculpin_phy.tre")
#read.tree("sculpin_phy.tre")

#sculpin_phy <-  read.nexus("bayes_tree_final_knope.nex")
#sculpin_phy <- read.tree("sculpin_phy.tre")
#write.tree(sculpin_phy) ## prints tree in parenthetical format
#write.tree(knope_pruned)

#Aphylo <- vcv.phylo(sculpin_phy, cor = TRUE)
#invAphylo <- inverseA(Aphylo)

#sculpin_phy <- read.tree("sculpin_phy_midpoint-rooted_FigTree.tre")
#sculpin_phy_prop <- read.tree("sculpin_phy_proporNotClado_FigTree.tre")
#sculpin_phy <- read.tree("knope_rooted_ultrametric")
#sculpin_phy_add2 <- bind.tip(sculpin_phy, "Oligocottus_maculosus_20", edge.length=NULL, where=18
#         )
#bind.tip(tree, tip.label, edge.length=NULL, where=NULL, position=0,
#         interactive=FALSE, ...)

#keepers <- c("Oligocottus_maculosus",
#             "Clinocottus_globiceps",
#             "Artedius_harringtoni",
#             "Artedius_lateralis",
#             "Artedius_fenestralis",
#             "Blepsias_cirrhosus",
#             "Scorpaenichthys_marmoratus",
#             "Enophrys_bison",
#             "Hemilepidotus_hemilepidotus")

#tree_root <- reroot(sculpin_phy, node.number = 9)
#tree_root <- midpoint.root(sculpin_phy) # Roots tree but creates branch length 0
#is.rooted(sculpin_phy)
#is.ultrametric(sculpin_phy)
#plot.phylo(sculpin_phy)
#ramon_phy <- chronopl(ramon_phy, 1)

#tree_root$edge.length <- tree_root$edge.length + 0.001

#tree_root$edge.length

#sculpin_phy$edge.length <-  mytree$edge.length + 0.001

## Pcrit and temp data group_by fish.id;
# calculate q10's for each fish.id
#md_ids <- md_pcrit_smr %>%
#  select(species, spps, date, temp, fish.id, mass.g, smr.best, pcrit.r, trial.no) %>%
#  filter(temp != 16) %>%
#  group_by(fish.id, temp) %>%
#  summarise(avg_pcrit = mean(pcrit.r), avg_smr = mean(smr.best))
