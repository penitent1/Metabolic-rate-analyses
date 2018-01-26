### Number of individuals used per species for Pcrit-Mo2-temp study
### 2017

library(tidyverse)

############################################
## Read sculpin pcrit and meta-data into R
############################################

md <- read.csv("SculpinPcritData_ComparativeAnalysisFormat.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
md <- as_tibble(md)

## Pcrit and temp data WITH NA's, only "trial.no == 1" inlcuded
md_groupid <- md %>%
  filter(trial.no == 1) %>%
  #filter(!is.na(fish.id)) %>%
  group_by(species, fish.id) %>%
  summarize(number_fish = length(spps))

length(md_groupid$fish.id[md_groupid$species=="Artedius_fenestralis"]) ## 14
length(md_groupid$fish.id[md_groupid$species=="Artedius_harringtoni"]) ## 6
length(md_groupid$fish.id[md_groupid$species=="Artedius_lateralis"]) ## 13
length(md_groupid$fish.id[md_groupid$species=="Oligocottus_maculosus"]) ## 18
length(md_groupid$fish.id[md_groupid$species=="Blepsias_cirrhosus"]) ## 6
length(md_groupid$fish.id[md_groupid$species=="Clinocottus_globiceps"]) ## 16
length(md_groupid$fish.id[md_groupid$species=="Scorpaenichthys_marmoratus"]) ## 6
length(md_groupid$fish.id[md_groupid$species=="Enophrys_bison"]) ## 8
length(md_groupid$fish.id[md_groupid$species=="Hemilepidotus_hemilepidotus"]) ## 6

species <- c("Artedius fenestralis",
             "Artedius harringtoni",
             "Artedius lateralis",
             "Oligocottus maculosus",
             "Blepsias cirrhosus",
             "Clinocottus globiceps",
             "Scorpaenichthys marmoratus",
             "Enophrys bison",
             "Hemilepidotus hemilepidotus")

number_fish <- c(14,
                6,
                13,
                18,
                6,
                16,
                6,
                8,
                6)

number_fish_df <- data.frame(species, number_fish)
number_fish_df <- number_fish_df[order(number_fish_df$species),]
write.csv(number_fish_df, "number_fish_used_sculpins_2017_Somo.csv")