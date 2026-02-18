# Alberto Rovellini
# 2/12/2026
# This script produces mFC values for Atlantis GOA projections based on historical F values from GOA assessments
# period 2016-2020 is used
# F time series were extracted from the most recent FULL assessment for each stock as of May 2025
# https://docs.google.com/spreadsheets/d/1RpWddJOC8l4aPhS85SQ615Qd3qMUhs6v/edit?usp=drive_link&ouid=105160877592505917667&rtpof=true&sd=true

# this script produces PRE-CALIBRATION mFC values
# it means that they need to undergo the 1-yr tuning, twice
# it focuses on groundfish. leaving the rest as it is
# eventually, we will need the mFC change setup to apply these F values after the burn-in
# Use the same logic used for the caps paper

library(readxl)
library(tidyverse)

# read in species info
grps <- read.csv("data/GOA_Groups.csv")
codes <- grps %>% pull(Code)

# FMP groundfish ----------------------------------------------------------

dat <- read_xlsx("data/F_from_assessments.xlsx", range = "A1:V49")

# pivot
dat_long <- dat %>%
  pivot_longer(-Year, names_to = "Code_tmp", values_to = "F")

# handle duplicate names:
# New names:
# • `FFS` -> `FFS...9`
# • `FFS` -> `FFS...10`
# • `RFS` -> `RFS...15`
# • `RFS` -> `RFS...16`

key <- data.frame("Code_tmp" = c("FFS...9", "FFS...10", "RFS...15", "RFS...16"),
                  "Species" = c("NRS", "SRS", "Northern", "REBS"),
                  "Code" = c("FFS","FFS","RFS","RFS"),
                  "Biomass_1990" = c(42569,96164,215131,43492))

dat_long1 <- dat_long %>% filter(!Code_tmp %in% key$Code_tmp) %>% rename(Code = Code_tmp)
dat_long2 <- dat_long %>% filter(Code_tmp %in% key$Code_tmp)

dat_long2 <- dat_long2 %>%
  left_join(key, by = "Code_tmp") %>%
  group_by(Year,Code) %>%
  summarise(`F` = weighted.mean(`F`,Biomass_1990))

dat_final <- rbind(dat_long1, dat_long2) %>%
  arrange(Year,Code)

# turn to mFC, but only where needed
exp_rate_stocks <- c("SKL","SKO","SKB","RFD","THO","DFD","SCU")

# Other groups ------------------------------------------------------------
# Lots of groups to fill in. The only groups without estimates for 1990s but with estimates for later decades are skates
# skates have exploitation rates
# I'd say use old 1/4 M values for consistency...

# Pacific halibut
# the target is SPR 46%
# we need to convert that to biomass as per the assessment and get the equivalent depletion rference point
# we then need to profile over F and get the F that leads us there
# the old model used a low F for HAL, 1/4 M
# F in the new model will be much higher, and as such we will need to calibrate HAL most likely
# to get the model closer to where it will be, use FMSY from the OY paper
# FMSY for HAL was 0.1151807, which corresponded to 44% depletion
# this is probably in the ballpark of F46%, so use that for calibration purposes
f_hal <- 0.1151807
mfc_hal <- 1-exp(-f_hal/365)

# Create PRM files --------------------------------------------------------
# this will be used for projection runs
oldfile <- "C:/Users/Alberto Rovellini/Documents/GOA/Paper3_projections/Atlantis_GOA_projections/AtlantisGOAV0_02357/GOA_harvest_background.prm"
ref_harvest <- readLines(oldfile)

# get harvest parameters outside loop
ref_mfc <- list()
for(sp in codes) {
  
  mfc_line <- ref_harvest[grep(paste0("mFC_", sp, " "), ref_harvest) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  ref_mfc[[sp]] <- list(mfc = mfc)
}

# this will have to be tuned with the one-yr runs
# this will be used for the projection runs
dat_recent <- dat_final %>%
  filter(Year %in% c(2016:2020)) %>%
  group_by(Code) %>%
  summarise(`F` = mean(`F`, na.rm = T))

# turn to mFC
dat_recent <- dat_recent %>%
  rowwise() %>%
  mutate(mFC = ifelse(Code %in% exp_rate_stocks,`F`/365,1-exp(-`F`/365))) %>%
  ungroup() %>%
  filter(!is.nan(mFC))

write.csv(dat_recent, "data/f_2016_2020.csv", row.names = F)

# NB: these values have to be calibrated

newfile <- paste0('mFC_tuning/mfc_recent_BEFORE_TUNING_feb2026.prm')
file.create(newfile)

for(i in 1:length(codes)){
  
  sp <- codes[i]
  parname <- paste0("mFC_", sp, " 33")
  
  print(sp)
  
  if(sp %in% dat_recent$Code){
    parval <- dat_recent %>% filter(Code == sp) %>% pull(mFC)
  } else {
    parval <- ref_mfc[[which(names(ref_mfc)==sp)]]$mfc
  }
  
  parval <- signif(parval, 8)
  
  parline <- paste(as.character(c(parval,rep(0,32))),collapse = " ")
  
  cat(parname, file=newfile, append=TRUE,'\n')
  cat(parline, file=newfile, append=TRUE, '\n')
  
}

# these values are pretty close to the 2015-2019 (pre-calibration) values