# Alberto Rovellini
# 02/17/2026
# This code builds the sets of input files needed to run the projection runs for the topic 2 paper
# These runs have several dimensions:

# 2 fishing setups: F=0 and F=recent avg
# 4 climate scenarios: 3 SSPs + 1 with fixed, historical climate (hindcast avg)
# total: 8 runs

# these dimensions will be defined by parameters in various files:

# change in f and start of HCR management at year 30 are all done in the harvest.prm files ( see older test runs to get this right )
# Climate scenarios will be determined by force.prm and will include salt, temp, hydro, and plankton scalar forcings
# runs will be 110 years: 30 years burn-in, 80 years projections
# The setup will take:
# 2 harvest files
# 4 force.prm files
# 8 sh files

library(tidyverse)

options(scipen = 999)

rm(list = ls())

grp <- read.csv("data/GOA_Groups.csv")
codes <- grp %>% pull(Code)
verts <- grp %>% filter(GroupType %in% c("MAMMAL","SHARK","BIRD","FISH")) %>% pull(Code)
hcr_grp <- c("POL", "COD", "POP", "SBF") # these are groups managed with HCR; they will need the FMSY proxy and have done the full mfc ramp
fmsy_grp <- c(hcr_grp, "HAL") # these are groups managed with HCR and halibut; they will need the FMSY proxy and have done the full mfc ramp
other_fmp_grp <- c("ATF", "FHS", "REX", "FFS", "FFD", "SKL", "SKB", "SKO", "RFS", "RFP", "RFD", "THO", "DFS", "DFD", "SCU")
oy_grp <- c(hcr_grp, other_fmp_grp) # halibut is not a OY group

# begin from the SS model that is used to run the F ramps and determine the reference points
# that model had been built on run 2357, and it uses 1990s mFC values and 1990s climatologies
ss_model <- "AtlantisGOA_SS/"

# copy to a new model
# Copy a folder and all its contents
if(!dir.exists("AtlantisGOA_Topic2")){dir.create("AtlantisGOA_Topic2/")}
file.copy(from = list.files("AtlantisGOA_SS/", full.names = TRUE), 
          to = "AtlantisGOA_Topic2/", 
          recursive = TRUE)

# remove all harvest.prm and sh files
files_to_remove <- list.files("AtlantisGOA_Topic2/", 
                              pattern = "run_atlantis|GOA_harvest", 
                              full.names = TRUE, 
                              recursive = TRUE)
file.remove(files_to_remove)

# Harvest.prm -------------------------------------------------------------

# create template harvest.prm with correct parameters, mfc modifiers, etc.
# the start point is the same harvest.prm file used to create the set of mfc ramps used in the SS runs
file.copy("AtlantisGOAV0_02357/GOA_harvest_background.prm", "AtlantisGOA_Topic2/GOA_harvest_2357.prm")
harvest_base <- readLines("AtlantisGOA_Topic2/GOA_harvest_2357.prm")

# need to do flagmfc change for everything that does not have the HCR
# the burn-in uses tuned mFC values hitting mean F from the 1990s
# the projection needs to use tuned mFC values hitting mean F from 2016-2020
# the run that has those properties is 2364
# the mFC in 2357 (and the PRM currently being modified) is the base value applied in the burn-in
# the mFC in 2364 is the target mFC to apply in proj
# set up a multiplier for the projection
# set all the necessary flags and switches for the mFC change

file.copy("AtlantisGOAV0_02364/GOA_harvest_background.prm", "AtlantisGOA_Topic2/GOA_harvest_2364.prm")
harvest_target <- readLines("AtlantisGOA_Topic2/GOA_harvest_2364.prm")

# get mFC from current template PRM (the 1990s mFC to scale to proj values)
burnin_mfc <- list()
for(sp in codes) {
  
  mfc_line <- harvest_base[grep(paste0("mFC_", sp, " "), harvest_base) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  burnin_mfc[[sp]] <- data.frame("Code" = sp, "burnin_mfc"=  mfc)
}
burnin_mfc <- bind_rows(burnin_mfc)

# get mFC from target harvest (the 2015-2019 mFC)
# NB: Halibut is going to need a different approach. The target mFC will be the mFC associated with the FMSY run
# NB2: the HCR species will not be using these values in projection because the HCR handles F dynamically, but let's still set it up in case we need to do runs without the HCR
proj_mfc <- list()
for(sp in codes) {
  
  mfc_line <- harvest_target[grep(paste0("mFC_", sp, " "), harvest_target) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  proj_mfc[[sp]] <- data.frame("Code" = sp, "proj_mfc"=  mfc)
}
proj_mfc <- bind_rows(proj_mfc)

# plug in target value for HAL, which won't be using the HCR
# NOT DOING THIS FOR TOPIC 2
# proj_mfc[proj_mfc$Code == "HAL",]$proj_mfc <- ref_points %>% filter(Code == "HAL") %>% pull(fref_mfc_for_check)

# get scalars
scalars <- burnin_mfc %>% 
  left_join(proj_mfc, by = "Code") %>%
  mutate(scalar = proj_mfc / burnin_mfc) %>%
  mutate(scalar = ifelse(is.nan(scalar), 1, scalar))

for (i in 1:length(codes)){
  
  sp <- codes[i]
  
  mfc_change <- scalars %>% filter(Code == sp) %>% pull(scalar)
  
  # flagchangeF 1
  old_line <- harvest_base[grep("flagchangeF", harvest_base)]
  new_line <- gsub("\t0", "\t1", old_line)
  harvest_base[grep("flagchangeF", harvest_base)] <- new_line
  
  # flagFchange_XXX - anchor with trailing space
  old_line <- harvest_base[grep(paste0("flagFchange_", sp, "\\t"), harvest_base) + 1]
  new_line <- sub("0", "1", old_line)
  harvest_base[grep(paste0("flagFchange_", sp, "\\t"), harvest_base) + 1] <- new_line
  
  # XXX_mFC_changes - sp is a prefix here so already unambiguous
  old_line <- harvest_base[grep(paste0(sp, "_mFC_changes"), harvest_base) + 1]
  new_line <- sub("0", "1", old_line)
  harvest_base[grep(paste0(sp, "_mFC_changes"), harvest_base) + 1] <- new_line
  
  # mFCchange_start_XXX - anchor with trailing space
  old_line <- harvest_base[grep(paste0("mFCchange_start_", sp, "\\t"), harvest_base) + 1]
  change_start <- 365 * 30
  new_line <- sub("0", as.character(change_start), old_line)
  harvest_base[grep(paste0("mFCchange_start_", sp, "\\t"), harvest_base) + 1] <- new_line
  
  # mFCchange_period_XXX - anchor with trailing space
  old_line <- harvest_base[grep(paste0("mFCchange_period_", sp, "\\t"), harvest_base) + 1]
  new_line <- sub("0", "1", old_line)
  harvest_base[grep(paste0("mFCchange_period_", sp, "\\t"), harvest_base) + 1] <- new_line
  
  # mFCchange_mult_XXX - anchor with trailing space
  old_line <- harvest_base[grep(paste0("mFCchange_mult_", sp, "\\t"), harvest_base) + 1]
  new_line <- as.character(round(mfc_change, digits = 5))
  harvest_base[grep(paste0("mFCchange_mult_", sp, "\\t"), harvest_base) + 1] <- new_line
}

# fix the mEff lines
harvest_base <- gsub("0.25 0.25 0.25 0.25", "0 0 0 0", harvest_base, fixed = T)

# write out
writeLines(harvest_base, "AtlantisGOA_Topic2/GOA_harvest_Frecent.prm")

# create zero-fishing variant: set all mFCchange_mult values to 0
harvest_nofishing <- harvest_base
harvest_nofishing[grep("mFCchange_mult_", harvest_nofishing) + 1] <- "0"
writeLines(harvest_nofishing, "AtlantisGOA_Topic2/GOA_harvest_F0.prm")

# Force.prm ---------------------------------------------------------------
# 4 climate scenarios
# one with no climate at all and stable conditions - this one is the base force.prm
# use that as template to create the other ones
force_file <- "AtlantisGOA_Topic2/GOA_force.prm"
force_template <- readLines(force_file)
scenarios <- c(126, 245, 585, "hind")
burnin_t <- 30
proj_t <- 80
runtime <- burnin_t + proj_t

for(i in 1:length(scenarios)){
  
  clim <- scenarios[i]
  filename <- paste0("AtlantisGOA_Topic2/GOA_force_ssp", clim,".prm")
  this_force <- force_template
  
  # hydro
  start_h <- grep("nhdfiles", this_force) - 1
  end_h <- grep("hd0.name", this_force) + 1
  
  head_file <- this_force[1:start_h]
  tail_file <- this_force[end_h:length(this_force)]
  
  h_section <- paste("nhdfiles", runtime)
  yr <- 2020
  for(j in 1:runtime){
    if(clim == "hind"){
      string <- paste0("hd", j-1, ".name forcings/hydro/goa_hydro_1999.nc")
    } else if(j <= burnin_t){
      string <- paste0("hd", j-1, ".name forcings/hydro/goa_hydro_1999.nc")
    } else {
      string <- paste0("hd", j-1, ".name ../../goa_calibration_runs/forcings_proj/ssp", clim, "/hydro/goa_hydro_", yr, ".nc")
      yr <- yr + 1
    }
    h_section <- c(h_section, string)
  }
  
  this_force <- c(head_file, h_section, tail_file)
  
  # temperature
  start_t <- grep("ntempfiles", this_force) - 1
  end_t <- grep("Temperature0.name", this_force) + 1
  
  head_file <- this_force[1:start_t]
  tail_file <- this_force[end_t:length(this_force)]
  
  t_section <- paste("ntempfiles", runtime)
  yr <- 2020
  for(j in 1:runtime){
    if(j <= burnin_t){
      string <- paste0("Temperature", j-1, ".name forcings/temp/mean_1990s_temperature.nc")
    } else if(clim == "hind"){
      string <- paste0("Temperature", j-1, ".name forcings/temp/mean_hindcast_temperature.nc")
    } else {
      string <- paste0("Temperature", j-1, ".name ../../goa_calibration_runs/forcings_proj/ssp", clim, "/temp/goa_roms_temp_", yr, ".nc")
      yr <- yr + 1
    }
    t_section <- c(t_section, string)
  }
  
  this_force <- c(head_file, t_section, tail_file)
  
  # salinity
  start_s <- grep("nsaltfiles", this_force) - 1
  end_s <- grep("Salinity0.name", this_force) + 1
  
  head_file <- this_force[1:start_s]
  tail_file <- this_force[end_s:length(this_force)]
  
  s_section <- paste("nsaltfiles", runtime)
  yr <- 2020
  for(j in 1:runtime){
    if(j <= burnin_t){
      string <- paste0("Salinity", j-1, ".name forcings/salt/mean_1990s_salinity.nc")
    } else if(clim == "hind"){
      string <- paste0("Salinity", j-1, ".name forcings/salt/mean_hindcast_salinity.nc")
    } else {
      string <- paste0("Salinity", j-1, ".name ../../goa_calibration_runs/forcings_proj/ssp", clim, "/salt/goa_roms_salt_", yr, ".nc")
      yr <- yr + 1
    }
    s_section <- c(s_section, string)
  }
  
  this_force <- c(head_file, s_section, tail_file)
  
  # plankton (SSP scenarios only)
  if(clim != "hind"){
    start_p <- grep("use_external_scaling", this_force) - 1
    
    head_file <- this_force[1:start_p]
    tail_file <- c(
      "use_external_scaling 1",
      "",
      "scale_all_mortality 0",
      "mortality_addition 0",
      "",
      paste0("externalBiologyForcingFile scalar_proj_ROMS_ssp", clim, "_110yr.nc"),
      "",
      "externalBiologyForcingFile_rewind 0"
    )
    
    this_force <- c(head_file, tail_file)
  }
  
  writeLines(this_force, filename)
}

# sh files ----------------------------------------------------------------

# now build 8 sh files that will be used to launch these runs
# Create shell scripts for all combinations of harvest and force files

# Define force files
force_files <- c("GOA_force_ssphind.prm", "GOA_force_ssp126.prm", "GOA_force_ssp245.prm", 
                 "GOA_force_ssp585.prm")

harvest_files <- c("GOA_harvest_F0", "GOA_harvest_Frecent")

# Create all combinations
combinations <- expand.grid(
  harvest_file = harvest_files,
  force_file = force_files,
  stringsAsFactors = FALSE
)

# Add run numbers (000-07)
combinations$run_id <- sprintf("%03d", 0:(nrow(combinations) - 1))

write.csv(combinations, "key.csv", row.names = F)

# Create shell scripts for each combination
for (i in 1:nrow(combinations)) {
  run_id <- combinations$run_id[i]
  harvest <- paste0(combinations$harvest_file[i], ".prm")
  force <- combinations$force_file[i]
  
  # Create shell script content matching the example format
  sh_content <- paste0(
    "atlantisMerged -i GOA_cb_summer.nc  0 -o outputGOA_", run_id, ".nc ",
    "-r GOA_run.prm -f ", force, " -p GOA_physics.prm -b GOAbioparam.prm ",
    "-h ", harvest, " -m GOAMigrations.csv -s GOA_Groups.csv ",
    "-q GOA_fisheries.csv -d outputFolder", run_id, " 2>outconsole\n"
  )
  
  # Write shell script file
  sh_filename <- paste0("AtlantisGOA_Topic2/run_atlantis_", run_id, ".sh")
  writeLines(sh_content, sh_filename)
  
  # Make it executable (on Unix-like systems)
  Sys.chmod(sh_filename, mode = "0755")
}
