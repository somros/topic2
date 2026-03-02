# extract outputs for Topic 2 runs
# make some plots to check that everything ran the way it should have
# expect similar results if paper runs are an indication

library(tidyverse)

setwd("~/topic2/")

grps <- read.csv("AtlantisGOA_Topic2/GOA_Groups.csv", header = T)
codes <- grps %>% pull(Code)
codes_t2 <- c("POL","COD","ATF")

# maturity at age
bio_prm <- paste0("AtlantisGOA_Topic2/GOAbioparam.prm")
bio <- readLines(bio_prm)

fspb_df <- data.frame()
for(i in 1:length(codes_t2)){
  
  sp <- codes_t2[i]
  fspb_line <- bio[grep(paste0("FSPB_", sp), bio) + 2]
  fspb <- as.numeric(strsplit(fspb_line, split = " ")[[1]])
  
  fspb_sp <- data.frame("Code" = rep(sp, length(fspb)), 
                        "age" = 0:(length(fspb)-1), 
                        "fspb" = fspb)
  
  fspb_df <- rbind(fspb_df, fspb_sp)
  
}


pull_fishery_info_t2 <- function(this_run) {
  
  print(this_run)
  
  # Build paths
  run_str  <- formatC(this_run, width = 3, flag = "0")
  wd       <- paste0("~/topic2/AtlantisGOA_Topic2_run", run_str,
                     "/outputFolder", run_str)
  biom_file  <- paste0("outputGOA_", run_str, "AgeBiomIndx.txt")
  catch_file <- paste0("outputGOA_", run_str, "Catch.txt")
  
  # ── Read ──────────────────────────────────────────────────────────────────
  biom  <- read.csv(file.path(wd, biom_file),  sep = " ", header = TRUE)
  catch <- read.csv(file.path(wd, catch_file), sep = " ", header = TRUE)
  
  # Force numeric (except Time)
  biom  <- biom  %>% mutate(across(-Time, ~as.numeric(as.character(.))))
  catch <- catch %>% mutate(across(-Time, ~as.numeric(as.character(.))))
  
  # Drop trailing incomplete rows
  # trim_complete <- function(df) {
  #   complete_rows <- apply(df[, -1], 1,
  #                          function(x) all(!is.na(suppressWarnings(as.numeric(as.character(x))))))
  #   df[1:max(which(complete_rows)), ]
  # }
  # biom  <- trim_complete(biom)
  # catch <- trim_complete(catch)
  
  # Optional time ceiling
  if (exists("yr_end")) {
    biom  <- biom  %>% filter(Time / 365 <= yr_end)
    catch <- catch %>% filter(Time / 365 <= yr_end)
  }
  
  # ── Biomass: pivot to long, split Code / Age ───────────────────────────
  biom_long <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") %>%
    separate(Code.Age, into = c("Code", "Age"), sep = "\\.", convert = TRUE)
  
  biom_filtered <- biom_long %>%
    filter(Code %in% codes_t2, Time < max(Time)) %>%
    filter(Time > 0) # drop initial condition biomass as not part of year 1 in principle
  
  # ── Helper: collapse across ages, then take max value per year ────────
  annual_max <- function(df, value_col) {
    df %>%
      group_by(Time, Code) %>%
      summarise(val = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
      mutate(Year = ceiling(Time / 365)) %>%
      group_by(Year, Code) %>%
      slice_max(val, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(-Time) %>%
      rename(!!value_col := val)
  }
  
  # ── Total biomass (sum all ages, max per year) ─────────────────────────
  biom_tot <- biom_filtered %>%
    annual_max("mt") %>%
    rename(biom_mt_tot = mt)
  
  # ── SSB (fspb-weighted biomass, max per year) ──────────────────────────
  biom_ssb <- biom_filtered %>%
    left_join(fspb_df, by = c("Code", "Age" = "age")) %>%
    mutate(ssb_mt = mt * fspb) %>%
    annual_max("ssb_mt")
  
  # ── Catch: keep only relevant species, annual total ────────────────────
  # Year = 0 is dropped, those catches equal 0
  catch_long <- catch %>%
    select(Time, any_of(codes_t2)) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt") %>%
    mutate(Year = ceiling(Time / 365)) %>%
    select(-Time) %>%
    filter(Year > 0) # drop zero catch
  
  # ── Combine ────────────────────────────────────────────────────────────
  result <- biom_tot %>%
    left_join(biom_ssb,   by = c("Year", "Code")) %>%
    left_join(catch_long, by = c("Year", "Code")) %>%
    mutate(
      run_id  = this_run
    )
  
  return(result)
}

# ── Run across 002–007 ────────────────────────────────────────────────────
run_ids <- 2:7

all_results <- map_dfr(run_ids, pull_fishery_info_t2)

# bring in run info
key <- read.csv("data/key.csv", header = T)
key <- key %>%
  mutate(
    f   = sub(".*_", "", harvest_file),
    env = sub(".*ssp(.*)\\.prm", "ssp\\1", force_file)
  )

all_results <- all_results %>% left_join(key, by = "run_id")

# view
all_results %>%
  select(Year, Code, biom_mt_tot, ssb_mt, catch_mt, f, env) %>%
  pivot_longer(-c(Year, Code, f, env), names_to = "var", values_to = "mt") %>%
  ggplot(aes(x = Year, y = mt, colour = env))+
  geom_line(aes(linetype = f))+
  scale_y_continuous(limits = c(0,NA))+
  ggh4x::facet_grid2(Code ~ var, scales = "free", independent = "y")

