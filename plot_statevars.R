# Alberto Rovellini
# 6/24/2022
# Produce plots of temp and salt

# List of packages for session
.packages = c("dplyr", "ggplot2", "purrr", "tidyr", "tidync", "sf", "ggh4x", "rbgm", "cowplot", "maps", "mapdata")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)

select <- dplyr::select

temp.file <- paste0("hydro_forcings_hindcast_revised/temp/mean_hindcast_temperature.nc")
salt.file <- paste0("hydro_forcings_hindcast_revised/salt/mean_hindcast_salinity.nc")
bgm.file <- "data/GOA_WGS84_V4_final.bgm" 
cum.depth <- c(1,30,100,200,500,1000,3969)

# make directory for plots
outdir <- 'statevar_plots_hindcast_climatology/'
dir.create(outdir)

# Atlantis model spatial domain
atlantis_bgm <- bgm.file %>% read_bgm()
atlantis_box <- atlantis_bgm %>% 
  box_sf() %>%
  st_transform(crs = 4326) %>%
  mutate(botz = -1*botz)

# Get state variables from NetCDF file ------------------------------------------
temp <- tidync(temp.file) 
salt <- tidync(salt.file)

these_vars <- hyper_grids(temp) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    temp %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

grids <- these_vars %>% filter(name=="temperature") %>% pluck('grd')

dat_temp <- temp %>% activate(grids) %>% hyper_tibble()
dat_salt <- salt %>% activate(grids) %>% hyper_tibble()

# numbering of b and z starts from 1 - change
dat_temp <- dat_temp %>% 
  mutate(b = b - 1, z = z - 1) 
dat_salt <- dat_salt %>% 
  mutate(b = b - 1, z = z - 1) 

dat <- dat_temp %>%
  left_join(dat_salt, by = c('t','b','z')) %>%
  pivot_longer(-c(t,b,z))

# here 0 is the bottom

# plot summary across model domain 
# dat %>%
#   left_join(atlantis_box %>% st_set_geometry(NULL) %>% select(box_id, boundary), by = c('b'='box_id')) %>%
#   filter(name == 'temperature', z == 0, boundary == FALSE) %>%
#   group_by(t) %>%
#   summarize(mean_val = mean(value)) %>%
#   mutate(year = 2017) %>%
#   write.csv('temp2017.csv', row.names = F)

# Loop for all boxes ------------------------------------------------------
all_boxes <- atlantis_box %>% pull(box_id) 

make_statevar_plot <- function(focal_box){
  
  p1 <- dat %>%
    filter(b == focal_box,
           z != 6) %>% # dropping the sediment
    group_by(t,b) %>%
    mutate(z_flip = max(z)-z) %>% # here is where 0 goes from being bottom to be surface
    ungroup() %>%
    mutate(t = as.POSIXct(t, origin = '2000-01-01', tz = 'UTC')) %>%
    ggplot(aes(x = t, y = value, color = factor(z_flip)))+
    geom_line(size = 1.5)+
    scale_color_brewer(palette = 'Dark2')+
    theme_bw()+
    labs(x = '', title = paste('State variables in box', focal_box, '(0 is the surface)'), color = 'Layer')+
    facet_wrap(~name, scales = 'free', 
               labeller = as_labeller( c(temperature = 'Temperature (C)', salinity = 'Salinity (PPT)') ))+
    theme(strip.placement = "outside")
  
  this_bbox <- st_bbox(atlantis_box %>% filter(box_id == focal_box))
  # add a small buffer around it?
  this_bbox[1] <- this_bbox[1] - 0.1 #xmin
  this_bbox[2] <- this_bbox[2] - 0.1 #ymin
  this_bbox[3] <- this_bbox[3] + 0.1 #xmax
  this_bbox[4] <- this_bbox[4] + 0.1 #ymax
  
  this_atlantis_box <- atlantis_box %>% 
    st_crop(this_bbox)%>% 
    mutate(lyrs = findInterval(botz, cum.depth)) 
  
  cc <- scales::seq_gradient_pal("yellow", "blue", "Lab")(seq(0,1,length.out = (max(this_atlantis_box$lyrs)+1)))
  
  coast <- maps::map("worldHires", regions = c("Canada", "USA"), plot = FALSE, fill = TRUE)
  coast_sf <- coast %>% st_as_sf() %>% st_transform(crs = 4326) %>% st_combine()
  
  p2 <- this_atlantis_box %>%
    ggplot()+
    geom_sf(aes(fill = factor(lyrs)), color = 'white')+
    scale_fill_manual(values = cc) +
    geom_sf(data=coast_sf)+
    coord_sf(xlim = c(this_bbox$xmin,
                      this_bbox$xmax),
             ylim = c(this_bbox$ymin,
                      this_bbox$ymax))+
    geom_sf_label(aes(label = box_id))+
    labs(x = '', y = '', fill = 'Depth layers')+
    theme_bw()
  
  # Produce a graphic output ------------------------------------------------
  
  p3 <- plot_grid(p2, p1, ncol = 1, nrow = 2, rel_heights = c(0.8, 1))
  save_plot(paste0(outdir, focal_box, '.png'), p3, base_height = 12, base_width = 18)
  
}

purrr::map(all_boxes, possibly(make_statevar_plot,NA))
  