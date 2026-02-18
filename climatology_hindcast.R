# Alberto Rovellini
# 07/03/2025
# This code takes salt and salt nc files for 1991-1999 and creates a climatology
# then it writes the files out to one forcing file each 

library(tidyverse)
library(ncdf4)
library(tidync)

# apply the function we have created for the delta correction as the averaging is the same process
# merge corrected projection files
shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" mergetime hydro_forcings_hindcast_revised\\hydro\\*.nc hydro_forcings_hindcast_revised\\hydro\\hind_clim_merged.nc")

# what variable are we handling?
this_variable <- "temperature"

# files
hindcast_file <- "hydro_forcings_hindcast_revised/temp/hind_clim_merged.nc"

# apply function
mean_hind <- make_delta_array(variable = this_variable,
                              sim_period = "hindcast",
                              mean = T,
                              leap = F)

# pack to netcdf
# for the purpose of not having to rebuild the netcdf from scratch, put the new variable to an exisiting file
mean_hind_file <- paste0("hydro_forcings_hindcast_revised/temp/mean_hindcast_", this_variable, ".nc")

# copy the original necdf file that we are correcting
file.copy("hydro_forcings_hindcast_revised/temp/goa_roms_temp_2000.nc", mean_hind_file)

# open the new file
mean_hind_nc <- nc_open(mean_hind_file, write = T)

# write the delta-corrected variable to it
ncvar_put(mean_hind_nc, varid = this_variable, vals = mean_hind)

# close the file
nc_close(mean_hind_nc)
