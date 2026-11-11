#######################################################
# Processing and formatting predictors from Tokranov et al.
# PFAS well model using hydrology variables that have been
# pre-processed and stored as raster databases

# I will read in rasters and then export a CWS-formatted set for
## avg_clay.asc ## average clay
## depth_m.asc ## well depth
## dtw_m.asc ## depth to water
## n_sep2010.asc ## nitrogen from septic 
## rech48_mmyr.asc ## recharge

## details in Tokranov et al.

#######################################################

library(tidyverse)
library(terra)
library(sf)


nj_cws <- st_read("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/CWS_new/nj_cws_bounds_new.shp")
nj_cws_v <- vect(nj_cws)           # convert sf → SpatVector
crs(nj_cws_v)  

tokranov_albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

asc_files <- list.files("C:/Users/ms3283/Box/PFAS Grant/Tokranov paper/national_pfas_gw/national_pfas_gw/Inputs/Domestic_rasters/", pattern ="\\.asc$", full.names=TRUE)


############################
## Clay
#############################

r_stack <- rast(asc_files)
## st_crs(clay_raster) ## missing

crs(r_stack ) <- tokranov_albers
tok_reproj <- project(r_stack, nj_cws_v)


## Extract and save results
tok_cws <- extract(tok_reproj, nj_cws_v, fun="mean", na.rm = TRUE)
tok_cws_joined <- nj_cws %>%
  mutate(ID = 1:n()) %>%
  left_join(tok_cws, by = "ID")

names(tok_cws_joined)

tok_cws_joined <- tok_cws_joined %>%
  dplyr::select(PWSID, avg_clay, depth_m, dist_pfas_sites, dtw_m, 
                n_sep2010, nwalt2002_pcturb, popdensity2000, rech48_mmyr) %>%
  st_drop_geometry()

save(tok_cws_joined, file = "C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/tok_cws_joined.Rda")