#######################################################
# CWS averaged PM2.5
# For PFAS predictions in NJ
# Created by: Hari Iyer
# Date: 11/24/2025
#######################################################

library(tidyverse)
library(sf)
library(terra)
library(yaml)

#need to change this
config <- yaml.load_file("C:/Users/ms3283/Box/PFAS Grant/XGBoost/scripts/Megan/01 w new cws file/pm25_config_njcws_MS.txt")##(Kindly change with your path)

## Get CWS boundaries

nj_cws <- st_read("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/CWS_new/nj_cws_bounds_new.shp")
nj_cws_v <- vect(nj_cws)           # convert sf → SpatVector
crs(nj_cws_v)                      # check target CRS

## Make sure newPath points to folder containing rasters (temporary directory nhair0a)
list_pm25 <- list.files(path=config$pm25_dir, full.names=TRUE, pattern = ".tif$")

config$pm25_dir

bStack <- rast(list_pm25)
## reproject
rstack_reproj <- project(bStack, nj_cws_v)


## Extract and save results
njpws_p25 <- extract(rstack_reproj, nj_cws_v, fun="mean", na.rm = TRUE)
njpws_pm25_joined <- nj_cws %>%
  mutate(ID = 1:n()) %>%
  left_join(njpws_p25, by = "ID")

save(njpws_pm25_joined, file=str_c("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/Pollution/",config$exp_filename))



