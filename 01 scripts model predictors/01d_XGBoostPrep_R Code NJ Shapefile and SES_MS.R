###################################################################
# Get NSES for NJ at Tract Level, Average to CWS
# Created by: Hari Iyer
# Date: 11/24/2025
###################################################################

## install.packages("tidycensus")
## install.packages("tigris")

library(tidyverse)
library(tidycensus)
library(tigris)
library(sf)
library(terra)

census_api_key("057ce44279d4e8fb68bc00ba20284dd386be7fb8")

options(tigris_use_cache = TRUE)

## 2010 census tracts

njtract <- tracts("NJ", year=2010) %>%
  rename(GEOID =  GEOID10) %>%
  ## st_drop_geometry() %>%
  mutate(area = ALAND10/2589988) %>% ## converting to sq miles
  select(GEOID, area)

njtract$GEOID <- as.numeric(njtract$GEOID)

nj_ses10 <- read.csv("//cinjasbp138.rad.rutgers.edu/pre2cap//Registries/data/States NSES/NJ_NSES_2010.csv")

nj_ses10 <- nj_ses10 %>%
  dplyr::select(-area)

njtract_nses_10 <- njtract %>%
  left_join(nj_ses10, by = "GEOID")

## read in CWS

nj_cws <- st_read("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/CWS_new/nj_cws_bounds_new.shp")

## Take zonal stat average for SES measures within CWS

st_crs(nj_cws)
st_crs(njtract_nses_10)

njtract_nses_10_p <- st_transform(njtract_nses_10, crs = st_crs(nj_cws))

vars <- c("edult810", "edugt8lths10", "eduhs10",  "edubach10","mdhval10", "mdinc10",       
  "percwhite10","percblack10","percforeign10","percpov10","pohse10","pfwchbf10",     
  "unemp10",  "munemp10",   "percfamint10", "percwelfr10",  "percsales10","percclerical10",
  "percrepair10", "perctransp10", "percnflabor10","percservice10","percchild10","percold10",     
  "mdrent10", "incice10", "rice10", "rcincice10")


cws_small   <- nj_cws %>% select(PWSID)
tract_small <- njtract_nses_10_p %>% select(all_of(vars))


# Intersect CWS and tracts
ints <- st_intersection(cws_small, tract_small)

# Compute area-weighted means for all 'vars'
cws_aw_means <- ints %>%
  mutate(area = as.numeric(st_area(.))) %>%     # area of each intersected piece
  st_drop_geometry() %>%                        # we just need attributes now
  group_by(PWSID) %>%
  summarise(
    across(
      all_of(vars),
      ~ sum(as.numeric(.x) * area, na.rm = TRUE) / sum(area[!is.na(.x)]),
      .names = "aw_{.col}"                      # prefix output cols with 'aw_'
    ),
    .groups = "drop"
  )

## save as data frame
save(cws_aw_means, file = "C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/cws_aw_means.Rda")

# stop here

st_write(njtract_nses_10, "Z:/Registries/data/New Jersey/data/NJ tract shapefile/2010/njtract_nses_10.shp")


############################################################
# Census Block Group
############################################################

nj_block <- get_acs(
  geography = "block group",
  variables = c("B19013_001","B01003_001"),
  state = "NJ",
  year = 2013,
  geometry = TRUE
) 

nj_block <- nj_block %>%
  filter(variable == "B01003_001") %>%
  st_drop_geometry()
  
nj_bg_tiger <- block_groups(state = "NJ", year = 2013, class = "sf")

v16 <- load_variables(2013, "acs5", cache = TRUE)

nj_bg <- nj_bg_tiger %>%
  left_join(nj_block, by = "GEOID") %>%
  mutate(
    area_sqmi = ALAND/2589988,
    pop_density_sqmi = estimate/area_sqmi
  )

nj_bg <- st_transform(nj_bg, 5070)

nj_bg <- nj_bg %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(cx = st_coordinates(centroid)[,1],
         cy = st_coordinates(centroid)[,2])


st_write(nj_bg, "C:/Users/hi97/Box/Registry transfer/NJ/NJ Blocks/nj_bg.shp")

