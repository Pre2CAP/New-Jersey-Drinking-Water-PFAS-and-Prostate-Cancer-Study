#############################################################################
# R Code to Process Distance from Centroid to Airports
# Created by: Hari Iyer
# Date: 12/1/2025

# I then calculated the distance from the census block group centroid to that facility, saving
# up to 3 locations.

# Here, I average across those distances and save to one row per CWS.

############################################################################

library(tidyverse)
library(sf)

nj_bg_airp <- st_read("C:/Users/ms3283/Box/PFAS Grant/NJDOT_Public_Airports/nj_bg_airport.shp") ## 18960

## average by block group first

nj_pfas_prep <- nj_bg_airp %>%
  dplyr::select(GEOID, NEAR_DIST, NEAR_RANK) %>%
  st_drop_geometry() %>%
  group_by(GEOID) %>%
  summarise(avg_dist_airp = mean(NEAR_DIST))

## merge back to census block shapefile

nj_bg <- st_read("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/NJ Blocks/nj_bg.shp")

nj_bg_pfas_prep <- nj_bg %>%
  left_join(nj_pfas_prep, by = "GEOID")


## aggregate to CWS


nj_cws <- st_read("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/CWS_new/nj_cws_bounds_new.shp")



st_crs(nj_cws)
st_crs(nj_bg_pfas_prep)

nj_bg_p <- st_transform(nj_bg_pfas_prep, crs = st_crs(nj_cws))

intersections <- st_intersection(
  nj_cws %>% dplyr::select(PWSID),
  nj_bg_p %>% dplyr::select(avg_dist_airp)
)

nj_cws_pfas_dist <- intersections %>%
  mutate(area = st_area(.)) %>%
  group_by(PWSID) %>%
  summarize(
    awairport = sum(as.numeric(avg_dist_airp * area)) / sum(as.numeric(area)),
    .groups = "drop"
  )


st_write(nj_cws_pfas_dist, "C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/nj_cws_airp_dist.shp", delete_dsn = TRUE)


