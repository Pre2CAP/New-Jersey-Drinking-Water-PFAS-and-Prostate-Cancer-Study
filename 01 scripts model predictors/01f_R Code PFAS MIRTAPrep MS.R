#############################################################################
# R Code to Process Distance from Centroid to PFAS industrial contamination
# Created by: Hari Iyer
# Date: 11/30/2025

# Steps - in ArcGIS Pro, I filtered Business location databases from 2006 to 2024
# for key industrial facility codes outlined in prior research (notes elsewhere).
# There was one row per business, and the business had to be open during the PFAS sampling
# data period.

# I then calculated the distance from the census block group centroid to that facility, saving
# up to 3 locations.

# Here, I average across those distances and save to one row per CWS. 

############################################################################

library(tidyverse)
library(sf)

nj_bg_mirta <- st_read("C:/Users/ms3283/Box/PFAS Grant/MIRTA/nj_bg_centroids_mirta.shp") ## 18960

## average by block group first

nj_pfas_prep <- nj_bg_mirta %>%
  dplyr::select(GEOID, NEAR_DIST, NEAR_RANK) %>%
  st_drop_geometry() %>%
  group_by(GEOID) %>%
  summarise(avg_dist_mirta = mean(NEAR_DIST))

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
  nj_bg_p %>% dplyr::select(avg_dist_mirta)
)

nj_cws_pfas_dist <- intersections %>%
  mutate(area = st_area(.)) %>%
  group_by(PWSID) %>%
  summarize(
    awmirta = sum(as.numeric(avg_dist_mirta * area)) / sum(as.numeric(area)),
    .groups = "drop"
  )


st_write(nj_cws_pfas_dist, "C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/nj_cws_mirta_dist.shp", delete_dsn = TRUE)


