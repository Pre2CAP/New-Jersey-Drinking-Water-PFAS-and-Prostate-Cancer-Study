###########################################################################
# Formatting and preparing for PFAS modeling
# Created by: Megan Shanahan, Hari Iyer
# Created inverse distance to MIRTA, firefighters, FUDS, airport
# Date: 4/24/26

# Updated to have full coverage of time periods
###########################################################################

library(tidyverse)
library(sf)

nj_all_centroids <- st_read("C:/Users/ms3283/Box/PFAS Grant/ArGIS Pro Projects/IW_combined.shp")
names(nj_all_centroids)

length(unique(nj_all_centroids$GEOID))

nj_prox_part1 <- nj_all_centroids %>%
  rename(IW_FIRE_5K_5 =	IW_FIRE_5K,
                IW_FIRE_5K_10 =	IW_FIRE__1,
                IW_FIRE_5K_15 =	IW_FIRE__2,
                IW_FIRE_3K_5 =	IW_FIRE_3K,
                IW_FIRE_3K_10 =	IW_FIRE__3,
                IW_FIRE_1K_5 =	IW_FIRE_1K,
                IW_AP_10K_5 =	IW_AP_10K_,
                IW_AP_5K_5 =	IW_AP_5K_5,
                IW_AP_3K_5 =	IW_AP_3K_5,
                IW_AP_1K_5 =	IW_AP_1K_5,
                IW_FUDS_10K_10 =	IW_FUDS_10,
                IW_FUDS_5K_5 =	IW_FUDS_5K,
                IW_FUDS_3K_5 =	IW_FUDS_3K,
                IW_FUDS_1K_5 =	IW_FUDS_1K,
                IW_MIRTA_10K_5 =	IW_MIRTA_1,
                IW_MIRTA_5K_5 =	IW_MIRTA_5,
                IW_MIRTA_3K_5 =	IW_MIRTA_3) %>%
  mutate(GEOID = as.numeric(GEOID))

length(unique(nj_prox_part1$GEOID)) 
  
nj_bus_locations <- read.csv("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/Business Locations - UPDATED/nj_bg_bus_all.csv")
names(nj_bus_locations)

nj_bus_locations <- nj_bus_locations %>%
  mutate(across(c(IDW_03_05, IDW_03_10, IDW_05_05, IDW_05_10, IDW_05_15, IDW_01_05), ~replace_na(., 0))) %>%
  filter(!is.na(GEOID))## assigning 0 to 
## all rows with no data because no proximity-based locations were found


length(unique(nj_bus_locations$GEOID)) 

## CHECKING DISTRIBUTION OF ADDRESSES

year_counts <- nj_bus_locations %>%
  group_by(GEOID) %>%
  summarise(n_years = n_distinct(year))

# Check if all GEOIDs have the same number of years
table(year_counts$n_years)

year_details <- nj_bus_locations %>%
  group_by(GEOID) %>%
  summarise(
    n_years = n_distinct(year),
    years = list(sort(unique(year)))
  )

year_details %>%
  filter(n_years < max(n_years))

year_details %>%
  mutate(years_str = sapply(years, paste, collapse = ", ")) %>%
  select(GEOID, n_years, years_str)

target_years <- c(2000, 2005, 2010, 2015, 2020)

nj_bus_expanded <- nj_bus_locations %>%
  mutate(year = as.numeric(year)) %>%
  group_by(GEOID) %>%
  complete(year = target_years) %>%
  arrange(GEOID, year) %>%
  fill(-GEOID, -year, .direction = "downup") %>%
  ungroup()

sum(is.na(nj_bus_locations))
sum(is.na(nj_bus_expanded))


summary(nj_bus_expanded)

length(unique(nj_bus_locations$GEOID))
length(unique(nj_prox_part1$GEOID))

nj_idw_merged <- nj_bus_expanded %>%
  full_join(nj_prox_part1, by = "GEOID") %>%
  dplyr::select(GEOID, year, IDW_03_05, IDW_03_10, IDW_05_05, IDW_05_10, IDW_05_15, IDW_01_05,
                pp_dns_, IW_FIRE_5K_5, IW_FIRE_5K_10, IW_FIRE_5K_15, IW_FIRE_3K_5,
                IW_FIRE_3K_10, IW_FIRE_1K_5, IW_AP_10K_5, IW_AP_5K_5,
                IW_AP_3K_5, IW_AP_1K_5, IW_FUDS_10K_10, IW_FUDS_5K_5, IW_FUDS_3K_5, IW_FUDS_1K_5,
                IW_MIRTA_10K_5, IW_MIRTA_5K_5, IW_MIRTA_3K_5) %>%
  st_drop_geometry() %>%
  filter(!is.na(GEOID)) 

length(unique(nj_idw_merged$GEOID)) 
summary(nj_idw_merged$year)
summary(nj_idw_merged$IDW_03_05)


year_counts <- nj_idw_merged %>%
  group_by(GEOID) %>%
  summarise(n_years = n_distinct(year))

# Check if all GEOIDs have the same number of years
table(year_counts$n_years)


nj_idw_merged <- nj_idw_merged %>%
  group_by(GEOID) %>%
  complete(year = target_years) %>%
  arrange(GEOID, year) %>%
  fill(-GEOID, -year, .direction = "downup") %>%
  ungroup() %>%
  filter(!is.na(year)) %>%
  mutate(across(c(IDW_03_05, IDW_03_10, IDW_05_05, IDW_05_10, IDW_05_15, IDW_01_05), ~replace_na(., 0)))


summary(nj_idw_merged)
table(nj_idw_merged$year)

save(nj_idw_merged, file="C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/nj_idw_merged.Rda")

###################################################################
# Merge with shape file for BG, then aggregate to PWS
###################################################################

nj_pfas_indust <- st_read("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/Business to Centroids/nj_dist_pfas.shp") ## 18960

## average by block group first

nj_pfas_prep <- nj_pfas_indust %>%
  dplyr::select(GEOID, NEAR_DIST, NEAR_RANK) %>%
  st_drop_geometry() %>%
  group_by(GEOID) %>%
  summarise(avg_dist_pfas = mean(NEAR_DIST))

## merge back to census block shapefile

nj_bg <- st_read("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/NJ Blocks/nj_bg.shp")

## aggregate to CWS

nj_cws <- st_read("C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/CWS_new/nj_cws_bounds_new.shp")

st_crs(nj_cws)
st_crs(nj_bg)

nj_bg_p <- st_transform(nj_bg, crs = st_crs(nj_cws))



crosswalk <- st_intersection(nj_bg_p, nj_cws)

crosswalk <- crosswalk %>%
  mutate(
    area_int = st_area(.),
    area_cbg = st_area(nj_bg_p[match(GEOID, nj_bg_p$GEOID), ]),
    w = as.numeric(area_int / area_cbg),
    GEOID = as.numeric(GEOID)
  ) %>%
  dplyr::select(GEOID, PWSID, area_int, area_cbg, w) %>%
  st_drop_geometry()

## THIS IS NORMALIZING THE WEIGHTS SO THAT THEY ADD TO 1 for each GEOID

crosswalk <- crosswalk %>%
  group_by(GEOID) %>%
  mutate(w = w / sum(w, na.rm = TRUE)) %>%
  ungroup()

crosswalk %>%
  group_by(GEOID) %>%
  summarise(w_sum = sum(w, na.rm = TRUE)) %>%
  summary()


length(unique(crosswalk$GEOID)) 
length(unique(crosswalk$PWSID))

length(unique(nj_idw_merged$GEOID)) 
names(nj_idw_merged)

nj_cws_idw <- nj_idw_merged %>%
  left_join(crosswalk, by = "GEOID") %>%
  group_by(PWSID, year) %>%
  summarise(
    across(
      .cols = c(IDW_03_05, IDW_03_10, IDW_05_05, IDW_05_10, IDW_05_15,
                pp_dns_, IW_FIRE_5K_5, IW_FIRE_5K_10, IW_FIRE_5K_15, IW_FIRE_3K_5,
                IW_FIRE_3K_10, IW_FIRE_1K_5, IW_AP_10K_5, IW_AP_5K_5,
                IW_AP_3K_5, IW_AP_1K_5, IW_FUDS_10K_10, IW_FUDS_5K_5, IW_FUDS_3K_5, IW_FUDS_1K_5,
                IW_MIRTA_10K_5, IW_MIRTA_5K_5, IW_MIRTA_3K_5),
      .fns = ~ sum(.x * w, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(!is.na(PWSID))

summary(nj_cws_idw)
length(unique(nj_cws_idw$PWSID))


save(nj_cws_idw, file="C:/Users/ms3283/Box/PFAS Grant/XGBoost/data/nj_cws_idw.Rda")


