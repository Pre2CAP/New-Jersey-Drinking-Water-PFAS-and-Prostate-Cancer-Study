#############################################################################
# Generate time-varying predictions of <ND, ND-<MCL, >=MCL
# For PFOA, PFOS, PFNA
# For all 566 CWS
# Date: 10/23/2026

# Created by: Hari
# Goal: Map with cutpoints defined by different PFAS to improve
# ability to see geographic patterns
#############################################################################


library(tidyverse)
library(sf)
library(terra)
library(xgboost)
library(DiagrammeR)
library(pROC)

## install.packages("xgboost")

## install.packages("tmap")
library(tmap)
library(caret)



## PFOA
## install.packages("viridis")
library(viridis)
library(grid)

## from C:\Users\hi97\Box\PFAS Grant\XGBoost\models\scored
load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/pfoa_model_score.Rda")
load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/pfos_model_score.Rda")
load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/pfna_model_score.Rda")

pfoa_0612 <- pfoa_model_score %>%
  filter(year >= 2006 & year <= 2012) %>%
  group_by(PWSID) %>%
  summarise(avg_pfoa = mean(pfoa_exp))


pfoa_1318 <- pfoa_model_score %>%
  filter(year >= 2013 & year <= 2018) %>%
  group_by(PWSID) %>%
  summarise(avg_pfoa = mean(pfoa_exp))

pfoa_1925 <- pfoa_model_score %>%
  filter(year >= 2019 & year <= 2025) %>%
  group_by(PWSID) %>%
  summarise(avg_pfoa = mean(pfoa_exp))



## READ IN CWS

## 2006-2012


nj_cws <- st_read("/extra/home/hi97/R/registry/nj_res_hx/nj_cws_new/nj_cws_bounds_new.shp")

nj_cws_pfoa_0612 <- nj_cws %>%
  left_join(pfoa_0612, by = "PWSID")

nj_cws_pfoa_1318 <- nj_cws %>%
  left_join(pfoa_1318, by = "PWSID")

nj_cws_pfoa_1925 <- nj_cws %>%
  left_join(pfoa_1925, by = "PWSID")

##

all_vals <- c(nj_cws_pfoa_0612$avg_pfoa,
              nj_cws_pfoa_1318$avg_pfoa,
              nj_cws_pfoa_1925$avg_pfoa)

# Jenks breaks (requires classInt)
## install.packages("classInt")
library(classInt)
#brks <- classIntervals(all_vals, n = 5, style = "jenks")$brks
pfoa_brks <- c(-Inf, 1.11, 1.50, 2.65, 5.29, Inf) ## match quintiles from paper
pfoa_cols <- viridis(length(pfoa_brks) - 1)

m1 <- tm_shape(nj_cws_pfoa_0612) +
  tm_polygons(
    col = "avg_pfoa",
    palette = pfoa_cols,
    style="fixed", breaks=pfoa_brks,
    title = "Predicted PFOA (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  tm_layout(
    main.title = "PFOA (ng/L) by CWS (2006-2012)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    ## legend.width = 0.8,
    ## legend.height = 0.5,
    frame = FALSE
  )

m2 <- tm_shape(nj_cws_pfoa_1318) +
  tm_polygons(
    col = "avg_pfoa",
    palette = pfoa_cols,
    style="fixed", breaks=pfoa_brks,
    title = "Predicted PFOA (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  tm_layout(
    main.title = "PFOA (ng/L) by CWS (2013-2018)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    ## legend.width = 0.8,
    ## legend.height = 0.5,
    frame = FALSE
  )


m3 <- tm_shape(nj_cws_pfoa_1925) +
  tm_polygons(
    col = "avg_pfoa",
    palette = pfoa_cols,
    style="fixed", breaks=pfoa_brks,
    title = "Predicted PFOA (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  tm_layout(
    main.title = "PFOA (ng/L) by CWS (2019-2025)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    ## legend.width = 0.8,
    ## legend.height = 0.5,
    frame = FALSE
  )


tmap_arrange(
  m1,m2,m3, 
  nrow=1)

## END
summary(pfoa_0612$avg_pfoa)
summary(pfoa_1318$avg_pfoa)
summary(pfoa_1925$avg_pfoa)


############################################################# PFOS ##################################################################


pfos_0612 <- pfos_model_score %>%
  filter(year >= 2006 & year <= 2012) %>%
  group_by(PWSID) %>%
  summarise(avg_pfos = mean(pfos_exp))


pfos_1318 <- pfos_model_score %>%
  filter(year >= 2013 & year <= 2018) %>%
  group_by(PWSID) %>%
  summarise(avg_pfos = mean(pfos_exp))

pfos_1925 <- pfos_model_score %>%
  filter(year >= 2019 & year <= 2025) %>%
  group_by(PWSID) %>%
  summarise(avg_pfos = mean(pfos_exp))



## READ IN CWS

## 2006-2012

nj_cws_pfos_0612 <- nj_cws %>%
  left_join(pfos_0612, by = "PWSID")

nj_cws_pfos_1318 <- nj_cws %>%
  left_join(pfos_1318, by = "PWSID")

nj_cws_pfos_1925 <- nj_cws %>%
  left_join(pfos_1925, by = "PWSID")

##

all_vals <- c(nj_cws_pfos_0612$avg_pfos,
              nj_cws_pfos_1318$avg_pfos,
              nj_cws_pfos_1925$avg_pfos)

# Jenks breaks (requires classInt)
## install.packages("classInt")
library(classInt)
#brks <- classIntervals(all_vals, n = 5, style = "jenks")$brks
pfos_brks <- c(-Inf, 1.03, 1.14, 1.53, 3.18, Inf)
pfos_cols <- pfoa_cols[1:(length(pfos_brks) - 1)]

m4 <- tm_shape(nj_cws_pfos_0612) +
  tm_polygons(
    col = "avg_pfos",
    palette = pfos_cols,
    style="fixed", breaks=pfos_brks,
    title = "Predicted PFOS (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  tm_layout(
    main.title = "PFOS (ng/L) by CWS (2006-2012)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    ## legend.width = 0.8,
    ## legend.height = 0.5,
    frame = FALSE
  )

m5 <- tm_shape(nj_cws_pfos_1318) +
  tm_polygons(
    col = "avg_pfos",
    palette = pfos_cols,
    style="fixed", breaks=pfos_brks,
    title = "Predicted PFOS (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  tm_layout(
    main.title = "PFOS (ng/L) by CWS (2013-2018)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    ## legend.width = 0.8,
    ## legend.height = 0.5,
    frame = FALSE
  )


m6 <- tm_shape(nj_cws_pfos_1925) +
  tm_polygons(
    col = "avg_pfos",
    palette = pfos_cols,
    style="fixed", breaks=pfos_brks,
    title = "Predicted PFOS (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  tm_layout(
    main.title = "PFOS (ng/L) by CWS (2019-2025)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    ## legend.width = 0.8,
    ## legend.height = 0.5,
    frame = FALSE
  )


tmap_arrange(m4,m5,m6, nrow=1)

summary(pfos_0612$avg_pfos)
summary(pfos_1318$avg_pfos)
summary(pfos_1925$avg_pfos)

################################################ PFNA ##################################################################################

pfna_0612 <- pfna_model_score %>%
  filter(year >= 2006 & year <= 2012) %>%
  group_by(PWSID) %>%
  summarise(avg_pfna = mean(pfna_exp))


pfna_1318 <- pfna_model_score %>%
  filter(year >= 2013 & year <= 2018) %>%
  group_by(PWSID) %>%
  summarise(avg_pfna = mean(pfna_exp))

pfna_1925 <- pfna_model_score %>%
  filter(year >= 2019 & year <= 2025) %>%
  group_by(PWSID) %>%
  summarise(avg_pfna = mean(pfna_exp))



## READ IN CWS

## 2006-2012

nj_cws_pfna_0612 <- nj_cws %>%
  left_join(pfna_0612, by = "PWSID")

nj_cws_pfna_1318 <- nj_cws %>%
  left_join(pfna_1318, by = "PWSID")

nj_cws_pfna_1925 <- nj_cws %>%
  left_join(pfna_1925, by = "PWSID")

##

all_vals <- c(nj_cws_pfna_0612$avg_pfna,
              nj_cws_pfna_1318$avg_pfna,
              nj_cws_pfna_1925$avg_pfna)




# Jenks breaks (requires classInt)
## install.packages("classInt")
library(classInt)
#brks <- classIntervals(all_vals, n = 5, style = "jenks")$brks
pfna_brks <- c(-Inf, 2, 3.8, Inf)
pfna_cols <- viridis(3)


# zoomed in code

high_pfna_0612 <- nj_cws_pfna_0612 %>%
  filter(!is.na(avg_pfna), avg_pfna > 5)

zoom_bb_raw <- st_bbox(high_pfna_0612)

pad <- 10000

zoom_bb <- st_bbox(
  c(
    xmin = as.numeric(zoom_bb_raw["xmin"]) - pad,
    xmax = as.numeric(zoom_bb_raw["xmax"]) + pad,
    ymin = as.numeric(zoom_bb_raw["ymin"]) - pad,
    ymax = as.numeric(zoom_bb_raw["ymax"]) + pad
  ),
  crs = st_crs(nj_cws_pfna_0612)
)

# Crop data for inset
nj_cws_pfna_0612_zoom <- st_crop(nj_cws_pfna_0612, zoom_bb)

# --------------------------------------------------
# 2. Main full map
# --------------------------------------------------

main_map <- tm_shape(nj_cws_pfna_0612) +
  tm_polygons(
    col = "avg_pfna",
    palette = pfna_cols,
    style = "fixed",
    breaks = pfna_brks,
    title = "Predicted PFNA (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  
  # red box showing zoom area
  tm_shape(st_as_sfc(zoom_bb)) +
  tm_borders(col = "red", lwd = 2) +
  
  tm_layout(
    main.title = "PFNA (ng/L) by CWS (2006-2012)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    frame = FALSE
  )

# --------------------------------------------------
# 3. Zoomed inset map
# --------------------------------------------------

zoom_map <- tm_shape(nj_cws_pfna_0612_zoom) +
  tm_polygons(
    col = "avg_pfna",
    palette = pfna_cols,
    style = "fixed",
    breaks = pfna_brks,
    title = ""
  ) +
  tm_borders(col = "black", lwd = 1) +
  tm_layout(
    legend.show = FALSE,
    frame = TRUE,
    inner.margins = 0
  )

# --------------------------------------------------
# 4. Print full map and overlay inset
# --------------------------------------------------

print(main_map)

print(
  zoom_map,
  vp = viewport(
    x = 0.22,
    y = 0.25,
    width = 0.32,
    height = 0.32
  )
)

m8 <- tm_shape(nj_cws_pfna_1318) +
  tm_polygons(
    col = "avg_pfna",
    palette = pfna_cols,
    style="fixed", breaks=pfna_brks,
    title = "Predicted PFNA (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  tm_layout(
    main.title = "PFNA (ng/L) by CWS (2013-2018)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    ## legend.width = 0.8,
    ## legend.height = 0.5,
    frame = FALSE
  )


m9 <- tm_shape(nj_cws_pfna_1925) +
  tm_polygons(
    col = "avg_pfna",
    palette = pfna_cols,
    style="fixed", breaks=pfna_brks,
    title = "Predicted PFNA (ng/L)"
  ) +
  tm_borders(col = "grey40", lwd = 0.5) +
  tm_layout(
    main.title = "PFNA (ng/L) by CWS (2019-2025)",
    main.title.position = "center",
    main.title.size = 1.2,
    legend.outside = TRUE,
    ## legend.width = 0.8,
    ## legend.height = 0.5,
    frame = FALSE
  )


grid.newpage()

pushViewport(viewport(
  layout = grid.layout(
    nrow=1, 
    ncol=4,
    widths = unit(c(1.2,3,3,3), "null")
                )
  ))

# Inset zoomed map in its own left-hand column
pushViewport(viewport(
  layout.pos.row = 1,
  layout.pos.col = 1
))

print(
  zoom_map,
  vp = viewport(
    x = 0.5,
    y = 0.5,
    width = 0.95,
    height = 0.45
  )
)

popViewport()

# Three main maps
main_maps <- list(main_map, m8, m9)

for (i in 1:3) {
  pushViewport(viewport(
    layout.pos.row = 1,
    layout.pos.col = i + 1
  ))
  
  print(main_maps[[i]], vp = viewport())
  
  popViewport()
}

popViewport()

summary(pfna_0612$avg_pfna)
summary(pfna_1318$avg_pfna)
summary(pfna_1925$avg_pfna)
