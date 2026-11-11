######################################################
# R Programming for XGBoost - Functions
# Created by: Hari Iyer
# Date: 3/19/26

# Adding IDW measures


######################################################

library(tidyverse)
## install.packages("caret")
library(caret)
library(xgboost)
library(pROC)

library(tidyverse)
library(sf)
## install.packages("terra")
library(terra)
## install.packages("xgboost")
library(xgboost)
## install.packages("DiagrammeR")
library(DiagrammeR)
## install.packages("pROC")
library(pROC)
## install.packages("caret")
library(caret)
## install.packages("ModelMetrics")
## install.packages("vtreat")
library(ModelMetrics)
library(vtreat)
## install.packages("MLmetrics")
library(MLmetrics)


## read in data

load("~/R/XGBoost/data/pfas_xg.Rda")
nj_cws_pfas_dist <- st_read("~/R/XGBoost/data/nj_cws_pfas_dist.shp")
load("~/R/XGBoost/data/cws_aw_means.Rda")
load("~/R/XGBoost/data/Pollution/njpws_pm25_joined.Rda")
load("~/R/XGBoost/data/Pollution/njpws_no2_joined.Rda")
load("~/R/XGBoost/data/tok_cws_joined.Rda")
nj_mirta <- st_read("~/R/XGBoost/data/nj_cws_mirta_dist.shp")
nj_fuds <- st_read("~/R/XGBoost/data/nj_cws_fuds_dist.shp")
nj_airport <- st_read("~/R/XGBoost/data/nj_cws_airp_dist.shp")
load("~/R/XGBoost/data/nj_cws_idw.Rda")

length(unique(cws_aw_means$PWSID)) 
length(unique(nj_cws_pfas_dist$PWSID)) 
length(unique(njpws_pm25_joined$PWSID)) 
length(unique(njpws_no2_joined$PWSID))
length(unique(tok_cws_joined$PWSID)) 
length(unique(nj_mirta$PWSID)) 
length(unique(nj_fuds$PWSID)) 
length(unique(nj_airport$PWSID)) 
length(unique(nj_cws_idw$PWSID)) 


length(unique(pfas_full$PWSID))

# filter out NTNCs from pfas full
pfas_full <- pfas_full |>
  filter(!(PWSID %in% c("NJ0105360", "NJ0108352", "NJ0329303", "NJ1351327", "NJ1924326", "NJ0102301", "NJ0102302", "NJ0111304", "NJ0112325", "NJ0808319", "NJ0817305", "NJ0820301", "NJ1106316", "NJ1106319", "NJ1332313", "NJ1406325", "NJ1511415", "NJ1904357", "NJ1908322", "NJ1910337", "NJ1922350", "NJ1924316" )))

# set up two inactive water systems

pfas_full <- pfas_full |>
  mutate(PWSID = case_when(PWSID== "NJ1003001" ~ "NJ2119001",
                           PWSID== "NJ1511300" ~ "NJ1513001",
                           TRUE ~ PWSID))

## clean and merge to CWS

nj_pfas_indust <- nj_cws_pfas_dist %>%
  st_drop_geometry() 

nj_fuds <- nj_fuds %>%
  st_drop_geometry()

nj_mirta <- nj_mirta %>%
  st_drop_geometry()

nj_airport <- nj_airport %>%
  st_drop_geometry()



njpws_no2 <- njpws_no2_joined %>%
  dplyr::select(PWSID, Pop_Cat_5, NO2_Annual_2000, NO2_Annual_2001, NO2_Annual_2002,
                NO2_Annual_2003, NO2_Annual_2004, NO2_Annual_2005, NO2_Annual_2006,
                NO2_Annual_2007, NO2_Annual_2008, NO2_Annual_2009, NO2_Annual_2010,
                NO2_Annual_2011, NO2_Annual_2012, NO2_Annual_2013, NO2_Annual_2014,
                NO2_Annual_2015, NO2_Annual_2016) %>%
  st_drop_geometry()

njpws_pm25 <- njpws_pm25_joined %>%
  dplyr::select(PWSID, US_pm25_2000, US_pm25_2001, US_pm25_2002,
                US_pm25_2003, US_pm25_2004, US_pm25_2005, US_pm25_2006,
                US_pm25_2007, US_pm25_2008, US_pm25_2009, US_pm25_2010,
                US_pm25_2011, US_pm25_2012, US_pm25_2013, US_pm25_2014,
                US_pm25_2015, US_pm25_2016) %>%
  st_drop_geometry()


table(nj_cws_idw$year) ## need to merge based on prior year

table(pfas_full$year)

year_map <- data.frame(
  year = 2006:2025
) %>%
  mutate(
    year_group = case_when(
      year >= 2006 & year <= 2010 ~ 2005,
      year >= 2011 & year <= 2015 ~ 2010,
      year >= 2016 & year <= 2020 ~ 2015,
      year >= 2020 ~ 2020
    )
  )

pfas_full <- pfas_full %>%
  left_join(year_map, by = "year")



nj_prep_xg <- pfas_full %>%
  left_join(nj_pfas_indust, by = "PWSID") %>%
  left_join(cws_aw_means, by = "PWSID") %>%
  left_join(njpws_pm25, by = "PWSID") %>%
  left_join(njpws_no2, by = "PWSID") %>%
  left_join(tok_cws_joined, by = "PWSID") %>%
  left_join(nj_fuds, by = "PWSID") %>%
  left_join(nj_mirta, by = "PWSID") %>%
  left_join(nj_airport, by = "PWSID") %>%
  left_join(nj_cws_idw, by = c("PWSID" = "PWSID","year_group" = "year"))

summary(nj_prep_xg)

check_prep <- nj_prep_xg |>
  filter(is.na(dist_pfas_sites))

missing_id_idw <- pfas_full |>
  filter(!(PWSID %in% nj_cws_idw$PWSID))


length(unique(missing_id_idw$PWSID))

missing_PWSID <- missing_id_idw |>
  distinct(PWSID, .keep_all = TRUE) |>
  select(PWSID)


missing_id_no2 <- pfas_full |>
  filter(!(PWSID %in% njpws_no2$PWSID))


length(unique(missing_id_no2$PWSID))



rm(nj_cws_pfas_dist, nj_pfas_indust, cws_aw_means,
   njpws_no2, njpws_no2_joined, njpws_pm25, njpws_pm25_joined,
   tok_cws_joined, nj_fuds, nj_mirta, nj_airport, nj_cws_idw)

nj_prep_xg$pm25 <- ifelse(nj_prep_xg$year == 2006, nj_prep_xg$US_pm25_2006,
                          ifelse(nj_prep_xg$year == 2007,nj_prep_xg$US_pm25_2007,
                                 ifelse(nj_prep_xg$year == 2008,nj_prep_xg$US_pm25_2008,
                                        ifelse(nj_prep_xg$year == 2009,nj_prep_xg$US_pm25_2009,
                                               ifelse(nj_prep_xg$year == 2010,nj_prep_xg$US_pm25_2010,
                                                      ifelse(nj_prep_xg$year == 2011,nj_prep_xg$US_pm25_2011,
                                                             ifelse(nj_prep_xg$year == 2012,nj_prep_xg$US_pm25_2012,
                                                                    ifelse(nj_prep_xg$year == 2013,nj_prep_xg$US_pm25_2013,
                                                                           ifelse(nj_prep_xg$year == 2014,nj_prep_xg$US_pm25_2014,
                                                                                  ifelse(nj_prep_xg$year == 2015,nj_prep_xg$US_pm25_2015,nj_prep_xg$US_pm25_2016))))))))))


nj_prep_xg$no2 <- ifelse(nj_prep_xg$year == 2006, nj_prep_xg$NO2_Annual_2006,
                         ifelse(nj_prep_xg$year == 2007,nj_prep_xg$NO2_Annual_2007,
                                ifelse(nj_prep_xg$year == 2008,nj_prep_xg$NO2_Annual_2008,
                                       ifelse(nj_prep_xg$year == 2009,nj_prep_xg$NO2_Annual_2009,
                                              ifelse(nj_prep_xg$year == 2010,nj_prep_xg$NO2_Annual_2010,
                                                     ifelse(nj_prep_xg$year == 2011,nj_prep_xg$NO2_Annual_2011,
                                                            ifelse(nj_prep_xg$year == 2012,nj_prep_xg$NO2_Annual_2012,
                                                                   ifelse(nj_prep_xg$year == 2013,nj_prep_xg$NO2_Annual_2013,
                                                                          ifelse(nj_prep_xg$year == 2014,nj_prep_xg$NO2_Annual_2014,
                                                                                 ifelse(nj_prep_xg$year == 2015,nj_prep_xg$NO2_Annual_2015,nj_prep_xg$NO2_Annual_2016))))))))))

nj_prep_xg <- nj_prep_xg %>%
  dplyr::select(-c(NO2_Annual_2000, NO2_Annual_2001, NO2_Annual_2002,
                   NO2_Annual_2003, NO2_Annual_2004, NO2_Annual_2005, NO2_Annual_2006,
                   NO2_Annual_2007, NO2_Annual_2008, NO2_Annual_2009, NO2_Annual_2010,
                   NO2_Annual_2011, NO2_Annual_2012, NO2_Annual_2013, NO2_Annual_2014,
                   NO2_Annual_2015, NO2_Annual_2016, US_pm25_2000, US_pm25_2001, US_pm25_2002,
                   US_pm25_2003, US_pm25_2004, US_pm25_2005, US_pm25_2006,
                   US_pm25_2007, US_pm25_2008, US_pm25_2009, US_pm25_2010,
                   US_pm25_2011, US_pm25_2012, US_pm25_2013, US_pm25_2014,
                   US_pm25_2015, US_pm25_2016))

nj_prep_xg$lnres <- log(nj_prep_xg$RESULT_NG + 0.01)

summary(nj_prep_xg)
colnames(nj_prep_xg)
# [1] "PWSID"             "pfas_name"         "pfas_label"        "pfas_label_num"    "season"            "year"             
# [7] "RESULT_NG"         "year_group"        "aw_mn__"           "aw_edult810"       "aw_edugt8lths10"   "aw_eduhs10"       
# [13] "aw_edubach10"      "aw_mdhval10"       "aw_mdinc10"        "aw_percwhite10"    "aw_percblack10"    "aw_percforeign10" 
# [19] "aw_percpov10"      "aw_pohse10"        "aw_pfwchbf10"      "aw_unemp10"        "aw_munemp10"       "aw_percfamint10"  
# [25] "aw_percwelfr10"    "aw_percsales10"    "aw_percclerical10" "aw_percrepair10"   "aw_perctransp10"   "aw_percnflabor10" 
# [31] "aw_percservice10"  "aw_percchild10"    "aw_percold10"      "aw_mdrent10"       "aw_incice10"       "aw_rice10"        
# [37] "aw_rcincice10"     "Pop_Cat_5"         "avg_clay"          "depth_m"           "dist_pfas_sites"   "dtw_m"            
# [43] "n_sep2010"         "nwalt2002_pcturb"  "popdensity2000"    "rech48_mmyr"       "awfuds"            "awmirta"          
# [49] "awairport"         "IDW_03_05"         "IDW_03_10"         "IDW_05_05"         "IDW_05_10"         "IDW_05_15"        
# [55] "pp_dns_"           "IW_FIRE_5K_5"      "IW_FIRE_5K_10"     "IW_FIRE_5K_15"     "IW_FIRE_3K_5"      "IW_FIRE_3K_10"    
# [61] "IW_FIRE_1K_5"      "IW_AP_10K_5"       "IW_AP_5K_5"        "IW_AP_3K_5"        "IW_AP_1K_5"        "IW_FUDS_10K_10"   
# [67] "IW_FUDS_5K_5"      "IW_FUDS_3K_5"      "IW_FUDS_1K_5"      "IW_MIRTA_10K_5"    "IW_MIRTA_5K_5"     "IW_MIRTA_3K_5"    
# [73] "pm25"              "no2"               "lnres"       



## CREATE BINARY LABEL - DETECT VS NON-DETECT

nj_prep_xg <- nj_prep_xg %>%
  mutate(pfas_bin = case_when(
    pfas_label %in% c("<MCL", ">=MCL", "<20ng/L", ">=20ng/L",">=40ng/L") ~ "DET",
    pfas_label == "ND" ~ "ND",
    is.na(pfas_label) ~ "Missing",
    TRUE ~ "Missing"
  ),
  pfas_bin_num = case_when(
    pfas_bin == "DET" ~ 1,
    pfas_bin == "ND" ~ 0,
  ))


table(nj_prep_xg$pfas_bin)
table(nj_prep_xg$pfas_label)
table(nj_prep_xg$pfas_bin_num)

## missing category for Pop_Cat_5

nj_prep_xg[is.na(nj_prep_xg$Pop_Cat_5),"Pop_Cat_5"] <- "Miss"

# impute median

na_cols <- names(colSums(is.na(nj_prep_xg)) > 0)[colSums(is.na(nj_prep_xg)) > 0]
names(nj_prep_xg)

na_cols <- c("season", "year","year_group","aw_mn__","aw_edult810", "aw_edugt8lths10", "aw_eduhs10",       
"aw_edubach10",      "aw_mdhval10",       "aw_mdinc10",        "aw_percwhite10",     "aw_percblack10",    "aw_percforeign10", 
"aw_percpov10",      "aw_pohse10",        "aw_pfwchbf10",      "aw_unemp10",        "aw_munemp10",       "aw_percfamint10",  
"aw_percwelfr10",    "aw_percsales10",    "aw_percclerical10", "aw_percrepair10",   "aw_perctransp10",   "aw_percnflabor10", 
"aw_percservice10",  "aw_percchild10",    "aw_percold10",      "aw_mdrent10",       "aw_incice10",       "aw_rice10",        
"aw_rcincice10",     "Pop_Cat_5",         "avg_clay",          "depth_m",           "dist_pfas_sites",   "dtw_m",            
"n_sep2010",         "nwalt2002_pcturb",  "popdensity2000",    "rech48_mmyr",       "awfuds",            "awmirta",          
"awairport",         "IDW_03_05",         "IDW_03_10",         "IDW_05_05",         "IDW_05_10",         "IDW_05_15",        
"pp_dns_",           "IW_FIRE_5K_5",      "IW_FIRE_5K_10",     "IW_FIRE_5K_15",     "IW_FIRE_3K_5",      "IW_FIRE_3K_10",    
"IW_FIRE_1K_5",      "IW_AP_10K_5",       "IW_AP_5K_5",        "IW_AP_3K_5",        "IW_AP_1K_5",        "IW_FUDS_10K_10",   
"IW_FUDS_5K_5",      "IW_FUDS_3K_5",      "IW_FUDS_1K_5",      "IW_MIRTA_10K_5",    "IW_MIRTA_5K_5",     "IW_MIRTA_3K_5",    
"pm25",              "no2")

summary(nj_prep_xg)

nj_prep_xg <- nj_prep_xg %>%
  mutate(across(
    all_of(na_cols),
    ~ replace(., is.na(.), median(., na.rm = TRUE))
  ))

nj_prep_xg_pfoa <- nj_prep_xg %>%
  filter(pfas_name == "PFOA") 

nj_prep_xg_pfos <- nj_prep_xg %>%
  filter(pfas_name == "PFOS")

nj_prep_xg_pfna <- nj_prep_xg %>%
  filter(pfas_name == "PFNA")

summary(nj_prep_xg_pfoa)

################# saving for test/train sensitivity analyses ####################

save(
  nj_prep_xg_pfoa,
  file = "~/R/registry/nj_prep_xg_pfoa.Rda"
)

save(
  nj_prep_xg_pfos,
  file = "~/R/registry/nj_prep_xg_pfos.Rda"
)

save(
  nj_prep_xg_pfna,
  file = "~/R/registry/nj_prep_xg_pfna.Rda"
)





#######################################
# stats for abstract
#######################################

table(nj_prep_xg_pfoa$pfas_bin_num)
table(nj_prep_xg_pfos$pfas_bin_num)
table(nj_prep_xg_pfna$pfas_bin_num)

## only in detects

pfoa_dets <- nj_prep_xg_pfoa %>%
  filter(pfas_bin_num == 1)

pfos_dets <- nj_prep_xg_pfos %>%
  filter(pfas_bin_num == 1)

pfna_dets <- nj_prep_xg_pfna %>%
  filter(pfas_bin_num == 1)

summary(pfoa_dets$RESULT_NG)
summary(pfos_dets$RESULT_NG)
summary(pfna_dets$RESULT_NG)


library(dplyr)
library(caret)
library(xgboost)
library(pROC)
library(ggplot2)

#==============================================================
# Utility metrics
#==============================================================

r2_score <- function(y_true, y_pred) {
  ss_res <- sum((y_true - y_pred)^2, na.rm = TRUE)
  ss_tot <- sum((y_true - mean(y_true, na.rm = TRUE))^2, na.rm = TRUE)
  1 - ss_res / ss_tot
}

rmse_score <- function(y_true, y_pred) {
  sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
}

#==============================================================
# Preprocessing: fit dummy encoder on training only
#==============================================================

fit_encoder <- function(df,
                        drop_vars = c("pfas_name","pfas_bin","year","year_group",
                                      "pfas_label_num","pfas_label",
                                      "pfas_bin_num","lnres","RESULT_NG"),
                        cat_vars = c("Pop_Cat_5"),
                        id_vars = c("PWSID")) {
  
  df2 <- df
  
  for (v in cat_vars) {
    df2[[v]] <- as.character(df2[[v]])
    df2[[v]][is.na(df2[[v]])] <- "Miss"
    df2[[v]] <- factor(df2[[v]])
  }
  
  x_df <- df2 %>%
    dplyr::select(-dplyr::any_of(drop_vars))
  
  x_input <- x_df %>%
    dplyr::select(-dplyr::any_of(id_vars))
  
  dv <- caret::dummyVars(~ ., data = x_input, fullRank = FALSE)
  
  list(
    dv = dv,
    drop_vars = drop_vars,
    cat_vars = cat_vars,
    id_vars = id_vars
  )
}

apply_encoder <- function(df, encoder) {
  df2 <- df
  
  for (v in encoder$cat_vars) {
    df2[[v]] <- as.character(df2[[v]])
    df2[[v]][is.na(df2[[v]])] <- "Miss"
    df2[[v]] <- factor(df2[[v]])
  }
  
  x_df <- df2 %>%
    dplyr::select(-dplyr::any_of(encoder$drop_vars))
  
  x_input <- x_df %>%
    dplyr::select(-dplyr::any_of(encoder$id_vars))
  
  x_mat <- predict(encoder$dv, newdata = x_input)
  x_mat <- as.matrix(x_mat)
  
  colnames(x_mat) <- make.names(colnames(x_mat), unique = TRUE)
  x_mat
}

#==============================================================
# Hyperparameter tuning: binary classification
#==============================================================

tune_xgb_binary <- function(X_train, y_train, tune_grid,
                            nfold = 5,
                            early_stopping_rounds = 25,
                            seed = 123) {
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  results <- vector("list", nrow(tune_grid))
  
  set.seed(seed)
  
  for (i in seq_len(nrow(tune_grid))) {
    
    params <- list(
      booster = "gbtree",
      objective = "binary:logistic",
      eval_metric = c("auc", "logloss"),
      eta = tune_grid$eta[i],
      max_depth = tune_grid$max_depth[i],
      gamma = tune_grid$gamma[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      min_child_weight = tune_grid$min_child_weight[i],
      subsample = tune_grid$subsample[i],
      lambda = tune_grid$lambda[i],
      alpha = tune_grid$alpha[i]
    )
    
    cv_fit <- tryCatch(
      xgb.cv(
        params = params,
        data = dtrain,
        nrounds = tune_grid$nrounds[i],
        nfold = nfold,
        early_stopping_rounds = early_stopping_rounds,
        verbose = 0
      ),
      error = function(e) NULL
    )
    
    if (is.null(cv_fit)) {
      results[[i]] <- cbind(tune_grid[i, , drop = FALSE],
                            best_nrounds = NA_integer_,
                            cv_auc = NA_real_,
                            cv_logloss = NA_real_)
      next
    }
    
    ev <- cv_fit$evaluation_log
    best_iter <- cv_fit$best_iteration
    
    auc_col <- "test_auc_mean"
    logloss_col <- "test_logloss_mean"
    
    if (is.null(best_iter) || length(best_iter) == 0 || is.na(best_iter)) {
      if (logloss_col %in% names(ev)) {
        best_iter <- which.min(ev[[logloss_col]])
      } else {
        best_iter <- nrow(ev)
      }
    }
    
    best_auc <- if (auc_col %in% names(ev)) ev[[auc_col]][best_iter] else NA_real_
    best_logloss <- if (logloss_col %in% names(ev)) ev[[logloss_col]][best_iter] else NA_real_
    
    results[[i]] <- cbind(tune_grid[i, , drop = FALSE],
                          best_nrounds = best_iter,
                          cv_auc = best_auc,
                          cv_logloss = best_logloss)
  }
  
  results_df <- dplyr::bind_rows(results) %>%
    dplyr::arrange(cv_logloss, dplyr::desc(cv_auc))
  
  results_df
}

#==============================================================
# Hyperparameter tuning: regression among detects
#==============================================================

tune_xgb_regression <- function(X_train, y_train, tune_grid,
                                nfold = 5,
                                early_stopping_rounds = 25,
                                seed = 123) {
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  results <- vector("list", nrow(tune_grid))
  
  set.seed(seed)
  
  for (i in seq_len(nrow(tune_grid))) {
    
    params <- list(
      booster = "gbtree",
      objective = "reg:squarederror",
      eval_metric = "rmse",
      eta = tune_grid$eta[i],
      max_depth = tune_grid$max_depth[i],
      gamma = tune_grid$gamma[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      min_child_weight = tune_grid$min_child_weight[i],
      subsample = tune_grid$subsample[i],
      lambda = tune_grid$lambda[i],
      alpha = tune_grid$alpha[i]
    )
    
    cv_fit <- tryCatch(
      xgb.cv(
        params = params,
        data = dtrain,
        nrounds = tune_grid$nrounds[i],
        nfold = nfold,
        early_stopping_rounds = early_stopping_rounds,
        verbose = 0
      ),
      error = function(e) NULL
    )
    
    if (is.null(cv_fit)) {
      results[[i]] <- cbind(tune_grid[i, , drop = FALSE],
                            best_nrounds = NA_integer_,
                            cv_rmse = NA_real_)
      next
    }
    
    ev <- cv_fit$evaluation_log
    best_iter <- cv_fit$best_iteration
    rmse_col <- "test_rmse_mean"
    
    if (is.null(best_iter) || length(best_iter) == 0 || is.na(best_iter)) {
      best_iter <- if (rmse_col %in% names(ev)) which.min(ev[[rmse_col]]) else nrow(ev)
    }
    
    best_rmse <- if (rmse_col %in% names(ev)) ev[[rmse_col]][best_iter] else NA_real_
    
    results[[i]] <- cbind(tune_grid[i, , drop = FALSE],
                          best_nrounds = best_iter,
                          cv_rmse = best_rmse)
  }
  
  results_df <- dplyr::bind_rows(results) %>%
    dplyr::arrange(cv_rmse)
  
  results_df
}

#==============================================================
# Fit final binary model with watchlist
#==============================================================

fit_final_binary <- function(X_train, y_train, X_test, y_test, best_row) {
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest  <- xgb.DMatrix(data = X_test, label = y_test)
  
  params <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = c("auc", "logloss"),
    eta = best_row$eta,
    max_depth = best_row$max_depth,
    gamma = best_row$gamma,
    colsample_bytree = best_row$colsample_bytree,
    min_child_weight = best_row$min_child_weight,
    subsample = best_row$subsample,
    lambda = best_row$lambda,
    alpha = best_row$alpha
  )
  
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_row$best_nrounds,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )
  
  p_train <- predict(model, newdata = dtrain)
  p_test  <- predict(model, newdata = dtest)
  
  brier_train <- mean((p_train - y_train)^2)
  brier_test <- mean((p_test - y_test)^2)
  
  auc_train <- as.numeric(pROC::auc(pROC::roc(y_train, p_train, quiet = TRUE)))
  auc_test  <- as.numeric(pROC::auc(pROC::roc(y_test, p_test, quiet = TRUE)))
  
  out <- list(
    model = model,
    preds_train = p_train,
    preds_test = p_test,
    auc_train = auc_train,
    auc_test = auc_test,
    brier_train = brier_train,
    brier_test = brier_test,
    evaluation_log = model$evaluation_log,
    params = params
  )
  return(out)
}

#==============================================================
# Fit final regression model with watchlist
#==============================================================

fit_final_regression <- function(X_train, y_train, X_test, y_test, best_row) {
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest  <- xgb.DMatrix(data = X_test, label = y_test)
  
  params <- list(
    booster = "gbtree",
    objective = "reg:squarederror",
    eval_metric = "rmse",
    eta = best_row$eta,
    max_depth = best_row$max_depth,
    gamma = best_row$gamma,
    colsample_bytree = best_row$colsample_bytree,
    min_child_weight = best_row$min_child_weight,
    subsample = best_row$subsample,
    lambda = best_row$lambda,
    alpha = best_row$alpha
  )
  
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_row$best_nrounds,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )
  
  pred_train <- predict(model, newdata = dtrain)
  pred_test  <- predict(model, newdata = dtest)
  
  rmse_train <- rmse_score(y_train, pred_train)
  rmse_test  <- rmse_score(y_test, pred_test)
  
  r2_train <- r2_score(y_train, pred_train)
  r2_test  <- r2_score(y_test, pred_test)
  
  out <- list(
    model = model,
    preds_train = pred_train,
    preds_test = pred_test,
    rmse_train = rmse_train,
    rmse_test = rmse_test,
    r2_train = r2_train,
    r2_test = r2_test,
    evaluation_log = model$evaluation_log,
    params = params
  )
  
  return(out)
}

#==============================================================
# Plot stage 1 AUC curves
#==============================================================

plot_auc_curves <- function(eval_log) {
  
  ggplot(eval_log, aes(x = iter)) +
    geom_line(aes(y = train_auc, color = "Train AUC")) +
    geom_line(aes(y = test_auc, color = "Test AUC")) +
    labs(
      x = "Boosting Iteration",
      y = "AUC",
      color = NULL,
      title = "Stage 1: Training and Test AUC"
    ) +
    theme_minimal()
}

#==============================================================
# Main wrapper
#==============================================================

run_two_phase_xgb <- function(df,
                              cont_outcome = "lnres",
                              bin_var = "pfas_bin_num",
                              det_value = 1,
                              seed = 123,
                              train_prop = 0.7,
                              tune_grid_bin,
                              tune_grid_reg,
                              nfold = 5,
                              early_stopping_rounds = 25) {
  
  set.seed(seed)
  n <- nrow(df)
  idx_train <- sample.int(n, size = round(train_prop * n), replace = FALSE)
  in_train <- rep(FALSE, n)
  in_train[idx_train] <- TRUE
  
  df_train <- df[in_train, , drop = FALSE]
  df_test  <- df[!in_train, , drop = FALSE]
  
  #------------------------------------------------------------
  # Encode predictors using training only
  #------------------------------------------------------------
  encoder <- fit_encoder(df_train)
  
  X_train <- apply_encoder(df_train, encoder)
  X_test  <- apply_encoder(df_test, encoder)
  
  #------------------------------------------------------------
  # Stage 1: detect vs non-detect
  #------------------------------------------------------------
  y_bin_train <- ifelse(df_train[[bin_var]] == det_value, 1, 0)
  y_bin_test  <- ifelse(df_test[[bin_var]] == det_value, 1, 0)
  
  bin_cv <- tune_xgb_binary(
    X_train = X_train,
    y_train = y_bin_train,
    tune_grid = tune_grid_bin,
    nfold = nfold,
    early_stopping_rounds = early_stopping_rounds,
    seed = seed
  )
  
  best_bin <- bin_cv[1, , drop = FALSE]
  
  bin_fit <- fit_final_binary(
    X_train = X_train,
    y_train = y_bin_train,
    X_test = X_test,
    y_test = y_bin_test,
    best_row = best_bin
  )
  
  auc_plot <- plot_auc_curves(bin_fit$evaluation_log)
  
  #------------------------------------------------------------
  # Stage 2: concentration among detects only
  #------------------------------------------------------------
  train_detect <- df_train[[bin_var]] == det_value & !is.na(df_train[[cont_outcome]])
  test_detect  <- df_test[[bin_var]] == det_value & !is.na(df_test[[cont_outcome]])
  
  X_train_det <- X_train[train_detect, , drop = FALSE]
  X_test_det  <- X_test[test_detect, , drop = FALSE]
  
  y_train_det <- df_train[[cont_outcome]][train_detect]
  y_test_det  <- df_test[[cont_outcome]][test_detect]
  
  reg_cv <- tune_xgb_regression(
    X_train = X_train_det,
    y_train = y_train_det,
    tune_grid = tune_grid_reg,
    nfold = nfold,
    early_stopping_rounds = early_stopping_rounds,
    seed = seed
  )
  
  best_reg <- reg_cv[1, , drop = FALSE]
  
  reg_fit <- fit_final_regression(
    X_train = X_train_det,
    y_train = y_train_det,
    X_test = X_test_det,
    y_test = y_test_det,
    best_row = best_reg
  )
  
  list(
    idx_train = idx_train,
    encoder = encoder,
    
    stage1 = list(
      cv_results = bin_cv,
      best_hyperparameters = best_bin,
      model = bin_fit$model,
      auc_train = bin_fit$auc_train,
      auc_test = bin_fit$auc_test,
      train_auc = bin_fit$evaluation_log$train_auc,
      test_auc  = bin_fit$evaluation_log$test_auc, 
      min_cv_logloss = best_bin$cv_logloss,
      evaluation_log = bin_fit$evaluation_log,
      brier_train = bin_fit$brier_train,
      brier_test = bin_fit$brier_test,
      auc_plot = auc_plot
    ),
    
    stage2 = list(
      cv_results = reg_cv,
      best_hyperparameters = best_reg,
      model = reg_fit$model,
      rmse_train = reg_fit$rmse_train,
      rmse_test = reg_fit$rmse_test,
      r2_train = reg_fit$r2_train,
      r2_test = reg_fit$r2_test,
      evaluation_log = reg_fit$evaluation_log
    )
  )
}

## start by tuning max_deptth, min_child_weight, eta, rounds
## Identify best parameters
## then move to subsample, colsample_bytree
## then move to gamma, lambda, alpha

summary(nj_prep_xg_pfoa)

###############################################
# PFOA
###############################################

tune_grid_bin_pfoa <- expand.grid(
  nrounds = 400,
  eta = 0.08,
  gamma = 0.05,
  colsample_bytree = 0.7,
  min_child_weight = 2,
  subsample = 1,
  lambda = 1,
  alpha = 1,
  max_depth = 4
)


tune_grid_reg_pfoa <- expand.grid(
  nrounds = 300,
  eta = 0.03,
  max_depth = 7,
  gamma = 0.05,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  lambda = 0.5,
  alpha = 0.5
)


res_pfoa <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfoa,
  tune_grid_reg = tune_grid_reg_pfoa,
  nfold = 5,
  early_stopping_rounds = 25
)


res_pfoa$stage1$auc_train ## 0.938266
res_pfoa$stage1$auc_test ## 0.9112909
res_pfoa$stage1$min_cv_logloss ## 0.3759071
res_pfoa$stage1$best_hyperparameters
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
#     400 0.08  0.05              0.7                2         1      1     1         4          350 0.9114539  0.3759071

res_pfoa$stage1$brier_train ## 0.1008656
res_pfoa$stage1$brier_test ## 0.1207971

#this doesn't print - double check the code for this - maybe it's above - could be a train_auc vs auc_train error
print(res_pfoa$stage1$auc_plot)

# Phase 2 results

res_pfoa$stage2$rmse_train ## 0.3812723
res_pfoa$stage2$rmse_test ## 0.425658
res_pfoa$stage2$r2_train ## 0.7536007
res_pfoa$stage2$r2_test ## 0.6838884
res_pfoa$stage2$best_hyperparameters 
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#     300 0.03         7  0.05              0.8                1       0.9    0.5   0.5          142 0.4344264

## check importance

names(res_pfoa$stage1)
names(res_pfoa$stage2)
str(res_pfoa$stage1, max.level = 1)
str(res_pfoa$stage2, max.level = 1)

# Stage 1: binary classification model
imp_stage1 <- xgb.importance(model = res_pfoa$stage1$model)
print(imp_stage1)
xgb.plot.importance(imp_stage1)


# Stage 2: regression model
imp_stage2 <- xgb.importance(model = res_pfoa$stage2$model)
print(imp_stage2)
xgb.plot.importance(imp_stage2)


save(res_pfoa, file="/home/ms3283/R/XGBoost/models/res_pfoa_4.30.26.Rda")
load("/home/ms3283/R/XGBoost/models/res_pfoa_4.30.26.Rda")



#####################################REDO

#######################################################################
# PFOS
#######################################################################

## start by tuning max_depth, min_child_weight, eta, rounds
## Identify best parameters
## then move to subsample, colsample_bytree
## then move to gamma, lambda, alpha

###############################################
# PFOS 
###############################################

tune_grid_bin_pfos <- expand.grid(
  nrounds = 1000,
  eta = 0.05,
  max_depth = 7,
  gamma = 0,
  colsample_bytree = 0.6,
  min_child_weight = 2,
  subsample = 1,
  lambda = 1,
  alpha = 0
)


tune_grid_reg_pfos <- expand.grid(
  nrounds = 1000,
  eta = 0.04,
  max_depth = 6,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1,
  lambda = 0,
  alpha = 0
)

res_pfos <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos,
  nfold = 5,
  early_stopping_rounds = 25
)



##################### final test for PFOS ##########

res_pfos$stage1$auc_train ## 0.9362747
res_pfos$stage1$auc_test ## 0.9009997
res_pfos$stage1$min_cv_logloss ## 0.3525931
res_pfos$stage1$best_hyperparameters 
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.05         7     0              0.6                2         1      1     0          213 0.9113486  0.3525931
res_pfos$stage1$brier_train ## 0.09787134
res_pfos$stage1$brier_test ## 0.1191432

#print(res_pfos$stage1$auc_plot)

res_pfos$stage2$rmse_train ## 0.4026584
res_pfos$stage2$rmse_test ## 0.4687772
res_pfos$stage2$r2_train ## 0.6345787
res_pfos$stage2$r2_test ## 0.5272821
res_pfos$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#    1000 0.04         6   0.1                1                1         1      0     0           95 0.4562135


## check importance

names(res_pfos$stage1)
names(res_pfos$stage2)
str(res_pfos$stage1, max.level = 1)
str(res_pfos$stage2, max.level = 1)

# Stage 1: binary classification model
pfos.imp_stage1 <- xgb.importance(model = res_pfos$stage1$model)
print(pfos.imp_stage1)
xgb.plot.importance(pfos.imp_stage1)


# Stage 2: regression model
pfos.imp_stage2 <- xgb.importance(model = res_pfos$stage2$model)
print(pfos.imp_stage2)
xgb.plot.importance(pfos.imp_stage2)


save(res_pfos, file="/home/ms3283/R/XGBoost/models/res_pfos_4.30.26.Rda")
saveRDS(res_pfoa, "res_pfoa.rds")


#######################################################################
# PFNA
#######################################################################

tune_grid_bin_pfna <- expand.grid(
  nrounds = 1000,
  eta = 0.04,
  max_depth = 4,
  gamma = 0.1,
  colsample_bytree = 0.8,
  min_child_weight = 3,
  subsample = 0.9,
  lambda = 1,
  alpha = 0
)


tune_grid_reg_pfna <- expand.grid(
  nrounds = 1000,
  eta = 0.03,
  max_depth = 6,
  gamma = 0.05,
  colsample_bytree = 1,
  min_child_weight = 4,
  subsample = 0.8,
  lambda = 2,
  alpha = 0.5
)

res_pfna <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,   # fixed binary
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5,
  early_stopping_rounds = 25
)



res_pfna$stage1$auc_train ## 0.9481301
res_pfna$stage1$auc_test ## 0.9166623
res_pfna$stage1$min_cv_logloss ## 0.1572805
res_pfna$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.04         4   0.1              0.8                3       0.9      1     0          279 0.9238612  0.1572805

res_pfna$stage1$brier_train ## 0.04068846
res_pfna$stage1$brier_test ## 0.04333066

#print(res_pfna$stage1$auc_plot)

res_pfna$stage2$rmse_train ## 0.441681
res_pfna$stage2$rmse_test ## 0.560773
res_pfna$stage2$r2_train ## 0.8581829
res_pfna$stage2$r2_test ## 0.7464395
res_pfna$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#    1000 0.03         6  0.05                1                4       0.8      2   0.5          159 0.5781429


## check importance

names(res_pfna$stage1)
names(res_pfna$stage2)
str(res_pfna$stage1, max.level = 1)
str(res_pfna$stage2, max.level = 1)

# Stage 1: binary classification model
pfna.imp_stage1 <- xgb.importance(model = res_pfna$stage1$model)
print(pfna.imp_stage1)
xgb.plot.importance(pfna.imp_stage1)


# Stage 2: regression model
pfna.imp_stage2 <- xgb.importance(model = res_pfna$stage2$model)
print(pfna.imp_stage2)
xgb.plot.importance(pfna.imp_stage2)

save(res_pfna, file="/home/ms3283/R/XGBoost/models/res_pfna_4.30.26.Rda")
head(res_pfna$stage1$evaluation_log)


############### saving final models as .rds #####################
saveRDS(res_pfoa, "res_pfoa.rds")
saveRDS(res_pfos, "res_pfos.rds")
saveRDS(res_pfna, "res_pfna.rds")

################### recreating feature importance plot #####################

load("~/R/XGBoost/models/res_pfoa.Rda")
load("~/R/XGBoost/models/res_pfos.Rda")
load("~/R/XGBoost/models/res_pfna.Rda")

# PFOA
pfoa_imp_stage1 <- xgb.importance(
  model = res_pfoa$stage1$model
)

pfoa_imp_stage2 <- xgb.importance(
  model = res_pfoa$stage2$model
)

# PFOS
pfos_imp_stage1 <- xgb.importance(
  model = res_pfos$stage1$model
)

pfos_imp_stage2 <- xgb.importance(
  model = res_pfos$stage2$model
)

# PFNA
pfna_imp_stage1 <- xgb.importance(
  model = res_pfna$stage1$model
)

pfna_imp_stage2 <- xgb.importance(
  model = res_pfna$stage2$model
)

# ============================================================
# 1. Create lookup: old variable name = new display name
# ============================================================

variable_labels <- c(
  "pm25"                   = "PM2.5",
  "aw_percforeign10"       = "Percent foreign born",
  "dtw_m"                  = "Depth to water",
  "no2"                    = "NO2",
  "rech48_mmyr"            = "Groundwater recharge",
  "aw_rcincice10"          = "Income - Race ICE",
  "aw_pohse10"             = "Occupied Housing",
  "aw_mdhval10"            = "Median home value",
  "pp_dns_"                = "Population density",
  "aw_percnflabor10"       = "Percent female labor",
  "avg_clay"               = "Average clay",
  "aw_mn__"                = "PFAS industrial contamination (aw)",
  "nwalt2002_pcturb"       = "Urban land use",
  "IDW_03_05"              = "Business locations (idw, \n 3 km radium, 5 closest matches)",
  "dist_pfas_sites"        = "Distance to PFAS sites",
  "aw_percchild10"         = "Percent children under 18",
  "awairport"              = "Airports (aw)",
  "IW_AP_5K_5"             = "Airports (idw, \n 5 km radius 5 closest matches)",
  "season"                 = "Season",
  "IDW_05_05"              = "Business locations (idw, \n 5 km radius 5 closest matches)",
  "IDW_03_10"              = "Business locations (idw, \n 3 km radius 10 closest matches)",
  "awmirta"                = "Military sites (aw)",
  "aw_percpov10"           = "Percent poverty",
  "IW_AP_10K_5"            = "Airports (idw, \n 10 km radius 5 closest matches",
  "aw_mdrent10"            = "Median rent",
  "popdensity2000"         = "Population density (2000)",
  "IW_FUDS_10K_10"         = "FUDS (idw, \n 10 km radius 10 closest matches)",
  "aw_perctransp10"        = "Percent transportation occupation",
  "depth_m"                = "Well depth",
  "n_sep2010"              = "Nitrogen from septic",
  "awfuds"                 = "FUDS sites (aw)",
  "Pop_Cat_5.3.301.10.000" = "CWS Population size category",
  "aw_percwhite10"         = "Percent White",
  "aw_eduhs10"             = "Population between \n high school and college education",
  "aw_percclerical10"      = "Percent clerical occupation",
  "IW_AP_3K_5"             = "Airports (idw, \n 3 km radius 5 closest matches)",
  "IW_FIRE_5K_5"           = "Fire stations (idw, \n 5 km radius 5 closest matches)",
  "aw_percservice10"       = "Percent service occupation",
  "aw_edult810"            = "Population less than 8th grade",
  "aw_percwelfr10"         = "Percent families \n receiving public assistance",
  "IDW_05_10"              = "Business locations (idw, \n 5 km radius 10 closest matches)",
  "aw_percsales10"         = "Percent sales occupation",
  "aw_percrepair10"        = "Percent repair/construction \n occupation",
  "aw_edubach10"           = "Percent population college or more",
  "IDW_05_15"              = "Business locations (idw, \n 5 km radius 15 closest matches)",
  "IW_FIRE_1K_5"           = "Fire stations (idw, \n 1 km radius 5 closest matches)",
  "aw_edugt8lths10"        = "Percent population between \n 8th grade and high school",
  "aw_munemp10"            = "Male unemployment",
  "aw_unemp10"             = "Percent unemployed",
  "IW_FUDS_5K_5"           = "FUDS (idw, \n 5 km radius 5 closest matches)"
)

# ============================================================
# 2. Function to rename features
# ============================================================

rename_xgb_features <- function(importance_data, labels) {
  
  importance_data$Feature <- ifelse(
    importance_data$Feature %in% names(labels),
    labels[importance_data$Feature],
    importance_data$Feature
  )
  
  return(importance_data)
}


pfoa_imp_stage1_labeled <- rename_xgb_features(
  pfoa_imp_stage1, variable_labels
)

pfoa_imp_stage2_labeled <- rename_xgb_features(
  pfoa_imp_stage2, variable_labels
)

pfos_imp_stage1_labeled <- rename_xgb_features(
  pfos_imp_stage1, variable_labels
)

pfos_imp_stage2_labeled <- rename_xgb_features(
  pfos_imp_stage2, variable_labels
)

pfna_imp_stage1_labeled <- rename_xgb_features(
  pfna_imp_stage1, variable_labels
)

pfna_imp_stage2_labeled <- rename_xgb_features(
  pfna_imp_stage2, variable_labels
)


####### plots #############


# PFOA
xgb.plot.importance(
  pfoa_imp_stage1_labeled,
  top_n = 30,
  main = "PFOA - Stage 1: Detection",
  cex = 0.7
) 

xgb.plot.importance(
  pfoa_imp_stage2_labeled,
  top_n = 30,
  main = "PFOA - Stage 2: Concentration",
  cex = 0.7
)


# PFOS
xgb.plot.importance(
  pfos_imp_stage1_labeled,
  top_n = 30,
  main = "PFOS - Stage 1: Detection",
  cex = 0.7
)

xgb.plot.importance(
  pfos_imp_stage2_labeled,
  top_n = 30,
  main = "PFOS - Stage 2: Concentration",
  cex = 0.7
)


# PFNA
xgb.plot.importance(
  pfna_imp_stage1_labeled,
  top_n = 30,
  main = "PFNA - Stage 1: Detection",
  cex = 0.7
)

xgb.plot.importance(
  pfna_imp_stage2_labeled,
  top_n = 30,
  main = "PFNA - Stage 2: Concentration",
  cex = 0.7
)

####### names of variables to relabel ##########

# PFOA
head(pfoa_imp_stage1$Feature, 30)
# [1] "pm25"             "aw_percforeign10" "dtw_m"            "no2"              "rech48_mmyr"     
# [6] "aw_rcincice10"    "aw_pohse10"       "aw_mdhval10"      "pp_dns_"          "aw_percnflabor10"
# [11] "avg_clay"         "aw_mn__"          "nwalt2002_pcturb" "IDW_03_05"        "dist_pfas_sites" 
# [16] "aw_percchild10"   "awairport"        "IW_AP_5K_5"       "season"           "IDW_05_05"       
# [21] "IDW_03_10"        "awmirta"          "aw_percpov10"     "IW_AP_10K_5"      "aw_mdrent10"     
# [26] "popdensity2000"   "IW_FUDS_10K_10"   "aw_perctransp10"  "depth_m"          "n_sep2010" 


head(pfoa_imp_stage2$Feature, 30)
# [1] "nwalt2002_pcturb"       "no2"                    "pm25"                  
# [4] "IDW_03_05"              "dtw_m"                  "popdensity2000"        
# [7] "depth_m"                "IW_FUDS_10K_10"         "season"                
# [10] "IDW_03_10"              "aw_percchild10"         "aw_mdrent10"           
# [13] "awfuds"                 "IDW_05_05"              "awairport"             
# [16] "aw_percforeign10"       "aw_mn__"                "Pop_Cat_5.3.301.10.000"
# [19] "aw_percwhite10"         "aw_perctransp10"        "aw_eduhs10"            
# [22] "IW_AP_10K_5"            "n_sep2010"              "IW_AP_5K_5"            
# [25] "aw_percclerical10"      "aw_percnflabor10"       "pp_dns_"               
# [28] "IW_AP_3K_5"             "IW_FIRE_5K_5"           "aw_percservice10" 



# PFOS
head(pfos_imp_stage1$Feature, 30)
# [1] "pm25"             "dtw_m"            "no2"              "nwalt2002_pcturb" "avg_clay"        
# [6] "IDW_03_05"        "aw_percforeign10" "popdensity2000"   "aw_pohse10"       "IDW_03_10"       
# [11] "IW_AP_10K_5"      "IDW_05_05"        "dist_pfas_sites"  "awmirta"          "n_sep2010"       
# [16] "aw_percchild10"   "pp_dns_"          "aw_edult810"      "aw_percwelfr10"   "aw_rcincice10"   
# [21] "depth_m"          "IDW_05_10"        "aw_percsales10"   "awfuds"           "season"          
# [26] "aw_percpov10"     "rech48_mmyr"      "aw_percrepair10"  "aw_percwhite10"   "IW_FUDS_10K_10" 

head(pfos_imp_stage2$Feature, 30)
# [1] "popdensity2000"   "no2"              "pp_dns_"          "IDW_03_05"        "pm25"            
# [6] "aw_mn__"          "aw_pohse10"       "aw_percchild10"   "awmirta"          "season"          
# [11] "IDW_03_10"        "IW_FUDS_10K_10"   "aw_percsales10"   "aw_mdhval10"      "aw_mdrent10"     
# [16] "IW_AP_10K_5"      "aw_edubach10"     "aw_edult810"      "nwalt2002_pcturb" "n_sep2010"       
# [21] "aw_percwhite10"   "IW_FIRE_5K_5"     "aw_rcincice10"    "aw_percnflabor10" "IDW_05_15"       
# [26] "IW_FIRE_1K_5"     "aw_percforeign10" "rech48_mmyr"      "aw_edugt8lths10"  "IDW_05_05"  


# PFNA
head(pfna_imp_stage1$Feature, 30)
# [1] "IW_AP_10K_5"       "IDW_03_05"         "no2"               "pm25"             
# [5] "aw_percnflabor10"  "aw_mdhval10"       "season"            "awmirta"          
# [9] "aw_percclerical10" "IW_AP_5K_5"        "IW_FUDS_10K_10"    "rech48_mmyr"      
# [13] "dist_pfas_sites"   "nwalt2002_pcturb"  "awairport"         "aw_munemp10"      
# [17] "aw_unemp10"        "IDW_05_05"         "aw_percforeign10"  "aw_percchild10"   
# [21] "depth_m"           "awfuds"            "IDW_03_10"         "aw_edult810"      
# [25] "IW_FUDS_5K_5"      "aw_mn__"           "aw_percsales10"    "aw_percwelfr10"   
# [29] "avg_clay"          "pp_dns_"  



head(pfna_imp_stage2$Feature, 30)
# [1] "pm25"                   "awairport"              "aw_percforeign10"      
# [4] "aw_mdhval10"            "IDW_03_05"              "awmirta"               
# [7] "season"                 "Pop_Cat_5.3.301.10.000" "aw_edult810"           
# [10] "IW_FIRE_5K_5"           "aw_eduhs10"             "aw_mn__"               
# [13] "no2"                    "IDW_05_05"              "aw_percchild10"        
# [16] "IDW_05_15"              "aw_mdrent10"            "aw_percclerical10"     
# [19] "aw_percwhite10"         "IDW_03_10"              "pp_dns_"               
# [22] "depth_m"                "aw_percwelfr10"         "aw_pohse10"            
# [25] "IW_FIRE_1K_5"           "aw_munemp10"            "IW_FUDS_10K_10"        
# [28] "aw_percnflabor10"       "aw_edugt8lths10"        "rech48_mmyr"  

# all variables 
head(pfoa_imp_stage1$Feature, 100)
# [1] "pm25"                     "aw_percforeign10"         "dtw_m"                   
# [4] "no2"                      "rech48_mmyr"              "aw_rcincice10"           
# [7] "aw_pohse10"               "aw_mdhval10"              "pp_dns_"                 
# [10] "aw_percnflabor10"         "avg_clay"                 "aw_mn__"                 
# [13] "nwalt2002_pcturb"         "IDW_03_05"                "dist_pfas_sites"         
# [16] "aw_percchild10"           "awairport"                "IW_AP_5K_5"              
# [19] "season"                   "IDW_05_05"                "IDW_03_10"               
# [22] "awmirta"                  "aw_percpov10"             "IW_AP_10K_5"             
# [25] "aw_mdrent10"              "popdensity2000"           "IW_FUDS_10K_10"          
# [28] "aw_perctransp10"          "depth_m"                  "n_sep2010"               
# [31] "aw_eduhs10"               "awfuds"                   "aw_percwelfr10"          
# [34] "aw_munemp10"              "aw_incice10"              "aw_percclerical10"       
# [37] "aw_edugt8lths10"          "aw_pfwchbf10"             "IDW_05_10"               
# [40] "aw_edult810"              "aw_rice10"                "aw_percsales10"          
# [43] "IW_FIRE_3K_5"             "aw_percblack10"           "IW_FIRE_5K_5"            
# [46] "IW_FIRE_5K_10"            "IDW_05_15"                "aw_percold10"            
# [49] "aw_percrepair10"          "IW_FIRE_5K_15"            "aw_mdinc10"              
# [52] "aw_unemp10"               "aw_edubach10"             "aw_percwhite10"          
# [55] "aw_percservice10"         "IW_AP_3K_5"               "IW_FIRE_1K_5"            
# [58] "IW_MIRTA_10K_5"           "IW_FUDS_3K_5"             "IW_FUDS_5K_5"            
# [61] "aw_percfamint10"          "IW_FIRE_3K_10"            "IW_FUDS_1K_5"            
# [64] "Pop_Cat_5.3.301.10.000"   "Pop_Cat_5.10.001.100.000" "IW_MIRTA_5K_5"           
# [67] "Pop_Cat_5...500" 

# training/test Ns

# PFOA
n_train_pfoa <- length(res_pfoa$idx_train)
n_total_pfoa <- length(res_pfoa$stage1$preds_train) +
  length(res_pfoa$stage1$preds_test)

n_test_pfoa <- n_total_pfoa - n_train_pfoa


n_train_pfoa
n_total_pfoa
n_test_pfoa


# PFOA
pfoa_train_n <- length(res_pfoa$idx_train)
pfoa_total_n <- nrow(nj_prep_xg_pfoa)
pfoa_test_n  <- pfoa_total_n - pfoa_train_n

# PFOS
pfos_train_n <- length(res_pfos$idx_train)
pfos_total_n <- nrow(nj_prep_xg_pfos)
pfos_test_n  <- pfos_total_n - pfos_train_n

# PFNA
pfna_train_n <- length(res_pfna$idx_train)
pfna_total_n <- nrow(nj_prep_xg_pfna)
pfna_test_n  <- pfna_total_n - pfna_train_n

split_counts <- data.frame(
  PFAS = c("PFOA", "PFOS", "PFNA"),
  
  Training_N = c(
    length(res_pfoa$idx_train),
    length(res_pfos$idx_train),
    length(res_pfna$idx_train)
  ),
  
  Test_N = c(
    nrow(nj_prep_xg_pfoa) - length(res_pfoa$idx_train),
    nrow(nj_prep_xg_pfos) - length(res_pfos$idx_train),
    nrow(nj_prep_xg_pfna) - length(res_pfna$idx_train)
  )
)

split_counts
# PFAS Training_N Test_N
# 1 PFOA       6160   2597
# 2 PFOS       6072   2559
# 3 PFNA       6081   2563


