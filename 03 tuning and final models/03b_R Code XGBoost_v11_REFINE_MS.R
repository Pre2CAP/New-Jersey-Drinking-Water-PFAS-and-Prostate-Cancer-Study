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

load("/home/ms3283/R/XGBoost/data/pfas_xg.Rda")
nj_cws_pfas_dist <- st_read("/home/ms3283/R/XGBoost/data/nj_cws_pfas_dist.shp")
load("/home/ms3283/R/XGBoost/data/cws_aw_means.Rda")
load("/home/ms3283/R/XGBoost/data/Pollution/njpws_pm25_joined.Rda")
load("/home/ms3283/R/XGBoost/data/Pollution/njpws_no2_joined.Rda")
load("/home/ms3283/R/XGBoost/data/tok_cws_joined.Rda")
nj_mirta <- st_read("/home/ms3283/R/XGBoost/data/nj_cws_mirta_dist.shp")
nj_fuds <- st_read("/home/ms3283/R/XGBoost/data/nj_cws_fuds_dist.shp")
nj_airport <- st_read("/home/ms3283/R/XGBoost/data/nj_cws_airp_dist.shp")
load("/home/ms3283/R/XGBoost/data/nj_cws_idw.Rda")

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
  filter(pfas_name == "PFOA") ## 8800

nj_prep_xg_pfos <- nj_prep_xg %>%
  filter(pfas_name == "PFOS") ## 8674

nj_prep_xg_pfna <- nj_prep_xg %>%
  filter(pfas_name == "PFNA") ## 8687

summary(nj_prep_xg_pfoa)

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

# best hyperparameters from previous models 

tune_grid_bin_pfoa <- expand.grid(
  nrounds = 500,
  eta = 0.08,
  gamma = 0.05,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 1,
  lambda = 1,
  alpha = 1,
  max_depth = 4
)

tune_grid_reg_pfoa <- expand.grid(
  nrounds = 250,
  eta = 0.04,
  max_depth = 6,
  gamma = 0.05,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  lambda = 1,
  alpha = 0
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

res_pfoa$stage1$best_hyperparameters
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds  cv_auc cv_logloss
#     500 0.08  0.05              0.7                1         1      1     1         4          344 0.91016  0.3776461

# start tuniing 

# Base (previous best)
base_bin <- list(
  nrounds = 500,
  eta = 0.08,
  gamma = 0.05,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 1,
  lambda = 1,
  alpha = 1,
  max_depth = 4
)

# Step 1 grid: vary core tree + learning parameters
tune_grid_bin_step1 <- expand.grid(
  nrounds = c(400, 500, 600),
  eta = c(0.06, 0.08, 0.10),
  gamma = base_bin$gamma,
  colsample_bytree = base_bin$colsample_bytree,
  min_child_weight = c(1, 2, 4),
  subsample = base_bin$subsample,
  lambda = base_bin$lambda,
  alpha = base_bin$alpha,
  max_depth = c(3, 4, 5)
)

res_bin_step1 <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_step1,
  tune_grid_reg = tune_grid_reg_pfoa,  # leave as-is for now
  nfold = 5,
  early_stopping_rounds = 25
)

res_bin_step1$stage1$best_hyperparameters
#nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
#     400 0.08  0.05              0.7                2         1      1     1         4          356 0.9149322  0.3690314


# step 1 is an improvement to original model and all parameters are within boundaries so no refinement is needed


base_bin_pfoa_step2 <- list(
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



tune_grid_bin_pfoa_step2 <- expand.grid(
  nrounds = base_bin_pfoa_step2$nrounds,
  eta = base_bin_pfoa_step2$eta,
  gamma = base_bin_pfoa_step2$gamma,
  colsample_bytree = c(0.6, 0.7, 0.8),
  min_child_weight = base_bin_pfoa_step2$min_child_weight,
  subsample = c(0.8, 0.9, 1.0),
  lambda = base_bin_pfoa_step2$lambda,
  alpha = base_bin_pfoa_step2$alpha,
  max_depth = base_bin_pfoa_step2$max_depth
)


res_bin_step2 <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfoa_step2,
  tune_grid_reg = tune_grid_reg_pfoa, 
  nfold = 5,
  early_stopping_rounds = 25
)

res_bin_step2$stage1$best_hyperparameters
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
#     400 0.08  0.05              0.7                2         1      1     1         4          348 0.9127234  0.3726695

# performed worse, sampling did not help, step 1 model was already optimal do not use step 2 result, keep step 1 model and move on to regularization

base_bin_pfoa_step3 <- list(
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

tune_grid_bin_pfoa_step3 <- expand.grid(
  nrounds = base_bin_pfoa_step3$nrounds,
  eta = base_bin_pfoa_step3$eta,
  gamma = c(0.01, 0.05, 0.1),
  colsample_bytree = base_bin_pfoa_step3$colsample_bytree,
  min_child_weight = base_bin_pfoa_step3$min_child_weight,
  subsample = base_bin_pfoa_step3$subsample,
  lambda = c(0.5, 1, 2),
  alpha = c(0.5, 1, 2),
  max_depth = base_bin_pfoa_step3$max_depth
)

res_bin_step3 <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfoa_step3,
  tune_grid_reg = tune_grid_reg_pfoa, 
  nfold = 5,
  early_stopping_rounds = 25
)

res_bin_step3$stage1$best_hyperparameters
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
#     400 0.08  0.05              0.7                2         1      1   0.5         4          400 0.9139626  0.3705556


# this does not improve from step 1, locking in step 1 as final model 

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
  nrounds = 250,
  eta = 0.04,
  max_depth = 6,
  gamma = 0.05,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  lambda = 1,
  alpha = 0
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

res_pfoa$stage1$best_hyperparameters
#nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
#     400 0.08  0.05              0.7                2         1      1     1         4          350 0.9114539  0.3759071




## moving on to PFOA reg model 

# step 1 tune max depth, min_child_weight, eta, and nrounds

#first try with original parameters

tune_grid_reg_pfoa <- expand.grid(
  nrounds = 250,
  eta = 0.04,
  max_depth = 6,
  gamma = 0.05,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  lambda = 1,
  alpha = 0
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

res_pfoa$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#     250 0.04         6  0.05              0.8                1       0.9      1     0          111 0.4344468



## now tuning ##

base_reg_step1 <- list(
  nrounds = 250,
  eta = 0.04,
  max_depth = 6,
  gamma = 0.05,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  lambda = 1,
  alpha = 0
)

tune_grid_reg_step1 <- expand.grid(
  nrounds = c(200, 250, 300),
  eta = c(0.03, 0.04, 0.05),
  max_depth = c(5, 6, 7),
  gamma = base_reg_step1$gamma,
  colsample_bytree = base_reg_step1$colsample_bytree,
  min_child_weight = c(1, 2, 4),
  subsample = base_reg_step1$subsample,
  lambda = base_reg_step1$lambda,
  alpha = base_reg_step1$alpha
)

res_reg_step1 <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfoa,
  tune_grid_reg = tune_grid_reg_step1,
  nfold = 5,
  early_stopping_rounds = 25
)

res_reg_step1$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#     250 0.03         7  0.05              0.8                1       0.9      1     0          151 0.4262557


# tuning improved step one but parameters close to boundaries so need to refine further

tune_grid_reg_step1_refine <- expand.grid(
  nrounds = c(200, 300, 400),
  eta = c(0.02, 0.025, 0.03, 0.035),
  max_depth = c(6, 7, 8),
  gamma = 0.05,
  colsample_bytree = 0.8,
  min_child_weight = c(1, 1.5, 2),
  subsample = 0.9,
  lambda = 1,
  alpha = 0
)

res_reg_step1_refine <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfoa,
  tune_grid_reg = tune_grid_reg_step1_refine,
  nfold = 5,
  early_stopping_rounds = 25
)

res_reg_step1_refine$stage2$best_hyperparameters
#nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#     300 0.03         7  0.05              0.8                1       0.9      1     0          151 0.4262557

# the exact same parameters were selected, step 1 is complete, moving onto step 2 sampling

base_reg_pfoa_step2 <- list(
  nrounds = 300,   # cap; early stopping will handle
  eta = 0.03,
  max_depth = 7,
  gamma = 0.05,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  lambda = 1,
  alpha = 0
)

tune_grid_reg_pfoa_step2 <- expand.grid(
  nrounds = base_reg_pfoa_step2$nrounds,
  eta = base_reg_pfoa_step2$eta,
  max_depth = base_reg_pfoa_step2$max_depth,
  gamma = base_reg_pfoa_step2$gamma,
  colsample_bytree = c(0.7, 0.8, 0.9),
  min_child_weight = base_reg_pfoa_step2$min_child_weight,
  subsample = c(0.8, 0.9, 1.0),
  lambda = base_reg_pfoa_step2$lambda,
  alpha = base_reg_pfoa_step2$alpha
)

res_reg_step2 <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfoa,
  tune_grid_reg = tune_grid_reg_pfoa_step2,
  nfold = 5,
  early_stopping_rounds = 25
)

res_reg_step2$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#     300 0.03         7  0.05              0.7                1         1      1     0          146 0.4298426

# this does not perform better, lock step 1

base_reg_pfoa_step3 <- list(
  nrounds = 300,
  eta = 0.03,
  max_depth = 7,
  gamma = 0.05,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  lambda = 1,
  alpha = 0
)

tune_grid_reg_pfoa_step3 <- expand.grid(
  nrounds = base_reg_pfoa_step3$nrounds,
  eta = base_reg_pfoa_step3$eta,
  max_depth = base_reg_pfoa_step3$max_depth,
  gamma = c(0.01, 0.05, 0.1),
  colsample_bytree = base_reg_pfoa_step3$colsample_bytree,
  min_child_weight = base_reg_pfoa_step3$min_child_weight,
  subsample = base_reg_pfoa_step3$subsample,
  lambda = c(0.5, 1, 2),
  alpha = c(0, 0.25, 0.5, 1)
)

res_reg_step3 <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfoa,
  tune_grid_reg = tune_grid_reg_pfoa_step3,
  nfold = 5,
  early_stopping_rounds = 25
)

res_reg_step3$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#     300 0.03         7  0.05              0.8                1       0.9    0.5   0.5          186 0.4253728

# this actually does perform better, lock step 3 as final PFOA regression model

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
#print(res_pfoa$stage1$auc_plot)

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
load("/home/ms3283/R/XGBoost/models/final models/res_pfoa_4.30.26.Rda")



#####################################REDO

#######################################################################
# PFOS
#######################################################################

## start by tuning max_deptth, min_child_weight, eta, rounds
## Identify best parameters
## then move to subsample, colsample_bytree
## then move to gamma, lambda, alpha

###############################################
# PFOS 
###############################################

#starting with best model from before

tune_grid_bin_pfos <- expand.grid(
  nrounds = 800,
  eta = 0.05,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 0.6,
  min_child_weight = 1,
  subsample = 1,
  lambda = 1,
  alpha = 0
)

tune_grid_reg_pfos <- expand.grid(
  nrounds = 800,
  eta = 0.05,
  max_depth = 6,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1,
  lambda = 0,
  alpha = 0
)

res_pfos <- run_two_phase_xgb(
  df = nj_prep_xg_pfos, ## change
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos,
  nfold = 5, 
  early_stopping_rounds = 25 ##25
)

res_pfos$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#     800 0.05         6     0              0.6                1         1      1     0          336 0.9112761  0.3537492



# first step is to tune max_depth, min_child_weight, eta, and nrounds

# PFOS binary base from previous best model
base_bin_pfos_step1 <- list(
  nrounds = 800,
  eta = 0.05,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 0.6,
  min_child_weight = 1,
  subsample = 1,
  lambda = 1,
  alpha = 0
)

# Step 1: tune max_depth, min_child_weight, eta, nrounds
tune_grid_bin_pfos_step1 <- expand.grid(
  nrounds = c(600, 800, 1000),
  eta = c(0.04, 0.05, 0.06),
  max_depth = c(5, 6, 7),
  gamma = base_bin_pfos_step1$gamma,
  colsample_bytree = base_bin_pfos_step1$colsample_bytree,
  min_child_weight = c(1, 2, 4),
  subsample = base_bin_pfos_step1$subsample,
  lambda = base_bin_pfos_step1$lambda,
  alpha = base_bin_pfos_step1$alpha
)

res_pfos_step1 <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos_step1,
  tune_grid_reg = tune_grid_reg_pfos,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfos_step1$stage1$best_hyperparameters
#nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.05         7     0              0.6                2         1      1     0          250 0.9113957    0.35239



# tuning provides a small but meaningful improvement, but need to expand the boundaries 

tune_grid_bin_pfos_step1_refine <- expand.grid(
  nrounds = c(900, 1000, 1100),
  eta = c(0.045, 0.05, 0.055),
  max_depth = c(6, 7, 8),
  gamma = 0,
  colsample_bytree = 0.6,
  min_child_weight = c(1.5, 2, 3),
  subsample = 1,
  lambda = 1,
  alpha = 0
)

res_pfos_step1_refine <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos_step1_refine,
  tune_grid_reg = tune_grid_reg_pfos,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfos_step1_refine$stage1$best_hyperparameters
#nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1100 0.05         8     0              0.6                2         1      1     0          167 0.9113181  0.3524434

# this refinement did not help, do not use refine model, use step 1 model

base_bin_pfos_step2 <- list(
  nrounds = 1000,   # cap is fine; early stopping handles it
  eta = 0.05,
  max_depth = 7,
  gamma = 0,
  colsample_bytree = 0.6,
  min_child_weight = 2,
  subsample = 1,
  lambda = 1,
  alpha = 0
)

tune_grid_bin_pfos_step2 <- expand.grid(
  nrounds = base_bin_pfos_step2$nrounds,
  eta = base_bin_pfos_step2$eta,
  max_depth = base_bin_pfos_step2$max_depth,
  gamma = base_bin_pfos_step2$gamma,
  colsample_bytree = c(0.5, 0.6, 0.7),
  min_child_weight = base_bin_pfos_step2$min_child_weight,
  subsample = c(0.8, 0.9, 1.0),
  lambda = base_bin_pfos_step2$lambda,
  alpha = base_bin_pfos_step2$alpha
)

res_pfos_step2 <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos_step2,
  tune_grid_reg = tune_grid_reg_pfos,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfos_step2$stage1$best_hyperparameters
#nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.05         7     0              0.5                2         1      1     0          204 0.9101211  0.3539765

# no improvement, do not use step 2 result, keep step 1 model, move on to step 3 regularization

base_bin_pfos_step3 <- list(
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

tune_grid_bin_pfos_step3 <- expand.grid(
  nrounds = base_bin_pfos_step3$nrounds,
  eta = base_bin_pfos_step3$eta,
  max_depth = base_bin_pfos_step3$max_depth,
  gamma = c(0, 0.01, 0.05, 0.1),
  colsample_bytree = base_bin_pfos_step3$colsample_bytree,
  min_child_weight = base_bin_pfos_step3$min_child_weight,
  subsample = base_bin_pfos_step3$subsample,
  lambda = c(0.5, 1, 2),
  alpha = c(0, 0.25, 0.5, 1)
)

res_pfos_step3 <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos_step3,
  tune_grid_reg = tune_grid_reg_pfos,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfos_step3$stage1$best_hyperparameters
#nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_auc cv_logloss
#   1000 0.05         7     0              0.6                2         1    0.5     0          212 0.911196  0.3529128


# auc and cv log loss got worse, reject the new model and do not use step 3 result. Lock in step 1 mdoel as final

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
  nrounds = 800,
  eta = 0.05,
  max_depth = 6,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1,
  lambda = 0,
  alpha = 0
)

res_pfos <- run_two_phase_xgb(
  df = nj_prep_xg_pfos, ## change
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos,
  nfold = 5, 
  early_stopping_rounds = 25 ##25
)

res_pfos$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.05         7     0              0.6                2         1      1     0          213 0.9113486  0.3525931



###################### moving on to PFOS reg

tune_grid_bin_pfos <- expand.grid(
  nrounds = 500,
  eta = 0.07,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 0.6,
  min_child_weight = 1.5,
  subsample = 1,
  lambda = 1,
  alpha = 0
)

tune_grid_reg_pfos <- expand.grid(
  nrounds = 800,
  eta = 0.05,
  max_depth = 6,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1,
  lambda = 0,
  alpha = 0
)

res_pfos <- run_two_phase_xgb(
  df = nj_prep_xg_pfos, ## change
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos,
  nfold = 5,
  early_stopping_rounds = 25 ##25
)

res_pfos$stage2$best_hyperparameters
#nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#     800 0.05         6   0.1                1                1         1      0     0           76 0.4559054


# tuning

base_reg_pfos_step1 <- list(
  nrounds = 800,
  eta = 0.05,
  max_depth = 6,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1,
  lambda = 0,
  alpha = 0
)


tune_grid_reg_pfos_step1 <- expand.grid(
  nrounds = c(600, 800, 1000),
  eta = c(0.04, 0.05, 0.06),
  max_depth = c(5, 6, 7),
  gamma = base_reg_pfos_step1$gamma,
  colsample_bytree = base_reg_pfos_step1$colsample_bytree,
  min_child_weight = c(1, 2, 4),
  subsample = base_reg_pfos_step1$subsample,
  lambda = base_reg_pfos_step1$lambda,
  alpha = base_reg_pfos_step1$alpha
)

res_pfos_reg_step1 <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos_step1,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfos_reg_step1$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#    1000 0.04         6   0.1                1                1         1      0     0           80 0.4517598


# tuning improved the mode, but doing a small refinement 
tune_grid_reg_pfos_step1_refine <- expand.grid(
  nrounds = c(800, 1000),
  eta = c(0.03, 0.035, 0.04, 0.045),
  max_depth = c(5, 6, 7),
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = c(1, 2),
  subsample = 1,
  lambda = 0,
  alpha = 0
)

res_pfos_reg_step1_refine <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos_step1_refine,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfos_reg_step1_refine$stage2$best_hyperparameters
# nrounds   eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#    1000 0.035         6   0.1                1                1         1      0     0           99 0.4519594

# refinement step performed slightly worse, locking step 1 and moving onto step 2 samplling

base_reg_pfos_step2 <- list(
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


tune_grid_reg_pfos_step2 <- expand.grid(
  nrounds = base_reg_pfos_step2$nrounds,
  eta = base_reg_pfos_step2$eta,
  max_depth = base_reg_pfos_step2$max_depth,
  gamma = base_reg_pfos_step2$gamma,
  colsample_bytree = c(0.8, 0.9, 1.0),
  min_child_weight = base_reg_pfos_step2$min_child_weight,
  subsample = c(0.8, 0.9, 1.0),
  lambda = base_reg_pfos_step2$lambda,
  alpha = base_reg_pfos_step2$alpha
)


res_pfos_reg_step2 <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos_step2,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfos_reg_step2$stage2$best_hyperparameters
#nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#    1000 0.04         6   0.1              0.8                1       0.8      0     0           94 0.4544659

# this performs worse, reject dont use and keep step 1, moving onto regularization

base_reg_pfos_step3 <- list(
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

tune_grid_reg_pfos_step3 <- expand.grid(
  nrounds = base_reg_pfos_step3$nrounds,
  eta = base_reg_pfos_step3$eta,
  max_depth = base_reg_pfos_step3$max_depth,
  gamma = c(0.05, 0.1, 0.2),
  colsample_bytree = base_reg_pfos_step3$colsample_bytree,
  min_child_weight = base_reg_pfos_step3$min_child_weight,
  subsample = base_reg_pfos_step3$subsample,
  lambda = c(0, 0.5, 1, 2),
  alpha = c(0, 0.25, 0.5, 1)
)

res_pfos_reg_step3 <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos_step3,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfos_reg_step3$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds  cv_rmse
#    1000 0.04         6   0.1                1                1         1      1   0.5          103 0.454133

# this performs worse, don't use this model, use step 1 model as final

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
load("/home/ms3283/R/XGBoost/models/final models/res_pfos_4.30.26.Rda")


#######################################################################
# PFNA
#######################################################################

tune_grid_bin_pfna <- expand.grid(
  nrounds = 800,
  eta = 0.05,
  max_depth = 4,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 5,
  subsample = 1,
  lambda = 1,
  alpha = 0
)



tune_grid_reg_pfna <- expand.grid(
  nrounds = 800,
  eta = 0.01,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 5,
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)

res_pfna <- run_two_phase_xgb(
  df = nj_prep_xg_pfna, ## change
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5, 
  early_stopping_rounds = 25 ##25
)

res_pfna$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#     800 0.05         4   0.1                1                5         1      1     0          290 0.9241937  0.1567988




# starting tuning


base_bin_pfna_step1 <- list(
  nrounds = 800,
  eta = 0.05,
  max_depth = 4,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 5,
  subsample = 1,
  lambda = 1,
  alpha = 0
)

tune_grid_bin_pfna_step1 <- expand.grid(
  nrounds = c(600, 800, 1000),
  eta = c(0.04, 0.05, 0.06),
  max_depth = c(3, 4, 5),
  gamma = base_bin_pfna_step1$gamma,
  colsample_bytree = base_bin_pfna_step1$colsample_bytree,
  min_child_weight = c(3, 5, 7),
  subsample = base_bin_pfna_step1$subsample,
  lambda = base_bin_pfna_step1$lambda,
  alpha = base_bin_pfna_step1$alpha
)

res_pfna_step1 <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna_step1,
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_step1$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.04         4   0.1                1                3         1      1     0          344 0.9282637  0.1533818


# step 1 had improvement, but parameters close to boundaries, refining further

tune_grid_bin_pfna_step1_refine <- expand.grid(
  nrounds = c(800, 1000),
  eta = c(0.035, 0.04, 0.045),
  max_depth = c(3, 4, 5),
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = c(2, 3, 4),
  subsample = 1,
  lambda = 1,
  alpha = 0
)

res_pfna_step1_refine <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna_step1_refine,
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_step1_refine$stage1$best_hyperparameters
# nrounds   eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.045         4   0.1                1                3         1      1     0          329 0.9276812  0.1536493


# performs worse, using original step 1 and moving onto step 2, sampling

base_bin_pfna_step2 <- list(
  nrounds = 1000,
  eta = 0.04,
  max_depth = 4,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 3,
  subsample = 1,
  lambda = 1,
  alpha = 0
)

tune_grid_bin_pfna_step2 <- expand.grid(
  nrounds = base_bin_pfna_step2$nrounds,
  eta = base_bin_pfna_step2$eta,
  max_depth = base_bin_pfna_step2$max_depth,
  gamma = base_bin_pfna_step2$gamma,
  colsample_bytree = c(0.8, 0.9, 1.0),
  min_child_weight = base_bin_pfna_step2$min_child_weight,
  subsample = c(0.8, 0.9, 1.0),
  lambda = base_bin_pfna_step2$lambda,
  alpha = base_bin_pfna_step2$alpha
)

res_pfna_step2 <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna_step2,
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_step2$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.04         4   0.1              0.8                3       0.9      1     0          299 0.9289513  0.1533335

# this does show improvement, sampling does help slightly, but worth doing a step 2 refinement
tune_grid_bin_pfna_step2_refine <- expand.grid(
  nrounds = base_bin_pfna_step2$nrounds,
  eta = base_bin_pfna_step2$eta,
  max_depth = base_bin_pfna_step2$max_depth,
  gamma = base_bin_pfna_step2$gamma,
  colsample_bytree = c(0.75, 0.8, 0.85),
  min_child_weight = base_bin_pfna_step2$min_child_weight,
  subsample = c(0.85, 0.9, 0.95),
  lambda = base_bin_pfna_step2$lambda,
  alpha = base_bin_pfna_step2$alpha
)

res_pfna_step2_refine <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna_step2_refine,
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_step2_refine$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.04         4   0.1             0.75                3       0.9      1     0          300 0.9288022  0.1533761

# this performed worse, so using step 2 result and moving onto step 3 regularizatioin

base_bin_pfna_step3 <- list(
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

tune_grid_bin_pfna_step3 <- expand.grid(
  nrounds = base_bin_pfna_step3$nrounds,
  eta = base_bin_pfna_step3$eta,
  max_depth = base_bin_pfna_step3$max_depth,
  gamma = c(0.05, 0.1, 0.2),
  colsample_bytree = base_bin_pfna_step3$colsample_bytree,
  min_child_weight = base_bin_pfna_step3$min_child_weight,
  subsample = base_bin_pfna_step3$subsample,
  lambda = c(0.5, 1, 2),
  alpha = c(0, 0.25, 0.5, 1)
)

res_pfna_step3 <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna_step3,
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_step3$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.04         4  0.05              0.8                3       0.9      1     0          299 0.9288542  0.1533331

# performed worse or about the same, did not improve model

# final bin model

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
  nrounds = 800,
  eta = 0.01,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 5,
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)

res_pfna <- run_two_phase_xgb(
  df = nj_prep_xg_pfna, ## change
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5, 
  early_stopping_rounds = 25 ##25
)

res_pfna$stage1$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds    cv_auc cv_logloss
#    1000 0.04         4   0.1              0.8                3       0.9      1     0          279 0.9238612  0.1572805


# moving onto PFNA regression

tune_grid_bin_pfna <- expand.grid(
  nrounds = 600,
  eta = 0.05,
  max_depth = 4,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 7,
  subsample = 1,
  lambda = 1,
  alpha = 0
)

tune_grid_reg_pfna <- expand.grid(
  nrounds = 800,
  eta = 0.01,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 5,
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)

res_pfna <- run_two_phase_xgb(
  df = nj_prep_xg_pfna, ## change
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,
  tune_grid_reg = tune_grid_reg_pfna,
  nfold = 5, 
  early_stopping_rounds = 25 ##25
)

res_pfna$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds  cv_rmse
#     800 0.01         6     0                1                5       0.8      0     0          333 0.581896


## tuning

base_reg_pfna_step1 <- list(
  nrounds = 800,
  eta = 0.01,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 5,
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)

tune_grid_reg_pfna_step1 <- expand.grid(
  nrounds = c(600, 800, 1000),
  eta = c(0.008, 0.01, 0.015, 0.02),
  max_depth = c(5, 6, 7),
  gamma = base_reg_pfna_step1$gamma,
  colsample_bytree = base_reg_pfna_step1$colsample_bytree,
  min_child_weight = c(3, 5, 7),
  subsample = base_reg_pfna_step1$subsample,
  lambda = base_reg_pfna_step1$lambda,
  alpha = base_reg_pfna_step1$alpha
)

res_pfna_reg_step1 <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,   # fixed binary
  tune_grid_reg = tune_grid_reg_pfna_step1,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_reg_step1$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#     800 0.02         7     0                1                3       0.8      0     0          149 0.5429467

# tuned model performs way better, parameters at boundaries so further refining is needed

tune_grid_reg_pfna_step1_refine <- expand.grid(
  nrounds = c(600, 800, 1000),
  eta = c(0.015, 0.02, 0.025, 0.03),
  max_depth = c(6, 7, 8),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(2, 3, 4),
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)

res_pfna_reg_step1_refine <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,   # fixed binary
  tune_grid_reg = tune_grid_reg_pfna_step1_refine,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_reg_step1_refine$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#    1000 0.03         6     0                1                4       0.8      0     0          122 0.5416239


# improvement, but eta hit boundary again, further refinement for eta

tune_grid_reg_pfna_step1_final <- expand.grid(
  nrounds = 1000,
  eta = c(0.025, 0.03, 0.035, 0.04),
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 4,
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)


res_pfna_reg_step1_final <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,   # fixed binary
  tune_grid_reg = tune_grid_reg_pfna_step1_final,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_reg_step1_final$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#    1000 0.03         6     0                1                4       0.8      0     0          128 0.5722415


tune_grid_reg_check <- expand.grid(
  nrounds = 1000,
  eta = 0.03,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 4,
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)


res_pfna_reg_check <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,   # fixed binary
  tune_grid_reg = tune_grid_reg_check,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_reg_check$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds  cv_rmse
#    1000 0.03         6     0                1                4       0.8      0     0          105 0.577967

# the increase in rmse is likely from cross validation instability, trust the grid search result and use step 1 result


base_reg_pfna_step2 <- list(
  nrounds = 1000,
  eta = 0.03,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 4,
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)

tune_grid_reg_pfna_step2 <- expand.grid(
  nrounds = base_reg_pfna_step2$nrounds,
  eta = base_reg_pfna_step2$eta,
  max_depth = base_reg_pfna_step2$max_depth,
  gamma = base_reg_pfna_step2$gamma,
  colsample_bytree = c(0.8, 0.9, 1.0),
  min_child_weight = base_reg_pfna_step2$min_child_weight,
  subsample = c(0.7, 0.8, 0.9),
  lambda = base_reg_pfna_step2$lambda,
  alpha = base_reg_pfna_step2$alpha
)

res_pfna_reg_step2 <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,   # fixed binary
  tune_grid_reg = tune_grid_reg_pfna_step2,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_reg_step2$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds  cv_rmse
#   1000 0.03         6     0              0.9                4       0.7      0     0          126 0.564561


# this performs worse, move onto step 3 regularization


base_reg_pfna_step3 <- list(
  nrounds = 1000,
  eta = 0.03,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 4,
  subsample = 0.8,
  lambda = 0,
  alpha = 0
)

tune_grid_reg_pfna_step3 <- expand.grid(
  nrounds = base_reg_pfna_step3$nrounds,
  eta = base_reg_pfna_step3$eta,
  max_depth = base_reg_pfna_step3$max_depth,
  gamma = c(0, 0.05, 0.1),
  colsample_bytree = base_reg_pfna_step3$colsample_bytree,
  min_child_weight = base_reg_pfna_step3$min_child_weight,
  subsample = base_reg_pfna_step3$subsample,
  lambda = c(0, 0.5, 1, 2),
  alpha = c(0, 0.25, 0.5, 1)
)

res_pfna_reg_step3 <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  cont_outcome = "lnres",
  bin_var = "pfas_bin_num",
  det_value = 1,
  seed = 123,
  train_prop = 0.7,
  tune_grid_bin = tune_grid_bin_pfna,   # fixed binary
  tune_grid_reg = tune_grid_reg_pfna_step3,
  nfold = 5,
  early_stopping_rounds = 25
)

res_pfna_reg_step3$stage2$best_hyperparameters
# nrounds  eta max_depth gamma colsample_bytree min_child_weight subsample lambda alpha best_nrounds   cv_rmse
#    1000 0.03         6  0.05                1                4       0.8      2   0.5          212 0.5421506

# performs very similar to step 1, so the unconstrained model from step 1 was slightly too flexible, regularization improves generalization not always CV mean
# even though RMSE is essentially equal, gained more stable model, less overfitting risk, and better interpretability 

# final PFNA model


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
load("/home/ms3283/R/XGBoost/models/final models/res_pfna_4.30.26.Rda")

