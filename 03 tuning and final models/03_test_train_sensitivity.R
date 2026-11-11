# Two-stage XGBoost with temporal or population-density holdout
# Time: train 2006–2017; test 2018–2025.
# Geographic: train > 1000 people/mi^2; test < 1000 people/mi^2.
# Rows with missing split values or exactly 1000 density are excluded.
# Set density_var to the exact column name in your data.


library(dplyr)
library(caret)
library(xgboost)
library(pROC)
library(ggplot2)

############### load data ########################################
load("~/R/registry/nj_prep_xg_pfoa.Rda")

load("~/R/registry/nj_prep_xg_pfos.Rda")

load("~/R/registry/nj_prep_xg_pfna.Rda")

# checking population density variable

nj_prep_xg_pfoa <- nj_prep_xg_pfoa |>
  mutate(pop_density = case_when(popdensity2000 < 1000 ~ "<1000",
                                 popdensity2000 >= 1000 ~ "1000+"))

nj_prep_xg_pfos <- nj_prep_xg_pfos |>
  mutate(pop_density = case_when(popdensity2000 < 1000 ~ "<1000",
                                 popdensity2000 >= 1000 ~ "1000+"))

nj_prep_xg_pfna <- nj_prep_xg_pfna |>
  mutate(pop_density = case_when(popdensity2000 < 1000 ~ "<1000",
                                 popdensity2000 >= 1000 ~ "1000+"))


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
                        drop_vars = c("pfas_name","pfas_bin","year","year_group", "popdensity2000",
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
  
  # Geographic splits can leave categorical predictors with only one
  # training level. Such predictors cannot be dummy encoded. Record them
  # as dropped so the same columns are removed from the test data.
  single_level_vars <- names(x_input)[vapply(
    x_input,
    function(x) (is.factor(x) || is.character(x)) &&
      length(unique(x[!is.na(x)])) < 2L,
    logical(1)
  )]
  if (length(single_level_vars)) {
    message("Dropping single-level training predictors: ",
            paste(single_level_vars, collapse = ", "))
    x_input <- x_input %>%
      dplyr::select(-dplyr::any_of(single_level_vars))
    drop_vars <- unique(c(drop_vars, single_level_vars))
  }
  if (!ncol(x_input)) stop("No predictors remain after encoding exclusions.")
  
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
                              split_type = c("time", "geographic"),
                              year_var = "year",
                              density_var = "popdensity2000",
                              density_cutoff = 1000,
                              tune_grid_bin,
                              tune_grid_reg,
                              nfold = 5,
                              early_stopping_rounds = 25) {
  
  split_type <- match.arg(split_type)
  set.seed(seed)
  
  if (split_type == "time") {
    if (!year_var %in% names(df)) stop("Year column not found: ", year_var)
    years <- df[[year_var]]
    if (!is.numeric(years)) stop("The year column must be numeric/integer.")
    in_test <- !is.na(years) & years >= 2006 & years <= 2017
    in_train  <- !is.na(years) & years >= 2018 & years <= 2025
  } else {
    if (is.null(density_var) || !density_var %in% names(df)) {
      stop("Supply density_var as the name of the numeric people-per-square-mile column.")
    }
    density <- df[[density_var]]
    if (!is.numeric(density)) stop("The density column must be numeric.")
    in_train <- !is.na(density) & density < density_cutoff
    in_test  <- !is.na(density) & density > density_cutoff
    # Rows exactly at the cutoff (and missing density) belong to neither set.
  }
  
  idx_train <- which(in_train)
  idx_test <- which(in_test)
  if (!length(idx_train) || !length(idx_test)) {
    stop("The split produced an empty training or test set. Check the split variable and cutoff.")
  }
  df_train <- df[idx_train, , drop = FALSE]
  df_test  <- df[idx_test, , drop = FALSE]
  
  if (split_type == "geographic" && "PWSID" %in% names(df)) {
    shared_ids <- intersect(unique(df_train$PWSID), unique(df_test$PWSID))
    shared_ids <- shared_ids[!is.na(shared_ids)]
    if (length(shared_ids)) {
      warning(length(shared_ids), " PWSIDs occur in both density groups; ",
              "consider one density classification per PWSID for a site-level holdout.")
    }
  }
  message("Split: ", split_type, "; train = ", nrow(df_train),
          "; test = ", nrow(df_test),
          "; excluded = ", nrow(df) - nrow(df_train) - nrow(df_test))
  
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
    idx_test = idx_test,
    split_type = split_type,
    split_counts = c(train = nrow(df_train), test = nrow(df_test),
                     excluded = nrow(df) - nrow(df_train) - nrow(df_test)),
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


################################ TIME ##############################################


tune_grid_bin <- data.frame(
  nrounds = 400, eta = 0.08, max_depth = 4,
  gamma = 0.01, colsample_bytree = 0.8,
  min_child_weight = 2, subsample = 0.95,
  lambda = 0.5, alpha = 0
)

tune_grid_reg <- data.frame(
  nrounds = 250, eta = 0.045, max_depth = 6,
  gamma = 0.05, colsample_bytree = 0.8,
  min_child_weight = 1.5, subsample = 0.9,
  lambda = 1, alpha = 0
)


# Replace pfoa_data with your PFOA analysis data frame
#################### PFOA #############################

res_pfoa_time <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  split_type = "time",
  tune_grid_bin = tune_grid_bin,
  tune_grid_reg = tune_grid_reg
)

res_pfoa_time$split_counts
# train     test excluded 
# 6686     2071        0 


res_pfoa_time$stage1$auc_train # 0.9134486
res_pfoa_time$stage1$auc_test # 0.5576948
res_pfoa_time$stage1$min_cv_logloss # 0.4045748


res_pfoa_time$stage1$brier_train # 0.1180697
res_pfoa_time$stage1$brier_test # 0.3241192


res_pfoa_time$stage2$rmse_train # 0.3817276
res_pfoa_time$stage2$rmse_test # 1.198539
res_pfoa_time$stage2$r2_train # 0.669444
res_pfoa_time$stage2$r2_test # -1.40204



####################### TUNING ##############################
####################### STAGE 1 ##############################

# Use the one-row grids you already defined as starting settings.
base_bin <- tune_grid_bin[1, , drop = FALSE]
base_reg <- tune_grid_reg[1, , drop = FALSE]

make_grid_from_base <- function(base, vary) {
  vals <- list(
    nrounds = base$nrounds,
    eta = base$eta,
    gamma = base$gamma,
    colsample_bytree = base$colsample_bytree,
    min_child_weight = base$min_child_weight,
    subsample = base$subsample,
    lambda = base$lambda,
    alpha = base$alpha,
    max_depth = base$max_depth
  )
  for (nm in names(vary)) vals[[nm]] <- vary[[nm]]
  do.call(expand.grid, vals)
}

# Prepare the 2018–2025 training data once.
pfoa_train <- nj_prep_xg_pfoa[
  !is.na(nj_prep_xg_pfoa$year) &
    nj_prep_xg_pfoa$year >= 2018 &
    nj_prep_xg_pfoa$year <= 2025,
  ,
  drop = FALSE
]

encoder <- fit_encoder(pfoa_train)
X_train <- apply_encoder(pfoa_train, encoder)
y_bin <- ifelse(pfoa_train$pfas_bin_num == 1, 1, 0)

detect <- pfoa_train$pfas_bin_num == 1 &
  !is.na(pfoa_train$lnres)

X_train_det <- X_train[detect, , drop = FALSE]
y_reg <- pfoa_train$lnres[detect]

# Step 1: learning rate and tree size
bin_grid_1 <- make_grid_from_base(
  base_bin,
  list(
    nrounds = c(400, 600),
    eta = c(0.04, 0.08),
    max_depth = c(3, 4, 5),
    min_child_weight = c(2, 5)
  )
)

bin_step_1 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_1 <- bin_step_1[1, , drop = FALSE]

# Step 2: row and column sampling
bin_grid_2 <- make_grid_from_base(
  best_bin_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.7, 0.8, 0.9)
  )
)

bin_step_2 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_2 <- bin_step_2[1, , drop = FALSE]

# Step 3: regularization
bin_grid_3 <- make_grid_from_base(
  best_bin_2,
  list(
    gamma = c(0, 0.01, 0.05),
    lambda = c(0.5, 1),
    alpha = c(0, 1)
  )
)

bin_step_3 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_bin <- bin_step_3[1, , drop = FALSE]

best_bin
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
# 1     400 0.08     0              0.7                2       0.9      1     0         5          188 0.8936107  0.4035798


######### stage 2

reg_grid_1 <- make_grid_from_base(
  base_reg,
  list(
    nrounds = c(250, 500),
    eta = c(0.025, 0.045),
    max_depth = c(4, 5, 6),
    min_child_weight = c(1.5, 3)
  )
)

reg_step_1 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_1 <- reg_step_1[1, , drop = FALSE]

reg_grid_2 <- make_grid_from_base(
  best_reg_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.7, 0.8, 0.9)
  )
)

reg_step_2 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_2 <- reg_step_2[1, , drop = FALSE]

reg_grid_3 <- make_grid_from_base(
  best_reg_2,
  list(
    gamma = c(0, 0.05, 0.1),
    lambda = c(0.5, 1),
    alpha = c(0, 1)
  )
)

reg_step_3 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_reg <- reg_step_3[1, , drop = FALSE]

best_reg
# nrounds   eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds cv_rmse
# 1     500 0.025     0              0.7                3       0.8    0.5     0         6          177 0.41145


parameter_cols <- names(base_bin)

res_pfoa_time_tuned <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  split_type = "time",
  tune_grid_bin = best_bin[, parameter_cols, drop = FALSE],
  tune_grid_reg = best_reg[, parameter_cols, drop = FALSE]
)

res_pfoa_time_tuned$stage1$auc_train # 0.9137423
res_pfoa_time_tuned$stage1$auc_test # 0.5598221
res_pfoa_time_tuned$stage2$r2_train # 0.6648338
res_pfoa_time_tuned$stage2$r2_test # -1.326039
res_pfoa_time_tuned$stage2$rmse_train # 0.3843804
res_pfoa_time_tuned$stage2$rmse_test # 1.179426

######################## PFOS ######################################

############ main parameters ############################################
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

res_pfos_time <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  split_type = "time",
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos
)

res_pfos_time$split_counts
# train     test excluded 
# 6677     1954        0 


res_pfos_time$stage1$auc_train # 0.8999606
res_pfos_time$stage1$auc_test # 0.6200419
res_pfos_time$stage1$min_cv_logloss # 0.4267176


res_pfos_time$stage1$brier_train # 0.1263724
res_pfos_time$stage1$brier_test # 0.195565


res_pfos_time$stage2$rmse_train # 0.4063284
res_pfos_time$stage2$rmse_test # 1.131731
res_pfos_time$stage2$r2_train # 0.5681919
res_pfos_time$stage2$r2_test # -0.234278

####################### STAGE 1 ##############################

# Use the one-row grids you already defined as starting settings.
base_bin <- tune_grid_bin_pfos[1, , drop = FALSE]
base_reg <- tune_grid_reg_pfos[1, , drop = FALSE]

make_grid_from_base <- function(base, vary) {
  vals <- list(
    nrounds = base$nrounds,
    eta = base$eta,
    gamma = base$gamma,
    colsample_bytree = base$colsample_bytree,
    min_child_weight = base$min_child_weight,
    subsample = base$subsample,
    lambda = base$lambda,
    alpha = base$alpha,
    max_depth = base$max_depth
  )
  for (nm in names(vary)) vals[[nm]] <- vary[[nm]]
  do.call(expand.grid, vals)
}

# Prepare the 2018–2025 training data once.
pfos_train <- nj_prep_xg_pfos[
  !is.na(nj_prep_xg_pfos$year) &
    nj_prep_xg_pfos$year >= 2018 &
    nj_prep_xg_pfos$year <= 2025,
  ,
  drop = FALSE
]

encoder <- fit_encoder(pfos_train)
X_train <- apply_encoder(pfos_train, encoder)
y_bin <- ifelse(pfos_train$pfas_bin_num == 1, 1, 0)

detect <- pfos_train$pfas_bin_num == 1 &
  !is.na(pfos_train$lnres)

X_train_det <- X_train[detect, , drop = FALSE]
y_reg <- pfos_train$lnres[detect]

# Step 1: learning rate and tree size
bin_grid_1 <- make_grid_from_base(
  base_bin,
  list(
    nrounds = c(600, 800, 1000),
    eta = c(0.03, 0.04, 0.05),
    max_depth = c(5, 7),
    min_child_weight = c(1, 2, 3)
  )
)

bin_step_1 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_1 <- bin_step_1[1, , drop = FALSE]

# Step 2: row and column sampling
bin_grid_2 <- make_grid_from_base(
  best_bin_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.6, 0.8, 1)
  )
)

bin_step_2 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_2 <- bin_step_2[1, , drop = FALSE]

# Step 3: regularization
bin_grid_3 <- make_grid_from_base(
  best_bin_2,
  list(
    gamma = c(0, 0.05, 0.1),
    lambda = c(0, 1),
    alpha = c(0, 1)
  )
)

bin_step_3 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_bin <- bin_step_3[1, , drop = FALSE]

best_bin
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
# 1     600 0.04   0.1              0.8                3       0.9      1     1         7          161 0.8806237  0.4205016


######### stage 2

reg_grid_1 <- make_grid_from_base(
  base_reg,
  list(
    nrounds = c(600, 800, 1000),
    eta = c(0.025, 0.04),
    max_depth = c(4, 6),
    min_child_weight = c(1, 3, 5)
  )
)

reg_step_1 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_1 <- reg_step_1[1, , drop = FALSE]

reg_grid_2 <- make_grid_from_base(
  best_reg_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.7, 0.85, 1)
  )
)

reg_step_2 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_2 <- reg_step_2[1, , drop = FALSE]

reg_grid_3 <- make_grid_from_base(
  best_reg_2,
  list(
    gamma = c(0, 0.1, 0.2),
    lambda = c(0, 1),
    alpha = c(0, 1)
  )
)

reg_step_3 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_reg <- reg_step_3[1, , drop = FALSE]

best_reg
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds   cv_rmse
# 1     600 0.04   0.2                1                5       0.9      1     0         4          139 0.4421521

parameter_cols <- names(base_bin)

res_pfos_time_tuned <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  split_type = "time",
  tune_grid_bin = best_bin[, parameter_cols, drop = FALSE],
  tune_grid_reg = best_reg[, parameter_cols, drop = FALSE]
)

res_pfos_time_tuned$stage1$auc_train # 0.9003504
res_pfos_time_tuned$stage1$auc_test # 0.6325936
res_pfos_time_tuned$stage2$r2_train # 0.5493653
res_pfos_time_tuned$stage2$r2_test # -0.2179506
res_pfos_time_tuned$stage2$rmse_train # 0.4150917
res_pfos_time_tuned$stage2$rmse_test # 1.124221


#################### PFNA #####################################################

# pfna time

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

res_pfna_time <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  split_type = "time",
  tune_grid_bin = tune_grid_bin_pfna,
  tune_grid_reg = tune_grid_reg_pfna
)

res_pfna_time$split_counts
#   train     test excluded 
#    6875     1769        0 


res_pfna_time$stage1$auc_train # 0.9299092
res_pfna_time$stage1$auc_test # 0.82516
res_pfna_time$stage1$min_cv_logloss # 0.1750841


res_pfna_time$stage1$brier_train # 0.04702558
res_pfna_time$stage1$brier_test # 0.06944888


res_pfna_time$stage2$rmse_train # 0.4477489 
res_pfna_time$stage2$rmse_test # 0.8770229
res_pfna_time$stage2$r2_train # 0.7941887
res_pfna_time$stage2$r2_test # -0.1558687


####################### STAGE 1 ##############################

# Use the one-row grids you already defined as starting settings.
base_bin <- tune_grid_bin_pfna[1, , drop = FALSE]
base_reg <- tune_grid_reg_pfna[1, , drop = FALSE]

make_grid_from_base <- function(base, vary) {
  vals <- list(
    nrounds = base$nrounds,
    eta = base$eta,
    gamma = base$gamma,
    colsample_bytree = base$colsample_bytree,
    min_child_weight = base$min_child_weight,
    subsample = base$subsample,
    lambda = base$lambda,
    alpha = base$alpha,
    max_depth = base$max_depth
  )
  for (nm in names(vary)) vals[[nm]] <- vary[[nm]]
  do.call(expand.grid, vals)
}

# Prepare the 2018–2025 training data once.
pfna_train <- nj_prep_xg_pfna[
  !is.na(nj_prep_xg_pfna$year) &
    nj_prep_xg_pfna$year >= 2018 &
    nj_prep_xg_pfna$year <= 2025,
  ,
  drop = FALSE
]

encoder <- fit_encoder(pfna_train)
X_train <- apply_encoder(pfna_train, encoder)
y_bin <- ifelse(pfna_train$pfas_bin_num == 1, 1, 0)

detect <- pfna_train$pfas_bin_num == 1 &
  !is.na(pfna_train$lnres)

X_train_det <- X_train[detect, , drop = FALSE]
y_reg <- pfna_train$lnres[detect]

# Step 1: learning rate and tree size
bin_grid_1 <- make_grid_from_base(
  base_bin,
  list(
    nrounds = c(600, 800, 1000),
    eta = c(0.03, 0.04, 0.05),
    max_depth = c(3, 4, 5),
    min_child_weight = c(1, 3, 5)
  )
)

bin_step_1 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_1 <- bin_step_1[1, , drop = FALSE]

# Step 2: row and column sampling
bin_grid_2 <- make_grid_from_base(
  best_bin_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.6, 0.8, 1)
  )
)

bin_step_2 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_2 <- bin_step_2[1, , drop = FALSE]

# Step 3: regularization
bin_grid_3 <- make_grid_from_base(
  best_bin_2,
  list(
    gamma = c(0, 0.05, 0.1),
    lambda = c(0, 1),
    alpha = c(0, 1)
  )
)

bin_step_3 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_pfna <- bin_step_3[1, , drop = FALSE]

best_bin_pfna
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
# 1    1000 0.03   0.1                1                1       0.8      0     0         3          419 0.9123627  0.1720669

######### stage 2

reg_grid_1 <- make_grid_from_base(
  base_reg,
  list(
    nrounds = c(600, 800, 1000),
    eta = c(0.02, 0.03, 0.04),
    max_depth = c(4, 6, 8),
    min_child_weight = c(1, 3, 4, 5)
  )
)

reg_step_1 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_1 <- reg_step_1[1, , drop = FALSE]

reg_grid_2 <- make_grid_from_base(
  best_reg_1,
  list(
    subsample = c(0.7, 0.8, 0.9),
    colsample_bytree = c(0.8, 0.9, 1.0)
  )
)

reg_step_2 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_2 <- reg_step_2[1, , drop = FALSE]

reg_grid_3 <- make_grid_from_base(
  best_reg_2,
  list(
    gamma = c(0, 0.05, 0.1),
    lambda = c(0, 1, 2),
    alpha = c(0, 0.5, 1)
  )
)

reg_step_3 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_pfna <- reg_step_3[1, , drop = FALSE]

best_reg_pfna
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds   cv_rmse
#    1000 0.02   0.1                1                1       0.7      0     0         4          329 0.5134762


parameter_cols <- names(base_bin)

res_pfna_time_tuned <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  split_type = "time",
  tune_grid_bin = best_bin_pfna[, parameter_cols, drop = FALSE],
  tune_grid_reg = best_reg_pfna[, parameter_cols, drop = FALSE]
)

res_pfna_time_tuned$stage1$auc_train # 0.9280796
res_pfna_time_tuned$stage1$auc_test # 0.8210595
res_pfna_time_tuned$stage2$r2_train # 0.7919894
res_pfna_time_tuned$stage2$r2_test # -0.1344439
res_pfna_time_tuned$stage2$rmse_train # 0.4501348
res_pfna_time_tuned$stage2$rmse_test # 0.8688568

# Split: time; train = 6875; test = 1769; excluded = 0




################################ GEOGRAPHIC ##############################################


tune_grid_bin <- data.frame(
  nrounds = 400, eta = 0.08, max_depth = 4,
  gamma = 0.01, colsample_bytree = 0.8,
  min_child_weight = 2, subsample = 0.95,
  lambda = 0.5, alpha = 0
)

tune_grid_reg <- data.frame(
  nrounds = 250, eta = 0.045, max_depth = 6,
  gamma = 0.05, colsample_bytree = 0.8,
  min_child_weight = 1.5, subsample = 0.9,
  lambda = 1, alpha = 0
)


res_pfoa_geographic <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  split_type = "geographic",
  density_var = "popdensity2000",
  tune_grid_bin = tune_grid_bin,
  tune_grid_reg = tune_grid_reg
)


res_pfoa_geographic$split_counts
# train     test excluded 
# 5252     3505        0 



res_pfoa_geographic$stage1$auc_train # 0.9320856
res_pfoa_geographic$stage1$auc_test # 0.7416784
res_pfoa_geographic$stage1$min_cv_logloss # 0.3723137


res_pfoa_geographic$stage1$brier_train # 0.1055646
res_pfoa_geographic$stage1$brier_test # 0.2654442


res_pfoa_geographic$stage2$rmse_train # 0.3982559
res_pfoa_geographic$stage2$rmse_test # 0.932461
res_pfoa_geographic$stage2$r2_train # 0.6263337
res_pfoa_geographic$stage2$r2_test # -0.5110383


######################### TUNING #######################################
####################### STAGE 1 ##############################

# Use the one-row grids you already defined as starting settings.
base_bin_pfoa_geo <- tune_grid_bin[1, , drop = FALSE]
base_reg_pfoa_geo <- tune_grid_reg[1, , drop = FALSE]

make_grid_from_base <- function(base, vary) {
  vals <- list(
    nrounds = base$nrounds,
    eta = base$eta,
    gamma = base$gamma,
    colsample_bytree = base$colsample_bytree,
    min_child_weight = base$min_child_weight,
    subsample = base$subsample,
    lambda = base$lambda,
    alpha = base$alpha,
    max_depth = base$max_depth
  )
  for (nm in names(vary)) vals[[nm]] <- vary[[nm]]
  do.call(expand.grid, vals)
}

# Prepare the population density <1000 training data.
pfoa_train_geo <- nj_prep_xg_pfoa[
  !is.na(nj_prep_xg_pfoa$pop_density) &
    nj_prep_xg_pfoa$pop_density == "<1000",
  ,
  drop = FALSE
]

encoder <- fit_encoder(pfoa_train_geo)
X_train <- apply_encoder(pfoa_train_geo, encoder)
y_bin <- ifelse(pfoa_train_geo$pfas_bin_num == 1, 1, 0)

detect <- pfoa_train_geo$pfas_bin_num == 1 &
  !is.na(pfoa_train_geo$lnres)

X_train_det <- X_train[detect, , drop = FALSE]
y_reg <- pfoa_train_geo$lnres[detect]

# Step 1: learning rate and tree size
bin_grid_1 <- make_grid_from_base(
  base_bin,
  list(
    nrounds = c(200, 400, 600),
    eta = c(0.04, 0.06, 0.08),
    max_depth = c(2, 4, 6),
    min_child_weight = c(1, 2, 4)
  )
)

bin_step_1 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_1 <- bin_step_1[1, , drop = FALSE]

# Step 2: row and column sampling
bin_grid_2 <- make_grid_from_base(
  best_bin_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.7, 0.8, 0.9)
  )
)

bin_step_2 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_2 <- bin_step_2[1, , drop = FALSE]

# Step 3: regularization
bin_grid_3 <- make_grid_from_base(
  best_bin_2,
  list(
    gamma = c(0, 0.01, 0.05),
    lambda = c(0.5, 1),
    alpha = c(0, 1)
  )
)

bin_step_3 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_pfoa_geo <- bin_step_3[1, , drop = FALSE]

best_bin_pfoa_geo
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
# 1     400 0.08  0.01              0.9                4         1      1     1         4          400 0.9115873  0.3716619


######### stage 2

reg_grid_1 <- make_grid_from_base(
  base_reg,
  list(
    nrounds = c(250, 350, 500),
    eta = c(0.025, 0.045, 0.6),
    max_depth = c(4, 5, 6),
    min_child_weight = c(1.5, 3)
  )
)

reg_step_1 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_1 <- reg_step_1[1, , drop = FALSE]

reg_grid_2 <- make_grid_from_base(
  best_reg_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.7, 0.8, 0.9)
  )
)

reg_step_2 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_2 <- reg_step_2[1, , drop = FALSE]

reg_grid_3 <- make_grid_from_base(
  best_reg_2,
  list(
    gamma = c(0, 0.05, 0.1),
    lambda = c(0.5, 1),
    alpha = c(0, 1)
  )
)

reg_step_3 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_pfoa_geo <- reg_step_3[1, , drop = FALSE]

best_reg_pfoa_geo
# nrounds   eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds   cv_rmse
# 1     250 0.025     0              0.7                3       0.8    0.5     0         6          185 0.4293512


parameter_cols <- names(base_bin)

res_pfoa_geo_tuned <- run_two_phase_xgb(
  df = nj_prep_xg_pfoa,
  split_type = "geographic",
  tune_grid_bin = best_bin_pfoa_geo[, parameter_cols, drop = FALSE],
  tune_grid_reg = best_reg_pfoa_geo[, parameter_cols, drop = FALSE]
)

# Split: geographic; train = 5252; test = 3505; excluded = 0

res_pfoa_geo_tuned$stage1$auc_train # 0.9325596
res_pfoa_geo_tuned$stage1$auc_test # 0.7449056
res_pfoa_geo_tuned$stage2$r2_train # 0.6273497
res_pfoa_geo_tuned$stage2$r2_test # -0.4272798
res_pfoa_geo_tuned$stage2$rmse_train # 0.3977141
res_pfoa_geo_tuned$stage2$rmse_test # 0.9062489



########################################### PFOS ################################################################################
# pfos geographic
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

res_pfos_geographic <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  split_type = "geographic",
  density_var = "popdensity2000",
  tune_grid_bin = tune_grid_bin_pfos,
  tune_grid_reg = tune_grid_reg_pfos
)

res_pfos_geographic$split_counts
# train     test excluded 
# 5204     3427        0 



res_pfos_geographic$stage1$auc_train # 0.9032199
res_pfos_geographic$stage1$auc_test # 0.6986058
res_pfos_geographic$stage1$min_cv_logloss # 0.3822625 


res_pfos_geographic$stage1$brier_train # 0.1128566
res_pfos_geographic$stage1$brier_test #  0.3263921


res_pfos_geographic$stage2$rmse_train # 0.4055746
res_pfos_geographic$stage2$rmse_test # 0.691053
res_pfos_geographic$stage2$r2_train # 0.5429994
res_pfos_geographic$stage2$r2_test # -0.05208473

################################# TUNING #####################################
####################### STAGE 1 ##############################

# Use the one-row grids you already defined as starting settings.
base_bin <- tune_grid_bin_pfos[1, , drop = FALSE]
base_reg <- tune_grid_reg_pfos[1, , drop = FALSE]

make_grid_from_base <- function(base, vary) {
  vals <- list(
    nrounds = base$nrounds,
    eta = base$eta,
    gamma = base$gamma,
    colsample_bytree = base$colsample_bytree,
    min_child_weight = base$min_child_weight,
    subsample = base$subsample,
    lambda = base$lambda,
    alpha = base$alpha,
    max_depth = base$max_depth
  )
  for (nm in names(vary)) vals[[nm]] <- vary[[nm]]
  do.call(expand.grid, vals)
}

# Prepare the population density <1000 training data.
pfos_train_geo <- nj_prep_xg_pfos[
  !is.na(nj_prep_xg_pfos$pop_density) &
    nj_prep_xg_pfos$pop_density == "<1000",
  ,
  drop = FALSE
]
encoder <- fit_encoder(pfos_train_geo)
X_train <- apply_encoder(pfos_train_geo, encoder)
y_bin <- ifelse(pfos_train_geo$pfas_bin_num == 1, 1, 0)

detect <- pfos_train_geo$pfas_bin_num == 1 &
  !is.na(pfos_train_geo$lnres)

X_train_det <- X_train[detect, , drop = FALSE]
y_reg <- pfos_train_geo$lnres[detect]

# Step 1: learning rate and tree size
bin_grid_1 <- make_grid_from_base(
  base_bin,
  list(
    nrounds = c(600, 800, 1000),
    eta = c(0.03, 0.04, 0.05),
    max_depth = c(5, 7),
    min_child_weight = c(1, 2, 3)
  )
)

bin_step_1 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_1 <- bin_step_1[1, , drop = FALSE]

# Step 2: row and column sampling
bin_grid_2 <- make_grid_from_base(
  best_bin_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.6, 0.8, 1)
  )
)

bin_step_2 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_2 <- bin_step_2[1, , drop = FALSE]

# Step 3: regularization
bin_grid_3 <- make_grid_from_base(
  best_bin_2,
  list(
    gamma = c(0, 0.05, 0.1),
    lambda = c(0, 1),
    alpha = c(0, 1)
  )
)

bin_step_3 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_pfos_geo <- bin_step_3[1, , drop = FALSE]

best_bin_pfos_geo
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
# 1     800 0.05     0              0.8                3         1      1     1         5          264 0.8786539  0.3764921



######### stage 2

reg_grid_1 <- make_grid_from_base(
  base_reg,
  list(
    nrounds = c(600, 800, 1000),
    eta = c(0.025, 0.04),
    max_depth = c(4, 6),
    min_child_weight = c(1, 3, 5)
  )
)

reg_step_1 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_1 <- reg_step_1[1, , drop = FALSE]

reg_grid_2 <- make_grid_from_base(
  best_reg_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.7, 0.85, 1)
  )
)

reg_step_2 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_2 <- reg_step_2[1, , drop = FALSE]

reg_grid_3 <- make_grid_from_base(
  best_reg_2,
  list(
    gamma = c(0, 0.1, 0.2),
    lambda = c(0, 1),
    alpha = c(0, 1)
  )
)

reg_step_3 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_pfos_geo <- reg_step_3[1, , drop = FALSE]

best_reg_pfos_geo
# nrounds   eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds   cv_rmse
# 1     800 0.025   0.1             0.85                1       0.8      0     0         6          165 0.4401054



parameter_cols <- names(base_bin)

res_pfos_geo_tuned <- run_two_phase_xgb(
  df = nj_prep_xg_pfos,
  split_type = "geographic",
  tune_grid_bin = best_bin_pfos_geo[, parameter_cols, drop = FALSE],
  tune_grid_reg = best_reg_pfos_geo[, parameter_cols, drop = FALSE]
)

# Split: geographic; train = 5204; test = 3427; excluded = 0

res_pfos_geo_tuned$stage1$auc_train # 0.901067
res_pfos_geo_tuned$stage1$auc_test # 0.7643506
res_pfos_geo_tuned$stage2$r2_train # 0.5539834
res_pfos_geo_tuned$stage2$r2_test # 0.00492707
res_pfos_geo_tuned$stage2$rmse_train # 0.400671
res_pfos_geo_tuned$stage2$rmse_test # 0.6720683






######################### PFNA geographic ########################### 
####### original parameters #####################

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

res_pfna_geographic <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  split_type = "geographic",
  density_var = "popdensity2000",
  tune_grid_bin = tune_grid_bin_pfna,
  tune_grid_reg = tune_grid_reg_pfna
)

res_pfna_geographic$split_counts
#    train     test excluded 
#    5200     3444        0 


res_pfna_geographic$stage1$auc_train # 0.954125
res_pfna_geographic$stage1$auc_test # 0.6916039
res_pfna_geographic$stage1$min_cv_logloss # 0.1256843


res_pfna_geographic$stage1$brier_train # 0.03121274
res_pfna_geographic$stage1$brier_test # 0.1035199


res_pfna_geographic$stage2$rmse_train  # 0.5464541
res_pfna_geographic$stage2$rmse_test # 1.234835
res_pfna_geographic$stage2$r2_train # 0.8017575
res_pfna_geographic$stage2$r2_test # -0.4469595


################################# TUNING #####################################
####################### STAGE 1 ##############################

# Use the one-row grids you already defined as starting settings.
base_bin <- tune_grid_bin_pfna[1, , drop = FALSE]
base_reg <- tune_grid_reg_pfna[1, , drop = FALSE]

make_grid_from_base <- function(base, vary) {
  vals <- list(
    nrounds = base$nrounds,
    eta = base$eta,
    gamma = base$gamma,
    colsample_bytree = base$colsample_bytree,
    min_child_weight = base$min_child_weight,
    subsample = base$subsample,
    lambda = base$lambda,
    alpha = base$alpha,
    max_depth = base$max_depth
  )
  for (nm in names(vary)) vals[[nm]] <- vary[[nm]]
  do.call(expand.grid, vals)
}

# Prepare the population density <1000 training data.
pfna_train_geo <- nj_prep_xg_pfna[
  !is.na(nj_prep_xg_pfna$pop_density) &
    nj_prep_xg_pfna$pop_density == "<1000",
  ,
  drop = FALSE
]

encoder <- fit_encoder(pfna_train)
X_train <- apply_encoder(pfna_train, encoder)
y_bin <- ifelse(pfna_train$pfas_bin_num == 1, 1, 0)

detect <- pfna_train$pfas_bin_num == 1 &
  !is.na(pfna_train$lnres)

X_train_det <- X_train[detect, , drop = FALSE]
y_reg <- pfna_train$lnres[detect]

# Step 1: learning rate and tree size
bin_grid_1 <- make_grid_from_base(
  base_bin,
  list(
    nrounds = c(600, 800, 1000),
    eta = c(0.03, 0.04, 0.05),
    max_depth = c(3, 4, 5),
    min_child_weight = c(1, 3, 5)
  )
)

bin_step_1 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_1 <- bin_step_1[1, , drop = FALSE]

# Step 2: row and column sampling
bin_grid_2 <- make_grid_from_base(
  best_bin_1,
  list(
    subsample = c(0.8, 0.9, 1),
    colsample_bytree = c(0.6, 0.8, 1)
  )
)

bin_step_2 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_2 <- bin_step_2[1, , drop = FALSE]

# Step 3: regularization
bin_grid_3 <- make_grid_from_base(
  best_bin_2,
  list(
    gamma = c(0, 0.05, 0.1),
    lambda = c(0, 1),
    alpha = c(0, 1)
  )
)

bin_step_3 <- tune_xgb_binary(
  X_train, y_bin, bin_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_bin_pfna_geo <- bin_step_3[1, , drop = FALSE]

best_bin_pfna_geo
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds    cv_auc cv_logloss
# 1    1000 0.03   0.1                1                1       0.8      0     0         3          419 0.9123627  0.1720669



######### stage 2

reg_grid_1 <- make_grid_from_base(
  base_reg,
  list(
    nrounds = c(600, 800, 1000),
    eta = c(0.02, 0.03, 0.04),
    max_depth = c(4, 6, 8),
    min_child_weight = c(1, 3, 4, 5)
  )
)

reg_step_1 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_1,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_1 <- reg_step_1[1, , drop = FALSE]

reg_grid_2 <- make_grid_from_base(
  best_reg_1,
  list(
    subsample = c(0.7, 0.8, 0.9),
    colsample_bytree = c(0.8, 0.9, 1.0)
  )
)

reg_step_2 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_2,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_2 <- reg_step_2[1, , drop = FALSE]

reg_grid_3 <- make_grid_from_base(
  best_reg_2,
  list(
    gamma = c(0, 0.05, 0.1),
    lambda = c(0, 1, 2),
    alpha = c(0, 0.5, 1)
  )
)

reg_step_3 <- tune_xgb_regression(
  X_train_det, y_reg, reg_grid_3,
  nfold = 5, early_stopping_rounds = 25
)
best_reg_pfna_geo <- reg_step_3[1, , drop = FALSE]

best_reg_pfna_geo
# nrounds  eta gamma colsample_bytree min_child_weight subsample lambda alpha max_depth best_nrounds   cv_rmse
# 1    1000 0.02   0.1                1                1       0.7      0     0         4          329 0.5134762


parameter_cols <- names(base_bin)

res_pfna_geo_tuned <- run_two_phase_xgb(
  df = nj_prep_xg_pfna,
  split_type = "geographic",
  tune_grid_bin = best_bin_pfna_geo[, parameter_cols, drop = FALSE],
  tune_grid_reg = best_reg_pfna_geo[, parameter_cols, drop = FALSE]
)

res_pfna_geo_tuned$stage1$auc_train # 0.9496076
res_pfna_geo_tuned$stage1$auc_test # 0.5998436
res_pfna_geo_tuned$stage2$r2_train # 0.7920697
res_pfna_geo_tuned$stage2$r2_test # -0.2524294
res_pfna_geo_tuned$stage2$rmse_train # 0.5596472
res_pfna_geo_tuned$stage2$rmse_test # 1.148835

# Split: geographic; train = 5200; test = 3444; excluded = 0





