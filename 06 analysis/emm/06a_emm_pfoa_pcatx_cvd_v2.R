#################################################################
# R Code - EMM

# Hari Iyer
# Date: 09/23/2026

# Goal: Effect measure modification between pfoa and cvd for treatment
# Analyses here for effect modification on the multiplicative
# and additive scale

# Need to incorporate into the imputation models
# Run these for all PFAS, all outcomes using the original IQR change
# Explore additive
#################################################################


library(tidyverse)
library(mice)
library(survival)
library(forcats)
library(broom)
library(dplyr)
library(ggplot2)
library(sandwich)


load("/extra/home/hi97/R/registry/nj_res_hx/long_list_imp.Rda")

colnames(long_imp_list[[1]])

levels(long_imp_list[[1]]$txsurg.f)
levels(long_imp_list[[1]]$txrad.f)
levels(long_imp_list[[1]]$txhorm.f)
levels(long_imp_list[[1]]$psa.f)
levels(long_imp_list[[1]]$grade_combo.f)
contrasts(long_imp_list[[1]]$stage.f)



# adding 65 age variable 
long_imp_list <- lapply(
  long_imp_list,
  function(dat) {
    dat %>%
      mutate(
        dxage65 = factor(
          case_when(
            dxage < 65  ~ "<65",
            dxage >= 65 ~ "65+",
            TRUE ~ NA_character_
          ),
          levels = c("<65", "65+")
        )
      )
  }
)

colnames(long_imp_list)

#####################################################################
# PFOA -> CVD mortality
# Effect modification by:
#   1. Race/ethnicity
#   2. Age at diagnosis
#   3. Stage
#####################################################################

exposure <- "pfoa_iqr"
outcome  <- "cvd_death"

effect_modifiers <- c(
  surg  = "txsurg.f",
  rad   = "txrad.f",
  horm = "txhorm.f"
)


#####################################################################
# Function
#####################################################################

run_pfoa_interaction <- function(effect_modifier,
                                 modifier_name,
                                 long_imp_list) {
  
  message("Running effect modification by: ", modifier_name)
  
  # ---------------------------------------------------------------
  # 1. Covariates
  #
  # Remove the effect modifier from the adjustment covariates,
  # because it will be included explicitly in the interaction term.
  # ---------------------------------------------------------------
  
  covariates <- c(
    "age_cat.f",
    "dxyr_cat.f",
    "nhrace",
    "married.f",
    "insur.f",
    "stage.f",
    "txsurg.f",
    "txrad.f",
    "psa.f",
    "nses_quin",
    "no2_new",
    "pm25_new",
    "ndvi27_final"
  )
  
  # Remove the effect modifier from the ordinary covariates
  covariates <- setdiff(covariates, effect_modifier)
  
  # If binary age is the effect modifier,
  # also remove the original 5-category age variable
  if (effect_modifier == "dxage65") {
    covariates <- setdiff(covariates, "age_cat.f")
  }
  
  # ---------------------------------------------------------------
  # 2. Main-effects model
  # ---------------------------------------------------------------
  
  base_formula <- as.formula(
    paste0(
      "Surv(tstart, tstop, ", outcome, ") ~ ",
      exposure, " + ",
      effect_modifier, " + ",
      paste(covariates, collapse = " + ")
    )
  )
  
  
  fit_base <- lapply(long_imp_list, function(dat) {
    
    coxph(
      base_formula,
      data = dat
    )
    
  })
  
  fit_base_mira <- as.mira(fit_base)
  
  
  # ---------------------------------------------------------------
  # 3. Interaction model
  # ---------------------------------------------------------------
  
  interaction_formula <- as.formula(
    paste0(
      "Surv(tstart, tstop, ", outcome, ") ~ ",
      exposure, " * ", effect_modifier, " + ",
      paste(covariates, collapse = " + ")
    )
  )
  
  
  fit_int <- lapply(long_imp_list, function(dat) {
    
    coxph(
      interaction_formula,
      data = dat
    )
    
  })
  
  fit_int_mira <- as.mira(fit_int)
  
  
  # ---------------------------------------------------------------
  # 4. Overall P-interaction
  # ---------------------------------------------------------------
  
  loglik_base <- sapply(
    fit_base_mira$analyses,
    function(model) as.numeric(logLik(model))
  )
  
  loglik_int <- sapply(
    fit_int_mira$analyses,
    function(model) as.numeric(logLik(model))
  )
  
  
  LRT_stat <- 2 * (loglik_int - loglik_base)
  
  
  # Number of interaction terms
  #
  # This automatically determines the df rather than hard-coding 2.
  
  n_int_terms <- sum(
    grepl(
      ":",
      names(coef(fit_int_mira$analyses[[1]]))
    )
  )
  
  
  lrt_result <- miceadds::micombine.chisquare(
    LRT_stat,
    df = n_int_terms
  )
  
  
  p_interaction <- lrt_result[2]
  
  
  # ---------------------------------------------------------------
  # 5. Full pooled interaction model
  # ---------------------------------------------------------------
  
  pooled <- summary(pool(fit_int_mira))
  
  
  pooled <- pooled %>%
    mutate(
      HR = exp(estimate),
      HR_95LL = exp(estimate - 1.96 * std.error),
      HR_95UL = exp(estimate + 1.96 * std.error),
      effect_modifier = modifier_name,
      p_interaction = p_interaction
    )
  
  
  # ---------------------------------------------------------------
  # Return everything
  # ---------------------------------------------------------------
  
  list(
    
    modifier = modifier_name,
    
    base_model = fit_base_mira,
    
    interaction_model = fit_int_mira,
    
    pooled_model = pooled,
    
    p_interaction = p_interaction,
    
    df_interaction = n_int_terms
    
  )
}


#####################################################################
# Run all 3 effect modification analyses
#####################################################################

pfoa_cvd_interactions <- lapply(
  names(effect_modifiers),
  function(x) {
    
    run_pfoa_interaction(
      effect_modifier = effect_modifiers[[x]],
      modifier_name = x,
      long_imp_list = long_imp_list
    )
    
  }
)

names(pfoa_cvd_interactions) <- names(effect_modifiers)



#####################################################################
# P-interaction results
#####################################################################

p_interaction_results <- bind_rows(
  lapply(
    names(pfoa_cvd_interactions),
    function(x) {
      
      data.frame(
        effect_modifier = x,
        df = pfoa_cvd_interactions[[x]]$df_interaction,
        p_interaction =
          pfoa_cvd_interactions[[x]]$p_interaction
      )
      
    }
  )
)

p_interaction_results

write.csv(
  p_interaction_results,
  "~/p_interaction_results.csv",
  row.names = FALSE
)

 save(
  pfoa_cvd_interactions,
  file = "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/output/pfoa_cvd_interactions.Rda"
)


#####################################################################
# Function to calculate and pool linear contrasts from Cox models
#####################################################################

pool_contrast_cox <- function(fit_mira, terms) {
  
  # Extract models from each imputed dataset
  models <- fit_mira$analyses
  
  # Function for one imputed model
  get_contrast <- function(model) {
    
    beta <- coef(model)
    V <- vcov(model)
    
    # Create contrast vector
    L <- rep(0, length(beta))
    names(L) <- names(beta)
    
    # terms should be coefficient names
    if (!all(terms %in% names(beta))) {
      stop(
        "These terms were not found in the model: ",
        paste(terms[!terms %in% names(beta)], collapse = ", ")
      )
    }
    
    L[terms] <- 1
    
    # Linear combination of coefficients
    estimate <- sum(L * beta)
    
    # Variance of linear combination
    variance <- as.numeric(
      t(L) %*% V %*% L
    )
    
    data.frame(
      estimate = estimate,
      variance = variance
    )
  }
  
  # Calculate contrast in every imputation
  contrast_list <- lapply(models, get_contrast)
  
  estimates <- sapply(contrast_list, function(x) x$estimate)
  variances <- sapply(contrast_list, function(x) x$variance)
  
  # Pool using Rubin's rules
  pooled <- mice::pool.scalar(
    Q = estimates,
    U = variances
  )
  
  estimate <- pooled$qbar
  se <- sqrt(pooled$t)
  
  # HR and 95% CI
  HR <- exp(estimate)
  HR_95LL <- exp(estimate - 1.96 * se)
  HR_95UL <- exp(estimate + 1.96 * se)
  
  # Wald test
  z <- estimate / se
  p_value <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  data.frame(
    estimate = estimate,
    std.error = se,
    HR = HR,
    HR_95LL = HR_95LL,
    HR_95UL = HR_95UL,
    p.value = p_value
  )
}


#####################################################################
# TX Surg: Yes vs No
#####################################################################

fit_surg <- pfoa_cvd_interactions$surg$interaction_model

# No - reference group
txsurg_No <- pool_contrast_cox(
  fit_surg,
  terms = c(
    "pfoa_iqr"
  )
)

# Yes
txsurg_Yes <- pool_contrast_cox(
  fit_surg,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:txsurg.fYes"
  )
)

surg_results <- bind_rows(
  "No surg" = txsurg_No,
  "Surg" = txsurg_Yes,
  .id = "subgroup"
) |>
  mutate(
    effect_modifier = "Surger",
    p_interaction =
      pfoa_cvd_interactions$surg$p_interaction
  )

surg_results


#####################################################################
# TX Radiation: Yes vs No
#####################################################################

fit_rad <- pfoa_cvd_interactions$rad$interaction_model

# No - reference group
txrad_No <- pool_contrast_cox(
  fit_rad,
  terms = c(
    "pfoa_iqr"
  )
)

# Yes
txrad_Yes <- pool_contrast_cox(
  fit_rad,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:txrad.fYes"
  )
)

rad_results <- bind_rows(
  "No rad" = txrad_No,
  "rad" = txrad_Yes,
  .id = "subgroup"
) |>
  mutate(
    effect_modifier = "Rad",
    p_interaction =
      pfoa_cvd_interactions$rad$p_interaction
  )

rad_results



#####################################################################
# TX Hormone: Yes vs No
#####################################################################

fit_horm <- pfoa_cvd_interactions$horm$interaction_model

# No - reference group
txhorm_No <- pool_contrast_cox(
  fit_horm,
  terms = c(
    "pfoa_iqr"
  )
)

# Yes
txhorm_Yes <- pool_contrast_cox(
  fit_horm,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:txhorm.fYes"
  )
)

horm_results <- bind_rows(
  "No horm" = txhorm_No,
  "horm" = txhorm_Yes,
  .id = "subgroup"
) |>
  mutate(
    effect_modifier = "horm",
    p_interaction =
      pfoa_cvd_interactions$horm$p_interaction
  )

horm_results
names(coef(fit_horm$analyses[[1]]))



#####################################################################
# COMBINE RESULTS
#####################################################################

pfoa_cvd_effect_modification <- bind_rows(
  surg_results,
  rad_results,
  horm_results
)

pfoa_cvd_effect_modification

save(
  pfoa_cvd_interactions,
  file = "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/pfoa_cvd_interactions.Rda"
)

write.csv(
  pfoa_cvd_effect_modification,
  "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/output/pfoa_cvd_effect_modification.csv",
  row.names = FALSE
)


###########################################################################
# Absolute rate differences by PFAS, Mortality for ADT interactions
###########################################################################

library(dplyr)
library(tidyr)
library(purrr)

make_rate_diff <- function(data, pfas_var, mortality_var) {
  
  data %>%
    filter(
      stage.f != "Distant",
      !(stage.f == "Localized" &
          grade_combo.f == "<7" &
          psa.f == "<10")
    ) %>%
    mutate(
      interval_py = (tstop - tstart) / 12,
      
      pfas_high = ifelse(
        .data[[pfas_var]] >= median(.data[[pfas_var]], na.rm = TRUE),
        "High PFAS",
        "Low PFAS"
      )
    ) %>%
    group_by(txhorm.f, pfas_high) %>%
    summarise(
      n_patients = n_distinct(patient_study_id),
      
      deaths = n_distinct(
        patient_study_id[.data[[mortality_var]] == 1]
      ),
      
      person_years = sum(interval_py, na.rm = TRUE),
      
      rate_per_1000py = (deaths / person_years) * 1000,
      
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = pfas_high,
      values_from = c(
        n_patients,
        deaths,
        person_years,
        rate_per_1000py
      )
    ) %>%
    mutate(
      rate_difference =
        `rate_per_1000py_High PFAS` -
        `rate_per_1000py_Low PFAS`,
      
      pfas = pfas_var,
      outcome = mortality_var
    ) %>%
    select(
      pfas,
      outcome,
      txhorm.f,
      everything()
    )
}

pfas_vars <- c(
  "pfoa_dx",
  "pfos_dx",
  "pfna_dx"
)

mortality_vars <- c(
  "all_death",
  "prostate_death",
  "cvd_death"
)

rate_diff_results <- crossing(
  pfas_var = pfas_vars,
  mortality_var = mortality_vars
) %>%
  mutate(
    results = map2(
      pfas_var,
      mortality_var,
      ~ make_rate_diff(
        data = long_imp_list[[1]],
        pfas_var = .x,
        mortality_var = .y
      )
    )
  ) %>%
  select(results) %>%
  unnest(results)

interaction_rd <- rate_diff_results %>%
  select(
    pfas,
    outcome,
    txhorm.f,
    rate_difference
  ) %>%
  pivot_wider(
    names_from = txhorm.f,
    values_from = rate_difference
  ) %>%
  mutate(
    difference_in_RD =
      Yes - No
  )


###########################################################################################
# Run for Poisson
###########################################################################################

poisson_pfas_adt <- function(
    data,
    pfas_var,
    outcome_var,
    covars,
    rate_scale = 1000
) {
  
  #------------------------------------------------------------
  # 1. Restrict analytic population
  #------------------------------------------------------------
  
  dat <- data %>%
    filter(
      !is.na(txhorm.f),
      txhorm.f %in% c("No", "Yes"),
      stage.f != "Distant",
      !(stage.f == "Localized" &
          grade_combo.f == "<7" &
          psa.f == "<10")
    ) %>%
    mutate(
      interval_py = (tstop - tstart) / 12,
      adt = ifelse(txhorm.f == "Yes", 1, 0)
    ) %>%
    filter(
      interval_py > 0,
      !is.na(.data[[pfas_var]]),
      !is.na(.data[[outcome_var]])
    )
  
  
  #------------------------------------------------------------
  # 2. Patient-level median PFAS
  #
  # Important: calculate median using one observation/person
  # so patients with more interval rows don't receive more weight.
  #------------------------------------------------------------
  
  pfas_median <- dat %>%
    distinct(patient_study_id, .data[[pfas_var]]) %>%
    summarise(
      med = median(.data[[pfas_var]], na.rm = TRUE)
    ) %>%
    pull(med)
  
  
  dat <- dat %>%
    mutate(
      pfas_high = ifelse(
        .data[[pfas_var]] >= pfas_median,
        1, 0
      )
    )
  
  
  #------------------------------------------------------------
  # 3. Construct model formula
  #------------------------------------------------------------
  
  rhs <- c(
    "pfas_high * adt",
    covars
  )
  
  form <- as.formula(
    paste(
      outcome_var,
      "~",
      paste(rhs, collapse = " + ")
    )
  )
  
  
  #------------------------------------------------------------
  # 4. Fit Poisson model
  #------------------------------------------------------------
  
  fit <- glm(
    form,
    family = poisson(link = "log"),
    offset = log(interval_py),
    data = dat,
    na.action = na.omit
  )
  
  
  #------------------------------------------------------------
  # 5. Identify exact rows actually used by model
  #------------------------------------------------------------
  
  used_rows <- as.numeric(rownames(model.frame(fit)))
  
  model_dat <- dat[used_rows, ]
  
  
  #------------------------------------------------------------
  # 6. Cluster-robust covariance matrix
  #------------------------------------------------------------
  
  robust_vcov <- sandwich::vcovCL(
    fit,
    cluster = model_dat$patient_study_id,
    type = "HC0"
  )
  
  
  #------------------------------------------------------------
  # 7. Standardization population
  #
  # One row per person avoids giving more weight to people who
  # happen to have more counting-process intervals.
  #
  # Assumes covariates used below are baseline/person-level.
  #------------------------------------------------------------
  
  std_dat <- model_dat %>%
    arrange(patient_study_id, tstart) %>%
    distinct(patient_study_id, .keep_all = TRUE)
  
  
  beta <- coef(fit)
  
  
  #------------------------------------------------------------
  # Helper: standardized rate + delta-method gradient
  #------------------------------------------------------------
  
  get_rate <- function(pfas_value, adt_value) {
    
    nd <- std_dat %>%
      mutate(
        pfas_high = pfas_value,
        adt = adt_value
      )
    
    X <- model.matrix(
      delete.response(terms(fit)),
      data = nd
    )
    
    X <- X[, names(beta), drop = FALSE]
    
    lp <- drop(X %*% beta)
    
    individual_rates <- exp(lp)
    
    rate <- mean(individual_rates) * rate_scale
    
    # Gradient for delta-method variance
    grad <- colMeans(
      X * individual_rates
    ) * rate_scale
    
    list(
      estimate = rate,
      gradient = grad
    )
  }
  
  
  #------------------------------------------------------------
  # 8. Four standardized PFAS x ADT rates
  #------------------------------------------------------------
  
  low_no   <- get_rate(0, 0)
  high_no  <- get_rate(1, 0)
  low_yes  <- get_rate(0, 1)
  high_yes <- get_rate(1, 1)
  
  
  #------------------------------------------------------------
  # Helper for CI from gradient
  #------------------------------------------------------------
  
  make_estimate <- function(est, grad) {
    
    var <- drop(
      t(grad) %*% robust_vcov %*% grad
    )
    
    se <- sqrt(var)
    
    tibble(
      estimate = est,
      se = se,
      lower95 = est - 1.96 * se,
      upper95 = est + 1.96 * se
    )
  }
  
  
  #------------------------------------------------------------
  # 9. Rates
  #------------------------------------------------------------
  
  rates <- bind_rows(
    
    make_estimate(
      low_no$estimate,
      low_no$gradient
    ) %>%
      mutate(
        PFAS = "Low",
        ADT = "No"
      ),
    
    make_estimate(
      high_no$estimate,
      high_no$gradient
    ) %>%
      mutate(
        PFAS = "High",
        ADT = "No"
      ),
    
    make_estimate(
      low_yes$estimate,
      low_yes$gradient
    ) %>%
      mutate(
        PFAS = "Low",
        ADT = "Yes"
      ),
    
    make_estimate(
      high_yes$estimate,
      high_yes$gradient
    ) %>%
      mutate(
        PFAS = "High",
        ADT = "Yes"
      )
    
  ) %>%
    mutate(
      pfas = pfas_var,
      outcome = outcome_var,
      measure = paste0(
        "Adjusted rate per ",
        rate_scale,
        " PY"
      )
    ) %>%
    select(
      pfas,
      outcome,
      PFAS,
      ADT,
      measure,
      estimate,
      se,
      lower95,
      upper95
    )
  
  
  #------------------------------------------------------------
  # 10. PFAS rate difference among NO ADT
  #
  # High PFAS - Low PFAS
  #------------------------------------------------------------
  
  rd_no_est <-
    high_no$estimate -
    low_no$estimate
  
  rd_no_grad <-
    high_no$gradient -
    low_no$gradient
  
  
  #------------------------------------------------------------
  # 11. PFAS rate difference among YES ADT
  #------------------------------------------------------------
  
  rd_yes_est <-
    high_yes$estimate -
    low_yes$estimate
  
  rd_yes_grad <-
    high_yes$gradient -
    low_yes$gradient
  
  
  #------------------------------------------------------------
  # 12. Difference in rate differences
  #
  # (High - Low PFAS | ADT yes)
  #       -
  # (High - Low PFAS | ADT no)
  #------------------------------------------------------------
  
  did_est <-
    rd_yes_est -
    rd_no_est
  
  did_grad <-
    rd_yes_grad -
    rd_no_grad
  
  
  contrasts <- bind_rows(
    
    make_estimate(
      rd_no_est,
      rd_no_grad
    ) %>%
      mutate(
        contrast = "PFAS RD among No ADT"
      ),
    
    make_estimate(
      rd_yes_est,
      rd_yes_grad
    ) %>%
      mutate(
        contrast = "PFAS RD among Yes ADT"
      ),
    
    make_estimate(
      did_est,
      did_grad
    ) %>%
      mutate(
        contrast = "Difference in rate differences"
      )
    
  ) %>%
    mutate(
      pfas = pfas_var,
      outcome = outcome_var,
      units = paste0(
        "deaths per ",
        rate_scale,
        " PY"
      )
    ) %>%
    select(
      pfas,
      outcome,
      contrast,
      units,
      estimate,
      se,
      lower95,
      upper95
    )
  
  
  #------------------------------------------------------------
  # Return model + results
  #------------------------------------------------------------
  
  list(
    model = fit,
    robust_vcov = robust_vcov,
    pfas_median = pfas_median,
    rates = rates,
    contrasts = contrasts
  )
}


covars <- c(
  "age_cat.f",
  "dxyr_cat.f",
  "nhrace",
  "married.f",
  "insur.f",
  "stage.f",
  "grade_combo.f",
  "txsurg.f",
  "txrad.f",
  "psa.f",
  "nses_quin",
  "no2_new",
  "pm25_new",
  "ndvi27_final"
)

pfas_vars <- c(
  "pfoa_dx",
  "pfos_dx",
  "pfna_dx"
)

outcome_vars <- c(
  "all_death",
  "cvd_death",
  "prostate_death"
)

model_grid <- tidyr::crossing(
  pfas = pfas_vars,
  outcome = outcome_vars
)

results <- model_grid %>%
  mutate(
    fit = map2(
      pfas,
      outcome,
      ~ poisson_pfas_adt(
        data = long_imp_list[[1]],
        pfas_var = .x,
        outcome_var = .y,
        covars = covars
      )
    )
  )


adjusted_rates <- results %>%
  mutate(
    rates = map(fit, "rates")
  ) %>%
  select(rates) %>%
  unnest(rates)

adjusted_rates

rate_differences <- results %>%
  mutate(
    contrasts = map(fit, "contrasts")
  ) %>%
  select(contrasts) %>%
  unnest(contrasts)

rate_differences

#######################################################################
# Multiplicative Cox
#######################################################################


## keep nonmetastatic, non-low risk
analysis_dat <- long_imp_list[[1]] %>%
  filter(
    stage.f != "Distant",
    !(stage.f == "Localized" &
        grade_combo.f == "<7" &
        psa.f == "<10")
  )  %>%
  mutate(pfna_cat = factor(case_when(pfna < 2 ~ 1,
                                     pfna >=2 & pfna < 3.8 ~ 2, 
                                     pfna >= 3.8 ~ 3)) )

## CVD mortality

coxph(
  Surv(tstart, tstop, cvd_death) ~
    pfoa_quin * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
  dxyr_cat.f +
  nhrace +
  married.f +
  insur.f +
  nses_quin +
  no2_new +
  pm25_new +
  ndvi27_final,
  data = analysis_dat
) 


coxph(
  Surv(tstart, tstop, cvd_death) ~
    pfos_quin * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
    dxyr_cat.f +
    nhrace +
    married.f +
    insur.f +
    nses_quin +
    no2_new +
    pm25_new +
    ndvi27_final,
  data = analysis_dat
) 


coxph(
  Surv(tstart, tstop, cvd_death) ~
    pfna_cat * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
    dxyr_cat.f +
    nhrace +
    married.f +
    insur.f +
    nses_quin +
    no2_new +
    pm25_new +
    ndvi27_final,
  data = analysis_dat
) 


length(unique(analysis_dat$patient_study_id))

## Prostate mortality

coxph(
  Surv(tstart, tstop, prostate_death) ~
    pfoa_quin * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
    dxyr_cat.f +
    nhrace +
    married.f +
    insur.f +
    nses_quin +
    no2_new +
    pm25_new +
    ndvi27_final,
  data = analysis_dat
) 


coxph(
  Surv(tstart, tstop, prostate_death) ~
    pfos_quin * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
    dxyr_cat.f +
    nhrace +
    married.f +
    insur.f +
    nses_quin +
    no2_new +
    pm25_new +
    ndvi27_final,
  data = analysis_dat
) 


coxph(
  Surv(tstart, tstop, prostate_death) ~
    pfna_cat * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
    dxyr_cat.f +
    nhrace +
    married.f +
    insur.f +
    nses_quin +
    no2_new +
    pm25_new +
    ndvi27_final,
  data = analysis_dat
) 


## all-cause mortality


coxph(
  Surv(tstart, tstop, all_death) ~
    pfoa_quin * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
    dxyr_cat.f +
    nhrace +
    married.f +
    insur.f +
    nses_quin +
    no2_new +
    pm25_new +
    ndvi27_final,
  data = analysis_dat
) 


coxph(
  Surv(tstart, tstop, all_death) ~
    pfos_quin * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
    dxyr_cat.f +
    nhrace +
    married.f +
    insur.f +
    nses_quin +
    no2_new +
    pm25_new +
    ndvi27_final,
  data = analysis_dat
) 


coxph(
  Surv(tstart, tstop, all_death) ~
    pfna_cat * txhorm.f +
    grade_combo.f +
    psa.f +
    stage.f +
    age_cat.f +
    dxyr_cat.f +
    nhrace +
    married.f +
    insur.f +
    nses_quin +
    no2_new +
    pm25_new +
    ndvi27_final,
  data = analysis_dat
)

