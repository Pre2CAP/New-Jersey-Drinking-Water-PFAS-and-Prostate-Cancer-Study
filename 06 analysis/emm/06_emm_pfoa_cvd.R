#################################################################
# R Code - EMM

# Megan Shanahan
# Date: 06/02/2026

# Goal: Effect measure modification between pfoa and cvd for age, race, stage
#################################################################


library(tidyverse)
library(mice)
library(survival)
library(forcats)
library(broom)
library(dplyr)
library(ggplot2)

load("~/R/registry/long_list_imp.Rda")

levels(long_imp_list[[1]]$nhrace)

levels(long_imp_list[[1]]$age_cat.f)

levels(long_imp_list[[1]]$stage.f)

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
  race  = "nhrace",
  age   = "dxage65",
  stage = "stage.f"
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
  file = "~/R/registry/pfoa_cvd_interactions.Rda"
)

load("~/R/registry/pfoa_cvd_interactions.Rda")

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
# RACE / ETHNICITY
#####################################################################

fit_race <- pfoa_cvd_interactions$race$interaction_model

race_nhw <- pool_contrast_cox(
  fit_race,
  terms = c("pfoa_iqr")
)

race_nhb <- pool_contrast_cox(
  fit_race,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:nhraceNHB"
  )
)

race_nhaapi <- pool_contrast_cox(
  fit_race,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:nhraceNHAAPI"
  )
)

race_hispanic <- pool_contrast_cox(
  fit_race,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:nhraceHispanic"
  )
)

race_results <- bind_rows(
  NHW      = race_nhw,
  NHB      = race_nhb,
  NHAAPI   = race_nhaapi,
  Hispanic = race_hispanic,
  .id = "subgroup"
) |>
  mutate(
    effect_modifier = "Race/ethnicity",
    p_interaction = pfoa_cvd_interactions$race$p_interaction
  )

race_results
#####################################################################
# AGE AT DIAGNOSIS: <65 vs 65+
#####################################################################

fit_age <- pfoa_cvd_interactions$age$interaction_model

# <65 years - reference group
age_lt65 <- pool_contrast_cox(
  fit_age,
  terms = c(
    "pfoa_iqr"
  )
)

# 65+ years
age_65plus <- pool_contrast_cox(
  fit_age,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:dxage6565+"
  )
)

age_results <- bind_rows(
  "<65" = age_lt65,
  "65+" = age_65plus,
  .id = "subgroup"
) |>
  mutate(
    effect_modifier = "Age at diagnosis",
    p_interaction =
      pfoa_cvd_interactions$age$p_interaction
  )

age_results

# save(
#  age_results,
#  file = "~/age_results.Rda"
# )

#####################################################################
# STAGE
#####################################################################

fit_stage <- pfoa_cvd_interactions$stage$interaction_model

stage_localized <- pool_contrast_cox(
  fit_stage,
  terms = c("pfoa_iqr")
)

stage_regional <- pool_contrast_cox(
  fit_stage,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:stage.fRegional"
  )
)

stage_distant <- pool_contrast_cox(
  fit_stage,
  terms = c(
    "pfoa_iqr",
    "pfoa_iqr:stage.fDistant"
  )
)

stage_results <- bind_rows(
  Localized = stage_localized,
  Regional  = stage_regional,
  Distant   = stage_distant,
  .id = "subgroup"
) |>
  mutate(
    effect_modifier = "Stage",
    p_interaction = pfoa_cvd_interactions$stage$p_interaction
  )


#####################################################################
# COMBINE RESULTS
#####################################################################

pfoa_cvd_effect_modification <- bind_rows(
  race_results,
  age_results,
  stage_results
)

pfoa_cvd_effect_modification

save(
  pfoa_cvd_interactions,
  file = "~/pfoa_cvd_interactions.Rda"
)

write.csv(
  pfoa_cvd_effect_modification,
  "~/pfoa_cvd_effect_modification.csv",
  row.names = FALSE
)




