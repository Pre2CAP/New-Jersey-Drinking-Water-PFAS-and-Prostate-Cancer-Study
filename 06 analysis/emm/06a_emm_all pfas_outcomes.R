#################################################################
# R Code - EMM

# Hari Iyer
# Date: 09/23/2026

# Goal: Effect measure modification between pfoa and cvd for treatment
#################################################################


library(tidyverse)
library(mice)
library(survival)
library(forcats)
library(broom)
library(dplyr)
library(ggplot2)

load("/extra/home/hi97/R/registry/nj_res_hx/long_list_imp.Rda")


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

exposures <- c(
  PFOA = "pfoa_iqr",
  PFOS = "pfos_iqr",
  PFNA = "pfna_iqr"
)

outcomes <- c(
  All_cause = "all_death",
  Prostate  = "prostate_death",
  CVD       = "cvd_death"
)

effect_modifiers <- c(
  race  = "nhrace",
  age   = "dxage65",
  stage = "stage.f"
)


#####################################################################
# Function
#####################################################################

run_pfas_interaction <- function(exposure,
                                 outcome,
                                 effect_modifier,
                                 modifier_name,
                                 long_imp_list) {
  
  message( "Exposure: ", exposure,
           " | Outcome: ", outcome,
           " | Effect modifier: ", modifier_name)
  
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
  
  
  # Count only terms involving the exposure × modifier interaction
  coef_names <- names(coef(fit_int_mira$analyses[[1]]))
  
  interaction_terms <- coef_names[
    grepl(":", coef_names) &
      grepl(exposure, coef_names) &
      grepl(effect_modifier, coef_names)
  ]
  
  n_int_terms <- length(interaction_terms)
  
  
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
      exposure = exposure,
      outcome = outcome,
      effect_modifier = modifier_name,
      p_interaction = p_interaction
    )
  
  
  # ---------------------------------------------------------------
  # Return everything
  # ---------------------------------------------------------------
  
  list(
    exposure = exposure,
    outcome = outcome,
    modifier = modifier_name,
    base_model = fit_base_mira,
    interaction_model = fit_int_mira,
    pooled_model = pooled,
    p_interaction = p_interaction,
    df_interaction = n_int_terms
  )
}



pool_contrast_cox <- function(fit_mira, terms) {
  
  models <- fit_mira$analyses
  
  get_contrast <- function(model) {
    
    beta <- coef(model)
    V <- vcov(model)
    
    L <- rep(0, length(beta))
    names(L) <- names(beta)
    
    if (!all(terms %in% names(beta))) {
      stop(
        "These terms were not found in the model: ",
        paste(terms[!terms %in% names(beta)], collapse = ", ")
      )
    }
    
    L[terms] <- 1
    
    estimate <- sum(L * beta)
    
    variance <- as.numeric(
      t(L) %*% V %*% L
    )
    
    data.frame(
      estimate = estimate,
      variance = variance
    )
  }
  
  contrast_list <- lapply(models, get_contrast)
  
  estimates <- sapply(
    contrast_list,
    function(x) x$estimate
  )
  
  variances <- sapply(
    contrast_list,
    function(x) x$variance
  )
  
  pooled <- mice::pool.scalar(
    Q = estimates,
    U = variances
  )
  
  estimate <- pooled$qbar
  se <- sqrt(pooled$t)
  
  HR <- exp(estimate)
  HR_95LL <- exp(estimate - 1.96 * se)
  HR_95UL <- exp(estimate + 1.96 * se)
  
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

extract_subgroup_results <- function(result,
                                     exposure_label,
                                     outcome_label) {
  
  fit <- result$interaction_model
  exposure <- result$exposure
  modifier <- result$modifier
  
  coef_names <- names(
    coef(fit$analyses[[1]])
  )
  
  
  # ---------------------------------------------------------------
  # RACE
  # ---------------------------------------------------------------
  
  if (modifier == "race") {
    
    groups <- c(
      NHW = NA,
      NHB = "NHB",
      NHAAPI = "NHAAPI",
      Hispanic = "Hispanic"
    )
    
    results <- lapply(
      names(groups),
      function(group) {
        
        # Reference group
        if (group == "NHW") {
          
          terms <- exposure
          
        } else {
          
          # Find interaction term robustly
          int_term <- coef_names[
            grepl(":", coef_names) &
              grepl(exposure, coef_names) &
              grepl(
                paste0("nhrace", groups[[group]]),
                coef_names,
                fixed = TRUE
              )
          ]
          
          terms <- c(
            exposure,
            int_term
          )
        }
        
        pool_contrast_cox(
          fit,
          terms = terms
        ) %>%
          mutate(subgroup = group)
      }
    )
    
    results <- bind_rows(results)
  }
  
  
  # ---------------------------------------------------------------
  # AGE
  # ---------------------------------------------------------------
  
  if (modifier == "age") {
    
    age_ref <- pool_contrast_cox(
      fit,
      terms = exposure
    ) %>%
      mutate(subgroup = "<65")
    
    
    int_term <- coef_names[
      grepl(":", coef_names) &
        grepl(exposure, coef_names) &
        grepl("dxage6565+", coef_names, fixed = TRUE)
    ]
    
    
    age_65 <- pool_contrast_cox(
      fit,
      terms = c(
        exposure,
        int_term
      )
    ) %>%
      mutate(subgroup = "65+")
    
    
    results <- bind_rows(
      age_ref,
      age_65
    )
  }
  
  
  # ---------------------------------------------------------------
  # STAGE
  # ---------------------------------------------------------------
  
  if (modifier == "stage") {
    
    stage_ref <- pool_contrast_cox(
      fit,
      terms = exposure
    ) %>%
      mutate(subgroup = "Localized")
    
    
    regional_term <- coef_names[
      grepl(":", coef_names) &
        grepl(exposure, coef_names) &
        grepl("stage.fRegional", coef_names, fixed = TRUE)
    ]
    
    
    distant_term <- coef_names[
      grepl(":", coef_names) &
        grepl(exposure, coef_names) &
        grepl("stage.fDistant", coef_names, fixed = TRUE)
    ]
    
    
    stage_regional <- pool_contrast_cox(
      fit,
      terms = c(
        exposure,
        regional_term
      )
    ) %>%
      mutate(subgroup = "Regional")
    
    
    stage_distant <- pool_contrast_cox(
      fit,
      terms = c(
        exposure,
        distant_term
      )
    ) %>%
      mutate(subgroup = "Distant")
    
    
    results <- bind_rows(
      stage_ref,
      stage_regional,
      stage_distant
    )
  }
  
  
  # ---------------------------------------------------------------
  # Add identifiers
  # ---------------------------------------------------------------
  
  results %>%
    mutate(
      exposure = exposure_label,
      outcome = outcome_label,
      effect_modifier = modifier,
      p_interaction = result$p_interaction
    ) %>%
    select(
      exposure,
      outcome,
      effect_modifier,
      subgroup,
      HR,
      HR_95LL,
      HR_95UL,
      p.value,
      p_interaction,
      everything()
    )
}

###################################################################
# Analysis grid
##################################################################

analysis_grid <- expand.grid(
  exposure_name = names(exposures),
  outcome_name = names(outcomes),
  modifier_name = names(effect_modifiers),
  stringsAsFactors = FALSE
)

#####################################################################
# Run one analysis at a time and retain only summarized results
#####################################################################

results_list <- lapply(
  seq_len(nrow(analysis_grid)),
  function(i) {
    
    exposure_name <- analysis_grid$exposure_name[i]
    outcome_name  <- analysis_grid$outcome_name[i]
    modifier_name <- analysis_grid$modifier_name[i]
    
    # -------------------------------------------------------------
    # Fit one analysis
    # -------------------------------------------------------------
    
    result <- run_pfas_interaction(
      exposure = exposures[[exposure_name]],
      outcome = outcomes[[outcome_name]],
      effect_modifier = effect_modifiers[[modifier_name]],
      modifier_name = modifier_name,
      long_imp_list = long_imp_list
    )
    
    
    # -------------------------------------------------------------
    # Immediately calculate subgroup-specific HRs
    # pool_contrast_cox() gets called inside this function
    # -------------------------------------------------------------
    
    subgroup_results <- extract_subgroup_results(
      result = result,
      exposure_label = exposure_name,
      outcome_label = outcome_name
    )
    
    
    # -------------------------------------------------------------
    # Also retain the overall interaction test
    # -------------------------------------------------------------
    
    p_int_result <- data.frame(
      exposure = exposure_name,
      outcome = outcome_name,
      effect_modifier = modifier_name,
      df = result$df_interaction,
      p_interaction = result$p_interaction
    )
    
    
    # -------------------------------------------------------------
    # Return ONLY small result tables
    # -------------------------------------------------------------
    
    list(
      subgroup = subgroup_results,
      p_interaction = p_int_result
    )
  }
)

#####################################################################
# Combine subgroup results
#####################################################################

all_subgroup_results <- bind_rows(
  lapply(results_list, function(x) x$subgroup)
)

p_interaction_results <- bind_rows(
  lapply(results_list, function(x) x$p_interaction)
)

## Inspect

all_subgroup_results
p_interaction_results


#####################################################################
# Combine interaction P-values
#####################################################################



write.csv(
  all_subgroup_results,
  "~/pfas_mortality_effect_modification.csv",
  row.names = FALSE
)

write.csv(
  p_interaction_results,
  "~/pfas_mortality_p_interactions.csv",
  row.names = FALSE
)

