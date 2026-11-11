#################################################################
# R Code - multiple imputation

# Megan Shanahan, Hari Iyer
# Date: 11/2/2026

# Goal: Use multiple imputation to resolve missingness
#################################################################


library(tidyverse)
library(mice)
library(survival)
library(forcats)
library(broom)
library(dplyr)
library(ggplot2)




############################# Cox models ###############################

# load("~/R/registry/long_list_imp.Rda")
load("/extra/home/hi97/R/registry/nj_res_hx/long_imp_list_lagged.Rda")
names(long_imp_list_lagged[[1]])

# defining exposures and outcomes 
pfas_quint <- c("pfoa_lag_quin", "pfos_lag_quin", "pfna_lag_quin")
pfas_trend <- c("pfoa_exp_lag_t", "pfos_exp_lag_t", "pfna_exp_lag_t")
pfas_iqr <- c("pfoa_exp_lag_per_iqr", "pfos_exp_lag_per_iqr", "pfna_exp_lag_per_iqr")
pfas_mcl   <- c("pfoa_mcl", "pfos_mcl", "pfna_mcl")

outcomes <- c(
  all_cause = "all_death",
  cvd       = "cvd_death",
  prostate  = "prostate_death"
)
run_pfas_cox <- function(exposure, outcome, long_imp_list_lagged, exposure_type) {
  
  fit_list <- lapply(long_imp_list_lagged, function(dat) {
    
    form <- as.formula(paste0(
      "Surv(tstart, tstop, ", outcome, ") ~ ",
      exposure,
      " + age_cat.f + dxyr_cat.f + ",
      "nhrace + married.f + insur.f + stage.f + txsurg.f + txrad.f + ",
      "psa.f + as.factor(nses_lag_quin) + no2_new_lag + pm25_new_lag + ndvi27_final_lag"
    ))
    
    coxph(form, data = dat)
  })
  
  tab <- summary(pool(as.mira(fit_list)))
  
  tab <- tab %>%
    dplyr::filter(grepl(exposure, term)) %>%
    dplyr::mutate(
      exposure = exposure,
      outcome = outcome,
      exposure_type = exposure_type,
      HR = exp(estimate),
      HR_95LL = exp(estimate - 1.96 * std.error),
      HR_95UL = exp(estimate + 1.96 * std.error)
    )
  
  if (exposure_type == "quintile") {
    tab <- tab %>%
      dplyr::mutate(
        contrast = dplyr::case_when(
          row_number() == 1 ~ "Q2 vs Q1",
          row_number() == 2 ~ "Q3 vs Q1",
          row_number() == 3 ~ "Q4 vs Q1",
          row_number() == 4 ~ "Q5 vs Q1",
          TRUE ~ as.character(term)
        )
      )
  } else if (exposure_type == "trend") {
    tab <- tab %>%
      dplyr::mutate(
        contrast = "trend"
      )
  } else if (exposure_type == "iqr") {
    tab <- tab %>%
      dplyr::mutate(
        contrast = "per IQR increase"
      )
  } else if (exposure_type == "iqr") {
    tab <- tab %>%
      dplyr::mutate(
        contrast = "per IQR increase"
      )
  } else if (exposure_type == "mcl") {
    
    tab <- tab |>
      mutate(contrast = "At or Above MCL vs Below MCL")
  }
  
  tab <- tab %>%
    dplyr::mutate(contrast = as.character(contrast))
  
  tab
}

# creating mcl variables


long_imp_list_lagged <- lapply(long_imp_list_lagged, function(dat) {
  
  dat <- dat |>
    mutate(
      pfoa_mcl = factor(
        if_else(
          is.na(pfoa_dx),
          NA_character_,
          if_else(pfoa_exp_lag >= 14, "At or Above MCL", "Below MCL")
        ),
        levels = c("Below MCL", "At or Above MCL")
      ),
      
      pfos_mcl = factor(
        if_else(
          is.na(pfos_dx),
          NA_character_,
          if_else(pfos_exp_lag >= 13, "At or Above MCL", "Below MCL")
        ),
        levels = c("Below MCL", "At or Above MCL")
      ),
      
      pfna_mcl = factor(
        if_else(
          is.na(pfna_dx),
          NA_character_,
          if_else(pfna_exp_lag >= 13, "At or Above MCL", "Below MCL")
        ),
        levels = c("Below MCL", "At or Above MCL")
      )
    )
  
  dat
})

long_imp_list_lagged[[1]] |>
  distinct(patient_study_id, pfoa_mcl, pfos_mcl, pfna_mcl) |>
  summarise(
    pfoa_below = sum(pfoa_mcl == "Below MCL", na.rm = TRUE),
    pfoa_above = sum(pfoa_mcl == "At or Above MCL", na.rm = TRUE),
    pfos_below = sum(pfos_mcl == "Below MCL", na.rm = TRUE),
    pfos_above = sum(pfos_mcl == "At or Above MCL", na.rm = TRUE),
    pfna_below = sum(pfna_mcl == "Below MCL", na.rm = TRUE),
    pfna_above = sum(pfna_mcl == "At or Above MCL", na.rm = TRUE)
  )

# seeing how many patients ever exceed the MCL
long_imp_list_lagged[[1]] |>
  group_by(patient_study_id) |>
  summarise(
    pfoa_above = any(pfoa_mcl == "At or Above MCL", na.rm = TRUE),
    pfos_above = any(pfos_mcl == "At or Above MCL", na.rm = TRUE),
    pfna_above = any(pfna_mcl == "At or Above MCL", na.rm = TRUE)
  ) |>
  summarise(
    pfoa_above = sum(pfoa_above),
    pfos_above = sum(pfos_above),
    pfna_above = sum(pfna_above)
  )


################################################# CASE COUNTS AND PERSON YEARS ###########################
load("/extra/home/hi97/R/registry/nj_res_hx/long_imp_list_lagged.Rda")

dat_counts <- long_imp_list_lagged[[1]]

get_quint_counts <- function(dat, exposure, outcome) {
  
  dat %>%
    group_by(.data[[exposure]]) %>%
    summarise(
      cases = sum(.data[[outcome]], na.rm = TRUE),
      person_years = sum((tstop - tstart) / 365.25, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(quintile_value = all_of(exposure)) %>%
    mutate(
      exposure = exposure,
      outcome = outcome,
      quintile = paste0("Q", as.character(quintile_value))
    )
}

quint_counts <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    do.call(
      rbind,
      lapply(pfas_quint, function(exposure) {
        
        get_quint_counts(
          dat = dat_counts,
          exposure = exposure,
          outcome = outcomes[[outcome_name]]
        ) %>%
          mutate(outcome_label = outcome_name)
        
      })
    )
  })
)

write.csv(
  quint_counts,
  "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/output/quint_counts_lagged.csv",
  row.names = FALSE
)


overall_counts <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    outcome_var <- outcomes[[outcome_name]]
    
    data.frame(
      outcome_label = outcome_name,
      outcome = outcome_var,
      cases = sum(dat_counts[[outcome_var]], na.rm = TRUE),
      person_years = sum((dat_counts$tstop - dat_counts$tstart) / 365.25, na.rm = TRUE)
    )
  })
)

write.csv(
  overall_counts,
  "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/output/overall_counts_lagged.csv",
  row.names = FALSE
)

# by pfas
pfas <- c("pfoa", "pfos", "pfna")

get_continuous_pfas_counts <- function(dat, pfas_var, outcome_var,
                                       pfas_label, outcome_label) {
  
  dat %>%
    filter(
      !is.na(.data[[pfas_var]]),
      !is.na(.data[[outcome_var]])
    ) %>%
    summarise(
      pfas = pfas_label,
      outcome_label = outcome_label,
      outcome = outcome_var,
      cases = sum(.data[[outcome_var]]),
      person_years = sum(
        (tstop - tstart) / 365.25,
        na.rm = TRUE
      )
    )
}

pfas_counts <- do.call(
  rbind,
  lapply(pfas, function(pfas_var) {
    
    do.call(
      rbind,
      lapply(names(outcomes), function(outcome_label) {
        
        get_continuous_pfas_counts(
          dat = dat_counts,
          pfas_var = pfas_var,
          outcome_var = outcomes[[outcome_label]],
          pfas_label = pfas_var,
          outcome_label = outcome_label
        )
      })
    )
  })
)

write.csv(
  pfas_counts,
  "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/output/lagged_pfas_continuous_counts.csv",
  row.names = FALSE
)

############### running models #####################################

results_quint <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_quint, function(exposure) {
    run_pfas_cox(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list_lagged = long_imp_list_lagged,
      exposure_type = "quintile"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))

results_trend <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_trend, function(exposure) {
    run_pfas_cox(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list_lagged = long_imp_list_lagged,
      exposure_type = "trend"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))

results_iqr <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_iqr, function(exposure) {
    run_pfas_cox(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list_lagged = long_imp_list_lagged,
      exposure_type = "iqr"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))


results_iqr <- results_iqr |>
  dplyr::select(
    outcome_label,
    outcome,
    exposure,
    exposure_type,
    contrast,
    HR,
    HR_95LL,
    HR_95UL,
    p.value
  )
  



summary(long_imp_list_lagged[[1]][, c("pfoa_dx", "pfos_dx", "pfna_dx")])

final_results <- bind_rows(results_quint, results_trend, results_iqr) %>%
  dplyr::select(
    outcome_label,
    outcome,
    exposure,
    exposure_type,
    contrast,
    HR,
    HR_95LL,
    HR_95UL,
    p.value
  )


# adding a footnote with covariates used for adjusting

# Create an empty row with exactly the same columns as df
footnote <- final_results[1, ]
footnote[] <- NA

# Put the note in the first column
footnote[[1]] <- "Models are adjusted for age, diagnosis year, race/ethnicity, marital status, insurance status, stage, treatment (surgery and or radiation), PSA,  nSES quintile, no2, pm2.5, and ndvi27"

# Add it to the bottom
final_results <- rbind(final_results, footnote)


write.csv(
  final_results,
  "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/output/lag_full_pfas_mortality_results.csv",
  row.names = FALSE
)


## age and year of diagnosis only 
run_pfas_cox_age_dxyr <- function(exposure, outcome, long_imp_list_lagged, exposure_type) {
  
  fit_list <- lapply(long_imp_list_lagged, function(dat) {
    
    form <- as.formula(paste0(
      "Surv(tstart, tstop, ", outcome, ") ~ ",
      exposure,
      " + age_cat.f + dxyr_cat.f  "
    ))
    
    coxph(form, data = dat)
  })
  
  tab <- summary(pool(as.mira(fit_list)))
  
  
  tab <- tab %>%
    dplyr::filter(grepl(exposure, term)) %>%
    dplyr::mutate(
      exposure = exposure,
      outcome = outcome,
      exposure_type = exposure_type,
      HR = exp(estimate),
      HR_95LL = exp(estimate - 1.96 * std.error),
      HR_95UL = exp(estimate + 1.96 * std.error)
    )
  
  if (exposure_type == "quintile") {
    tab <- tab %>%
      dplyr::mutate(
        contrast = dplyr::case_when(
          row_number() == 1 ~ "Q2 vs Q1",
          row_number() == 2 ~ "Q3 vs Q1",
          row_number() == 3 ~ "Q4 vs Q1",
          row_number() == 4 ~ "Q5 vs Q1",
          TRUE ~ as.character(term)
        )
      )
  } else if (exposure_type == "trend") {
    tab <- tab %>%
      dplyr::mutate(
        contrast = "trend"
      )
  } else if (exposure_type == "iqr") {
    tab <- tab %>%
      dplyr::mutate(
        contrast = "per IQR increase"
      )
  } else if (exposure_type == "mcl") {
    
    tab <- tab |>
      dplyr::mutate(
        contrast = "At or Above MCL vs below MCL"
      )
    
  } else {
    
    tab <- tab |>
      dplyr::mutate(
        contrast = as.character(term)
      )
  }
  tab <- tab %>%
    dplyr::mutate(contrast = as.character(contrast))
  
  tab
}


results_age_dxyr_quint <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_quint, function(exposure) {
    run_pfas_cox_age_dxyr(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list_lagged = long_imp_list_lagged,
      exposure_type = "quintile"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))

results_age_dxyr_trend <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_trend, function(exposure) {
    run_pfas_cox_age_dxyr(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list_lagged = long_imp_list_lagged,
      exposure_type = "trend"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))

results_iqr_age_dxyr <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_iqr, function(exposure) {
    run_pfas_cox_age_dxyr(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list_lagged = long_imp_list_lagged,
      exposure_type = "iqr"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))





final_results_age_dxyr <- bind_rows(
  results_age_dxyr_quint,
  results_age_dxyr_trend,
  results_iqr_age_dxyr,
) %>%
  dplyr::select(
    outcome_label,
    outcome,
    exposure,
    exposure_type,
    contrast,
    HR,
    HR_95LL,
    HR_95UL,
    p.value
  )



# adding a footnote with covariates used for adjusting

# Create an empty row with exactly the same columns as df
footnote <- final_results_age_dxyr[1, ]
footnote[] <- NA

# Put the note in the first column
footnote[[1]] <- "Models are adjusted for age and diagnosis year"

# Add it to the bottom
final_results_age_dxyr <- rbind(final_results_age_dxyr, footnote)


write.csv(
  final_results_age_dxyr,
  "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/output/lag_pfas_mortality_age_dxyr.csv",
  row.names = FALSE
)