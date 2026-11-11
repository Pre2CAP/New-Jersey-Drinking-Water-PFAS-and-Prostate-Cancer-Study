#################################################################
# R Code Cox models

# Megan Shanahan
# Date: 09/17/2026

# Goal: Cox models - PFNA
#################################################################


library(tidyverse)
library(mice)
library(survival)
library(forcats)
library(broom)
library(dplyr)
library(ggplot2)

# load data
load("/extra/home/hi97/R/registry/nj_res_hx/long_list_imp.Rda")



# Creating new categories for PFNA
long_imp_list <- lapply(
  long_imp_list,
  function(dat) {
    
    dat %>%
      mutate(
        pfna_cat = factor(
          case_when(
            pfna < 2 ~ 1,
            pfna >=2 & pfna < 3.8 ~ 2, 
            pfna >= 3.8 ~ 3
            )
          ) 
        ) %>%
      group_by(pfna_cat) %>%
      mutate(
        pfna_dx_t = median(pfna, na.rm = TRUE)
      ) %>%
      ungroup()
  }
)



# defining exposures and outcomes 
pfas_quint <- c("pfna_cat")
pfas_trend <- c( "pfna_dx_t")
pfas_iqr <- c("pfna_iqr")
pfas_mcl   <- c("pfna_mcl")

outcomes <- c(
  all_cause = "all_death",
  cvd       = "cvd_death",
  prostate  = "prostate_death"
)
run_pfas_cox <- function(exposure, outcome, long_imp_list, exposure_type) {
  
  fit_list <- lapply(long_imp_list, function(dat) {
    
    form <- as.formula(paste0(
      "Surv(tstart, tstop, ", outcome, ") ~ ",
      exposure,
      " + age_cat.f + dxyr_cat.f + ",
      "nhrace + married.f + insur.f + stage.f + txsurg.f + txrad.f + ",
      "psa.f + nses_quin + no2_new + pm25_new + ndvi27_final"
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


long_imp_list <- lapply(long_imp_list, function(dat) {
  
  dat <- dat |>
    mutate(
      pfoa_mcl = factor(
        if_else(
          is.na(pfoa_dx),
          NA_character_,
          if_else(pfoa_dx >= 14, "At or Above MCL", "Below MCL")
        ),
        levels = c("Below MCL", "At or Above MCL")
      ),
      
      pfos_mcl = factor(
        if_else(
          is.na(pfos_dx),
          NA_character_,
          if_else(pfos_dx >= 13, "At or Above MCL", "Below MCL")
        ),
        levels = c("Below MCL", "At or Above MCL")
      ),
      
      pfna_mcl = factor(
        if_else(
          is.na(pfna_dx),
          NA_character_,
          if_else(pfna_dx >= 13, "At or Above MCL", "Below MCL")
        ),
        levels = c("Below MCL", "At or Above MCL")
      )
    )
  
  dat
})

long_imp_list[[1]] |>
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
long_imp_list[[1]] |>
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



#################################################
# CASE COUNTS AND PERSON-YEARS FOR PFNA
#################################################

dat_counts <- long_imp_list[[1]]

get_quint_counts <- function(dat, exposure, outcome) {
  
  dat %>%
    filter(
      !is.na(.data[[exposure]]),
      !is.na(.data[[outcome]])
    ) %>%
    group_by(.data[[exposure]]) %>%
    summarise(
      cases = sum(.data[[outcome]]),
      person_years = sum(
        (tstop - tstart) / 365.25,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    rename(quintile_value = all_of(exposure)) %>%
    mutate(
      exposure = exposure,
      outcome = outcome,
      quintile = paste0("Q", as.character(quintile_value))
    )
}

#################################################
# PFNA CATEGORICAL: pfna_cat
#################################################

quint_counts_pfna <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    get_quint_counts(
      dat = dat_counts,
      exposure = "pfna_cat",
      outcome = outcomes[[outcome_name]]
    ) %>%
      mutate(outcome_label = outcome_name)
  })
)

write.csv(
  quint_counts_pfna,
  "~/quint_counts_pfna.csv",
  row.names = FALSE
)


#################################################
# PFNA CONTINUOUS: pfna
#################################################

overall_counts_pfna <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    outcome_var <- outcomes[[outcome_name]]
    
    dat_counts %>%
      filter(
        !is.na(pfna),
        !is.na(.data[[outcome_var]])
      ) %>%
      summarise(
        exposure = "pfna",
        outcome_label = outcome_name,
        outcome = outcome_var,
        cases = sum(.data[[outcome_var]]),
        person_years = sum(
          (tstop - tstart) / 365.25,
          na.rm = TRUE
        )
      )
  })
)

write.csv(
  overall_counts_pfna,
  "~/overall_counts_pfna.csv",
  row.names = FALSE
)

############### running models #####################################

results_quint_pfna <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_quint, function(exposure) {
    run_pfas_cox(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list = long_imp_list,
      exposure_type = "quintile"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))

results_trend_pfna <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_trend, function(exposure) {
    run_pfas_cox(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list = long_imp_list,
      exposure_type = "trend"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))

results_iqr_pfna <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_iqr, function(exposure) {
    run_pfas_cox(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list = long_imp_list,
      exposure_type = "iqr"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))



results_mcl_pfna <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    do.call(
      rbind,
      lapply(pfas_mcl, function(exposure) {
        
        run_pfas_cox(
          exposure = exposure,
          outcome = outcomes[[outcome_name]],
          long_imp_list = long_imp_list,
          exposure_type = "mcl"
        ) |>
          mutate(outcome_label = outcome_name)
        
      })
    )
  })
)


summary(long_imp_list[[1]][, c("pfna_dx")])

final_results_pfna <- bind_rows(results_quint_pfna, results_trend_pfna, results_mcl_pfna, results_iqr_pfna) %>%
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
footnote <- final_results_pfna[1, ]
footnote[] <- NA

# Put the note in the first column
footnote[[1]] <- "Models are adjusted for age, diagnosis year, race/ethnicity, marital status, insurance status, stage, treatment (surgery and or radiation), PSA,  nSES quintile, no2, pm2.5, and ndvi27"

# Add it to the bottom
final_results_pfna <- rbind(final_results_pfna, footnote)

final_results_pfna <- final_results_pfna |>
  mutate(exposure_type = case_when(exposure_type == "quintile" ~ "tertile",
                                   TRUE ~ exposure_type))



write.csv(
  final_results_pfna,
  "~/pfas_mortality_imputed_cox_results_main_pfna.csv",
  row.names = FALSE
)


## age and year of diagnosis only 
run_pfas_cox_age_dxyr <- function(exposure, outcome, long_imp_list, exposure_type) {
  
  fit_list <- lapply(long_imp_list, function(dat) {
    
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


results_age_dxyr_quint_pfna <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_quint, function(exposure) {
    run_pfas_cox_age_dxyr(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list = long_imp_list,
      exposure_type = "quintile"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))

results_age_dxyr_trend_pfna <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_trend, function(exposure) {
    run_pfas_cox_age_dxyr(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list = long_imp_list,
      exposure_type = "trend"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))

results_iqr_age_dxyr_pfna <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
  do.call(rbind, lapply(pfas_iqr, function(exposure) {
    run_pfas_cox_age_dxyr(
      exposure = exposure,
      outcome = outcomes[[outcome_name]],
      long_imp_list = long_imp_list,
      exposure_type = "iqr"
    ) %>%
      mutate(outcome_label = outcome_name)
  }))
}))



results_mcl_age_dxyr_pfna <- bind_rows(
  lapply(names(outcomes), function(outcome_name) {
    
    bind_rows(
      lapply(pfas_mcl, function(exposure) {
        
        run_pfas_cox_age_dxyr(
          exposure = exposure,
          outcome = outcomes[[outcome_name]],
          long_imp_list = long_imp_list,
          exposure_type = "mcl"
        ) |>
          dplyr::mutate(
            outcome_label = outcome_name
          )
      })
    )
  })
)



final_results_age_dxyr_pfna <- bind_rows(
  results_age_dxyr_quint_pfna,
  results_age_dxyr_trend_pfna,
  results_mcl_age_dxyr_pfna,
  results_iqr_age_dxyr_pfna,
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
footnote <- final_results_age_dxyr_pfna[1, ]
footnote[] <- NA

# Put the note in the first column
footnote[[1]] <- "Models are adjusted for age and diagnosis year"

# Add it to the bottom
final_results_age_dxyr_pfna <- rbind(final_results_age_dxyr_pfna, footnote)

final_results_age_dxyr_pfna <- final_results_age_dxyr_pfna |>
  mutate(exposure_type = case_when(exposure_type == "quintile" ~ "tertile",
                                   TRUE ~ exposure_type))


write.csv(
  final_results_age_dxyr_pfna,
  "~/pfas_mortality_imputed_cox_age_dxyr_results_pfna.csv",
  row.names = FALSE
)


############################################################
# MUTUALLY ADJUSTED PFAS MODELS
# Fully adjusted + quintiles of the other two PFAS
############################################################

run_pfas_cox_mutual <- function(
    exposure,
    outcome,
    long_imp_list,
    exposure_type
) {
  
  # ----------------------------------------------------------
  # Determine which PFAS is the exposure of interest
  # ----------------------------------------------------------
  
  if (grepl("^pfoa", exposure)) {
    main_pfas <- "pfoa"
  } else if (grepl("^pfos", exposure)) {
    main_pfas <- "pfos"
  } else if (grepl("^pfna", exposure)) {
    main_pfas <- "pfna"
  } else {
    stop("Could not identify PFAS from exposure name.")
  }
  
  
  # ----------------------------------------------------------
  # Other PFAS are ALWAYS included as quintiles
  # ----------------------------------------------------------
  
  all_pfas_quint <- c(
    "pfoa_quin",
    "pfos_quin",
    "pfna_quin"
  )
  
  main_quint <- paste0(main_pfas, "_quin")
  
  other_pfas_quint <- setdiff(
    all_pfas_quint,
    main_quint
  )
  
  
  # ----------------------------------------------------------
  # Fit model in each imputed dataset
  # ----------------------------------------------------------
  
  fit_list <- lapply(long_imp_list, function(dat) {
    
    form <- as.formula(
      paste0(
        
        "Surv(tstart, tstop, ", outcome, ") ~ ",
        
        # PFAS exposure of interest
        exposure,
        
        # Other 2 PFAS as quintiles
        " + ",
        paste(other_pfas_quint, collapse = " + "),
        
        # Fully adjusted covariates
        " + age_cat.f + dxyr_cat.f + ",
        "nhrace + married.f + insur.f + stage.f + ",
        "txsurg.f + txrad.f + psa.f + nses_quin + ",
        "no2_new + pm25_new + ndvi27_final"
        
      )
    )
    
    
    coxph(
      form,
      data = dat
    )
    
  })
  
  
  # ----------------------------------------------------------
  # Pool results across imputations
  # ----------------------------------------------------------
  
  tab <- summary(
    pool(
      as.mira(fit_list)
    )
  )
  
  
  # ----------------------------------------------------------
  # Keep only exposure of interest
  # ----------------------------------------------------------
  
  tab <- tab %>%
    dplyr::filter(grepl(exposure, term, fixed = TRUE)) %>%
    dplyr::mutate(
      exposure = exposure,
      outcome = outcome,
      exposure_type = exposure_type,
      HR = exp(estimate),
      HR_95LL = exp(estimate - 1.96 * std.error),
      HR_95UL = exp(estimate + 1.96 * std.error)
    )
  
  
  # ----------------------------------------------------------
  # Add contrast labels
  # ----------------------------------------------------------
  
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
    
    tab <- tab %>%
      dplyr::mutate(
        contrast = "At or Above MCL vs Below MCL"
      )
    
  }
  
  
  tab <- tab %>%
    dplyr::mutate(
      contrast = as.character(contrast)
    )
  
  
  tab
}

############################################################
# Run mutually adjusted quintile models
############################################################

results_quint_mutual <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    do.call(
      rbind,
      lapply(pfas_quint, function(exposure) {
        
        run_pfas_cox_mutual(
          exposure = exposure,
          outcome = outcomes[[outcome_name]],
          long_imp_list = long_imp_list,
          exposure_type = "quintile"
        ) %>%
          mutate(
            outcome_label = outcome_name
          )
        
      })
    )
  })
)

results_trend_mutual <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    do.call(
      rbind,
      lapply(pfas_trend, function(exposure) {
        
        run_pfas_cox_mutual(
          exposure = exposure,
          outcome = outcomes[[outcome_name]],
          long_imp_list = long_imp_list,
          exposure_type = "trend"
        ) %>%
          mutate(
            outcome_label = outcome_name
          )
        
      })
    )
  })
)

results_iqr_mutual <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    do.call(
      rbind,
      lapply(pfas_iqr, function(exposure) {
        
        run_pfas_cox_mutual(
          exposure = exposure,
          outcome = outcomes[[outcome_name]],
          long_imp_list = long_imp_list,
          exposure_type = "iqr"
        ) %>%
          mutate(
            outcome_label = outcome_name
          )
        
      })
    )
  })
)


results_mcl_mutual <- do.call(
  rbind,
  lapply(names(outcomes), function(outcome_name) {
    
    do.call(
      rbind,
      lapply(pfas_mcl, function(exposure) {
        
        run_pfas_cox_mutual(
          exposure = exposure,
          outcome = outcomes[[outcome_name]],
          long_imp_list = long_imp_list,
          exposure_type = "mcl"
        ) %>%
          mutate(
            outcome_label = outcome_name
          )
        
      })
    )
  })
)



final_results_mutual <- bind_rows(
  results_quint_mutual,
  results_trend_mutual,
  results_mcl_mutual,
  results_iqr_mutual
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

write.csv(
  final_results_mutual,
  "~/pfas_mortality_mutually_adjusted_results.csv",
  row.names = FALSE
)
