#################################################################
# R Code - multiple imputation

# Megan Shanahan
# Date: 06/02/2026

# Goal: Use multiple imputation to resolve missingness
#################################################################


library(tidyverse)
library(mice)
library(survival)
library(forcats)
library(broom)
library(dplyr)
library(ggplot2)



# load data
load("~/R/registry/df_table1_diagnosis.Rda")

df_impute <- df_table1_diagnosis |>
  mutate(all_death.f = as.factor(all_death),
         cvd_death.f = as.factor(cvd_death),
         prostate_death.f = as.factor(prostate_death),
         fupyrs = tstop/12,
         pfoa_dx = pfoa, 
         pfos_dx = pfos) 

length(unique(df_table1_diagnosis$patient_study_id)) 


# recoding "Missing" categories to NA

df_impute <- df_impute |>
  mutate(
    married.f = case_when(
      married.f == "Missing" ~ NA_character_,
      TRUE ~ as.character(married.f)
    ),
    psa.f = case_when(
      psa.f == "Missing" ~ NA_character_,
      TRUE ~ as.character(psa.f)
    ),
    married.f = factor(married.f),
    psa.f = factor(psa.f)
  )


names(df_impute)

table(df_impute$dxyr_cat.f, useNA="always")
class(df_impute$dxyr_cat.f)

table(df_impute$psa_perc_cat.f)

## create trend variables for the exposure variables:
# Function to create trend variable based on quintiles
create_trend_variable <- function(df, cont_var, quint_var, new_trend_var) {
  
  df[[quint_var]] <- as.factor(df[[quint_var]])
  
  median_vals <- aggregate(
    df[[cont_var]],
    by = list(df[[quint_var]]),
    FUN = median,
    na.rm = TRUE
  )
  
  colnames(median_vals) <- c(quint_var, "median_value")
  
  df <- merge(df, median_vals, by = quint_var, all.x = TRUE, sort = FALSE)
  
  df[[new_trend_var]] <- df$median_value
  df$median_value <- NULL
  
  return(df)
}

pfas_vars <- c("pfoa", "pfos")

df_impute_exposure <- df_impute

for (x in pfas_vars) {
  cont_var <- paste0(x, "_dx")
  quint_var <- paste0(x, "_quin")
  trend_var <- paste0(x, "_dx_t")
  
  df_impute_exposure <- create_trend_variable(
    df_impute_exposure,
    cont_var,
    quint_var,
    trend_var
  )
}




## one row per person
pfas_imp_dat <- df_impute[, c(
  "patient_study_id",
  "all_death.f", "cvd_death.f", "prostate_death.f",
  "age_cat.f", "nhrace", "dxyr_cat.f",
  "married.f", "insur.f", "grade_combo.f", "stage.f",
  "txsurg.f", "txrad.f",
  "psa.f",
  "nses_quin",
  "no2_new", "pm25_new", "ndvi27_final"
)]

md.pattern(pfas_imp_dat)

meth <- make.method(pfas_imp_dat)

## do not impute ID or fully observed/non-target variables
meth["patient_study_id"] <- ""
meth[c("all_death.f", "cvd_death.f", "prostate_death.f")] <- ""

## continuous
meth[c("no2_new", "pm25_new", "ndvi27_final")] <- "pmm"

## categorical: choose based on type
meth[c("age_cat.f","dxyr_cat.f", "married.f", "insur.f", "grade_combo.f", "stage.f", "txsurg.f", "txrad.f",
          "psa.f", "nhrace", "nses_quin")] <- "lda" # Categorical variables

pred <- make.predictorMatrix(pfas_imp_dat)

## ID should not be used to predict anything
pred[, "patient_study_id"] <- 0
pred["patient_study_id", ] <- 0

imputed <- mice(
  pfas_imp_dat,
  m = 5,
  method = meth,
  predictorMatrix = pred,
  seed = 123
)

View(imputed$loggedEvents)

# after this merge each impute set back to long

## load data
load("~/R/registry/df_long_final.Rda") 

# removing the variables that were imputed from the long data
df_long <- df_long_final |>
  dplyr::select(-c("age_cat.f","dxyr_cat.f", "married.f", "insur.f", "grade_combo.f", "stage.f", "txsurg.f", "txrad.f",
                   "psa.f", "nhrace", "nses_quin", "nses_quin",
                   "no2_new", "pm25_new", "ndvi27_final")) |>
  filter( tstop > tstart) |> 
  mutate(pfoa_dx = pfoa, 
         pfos_dx = pfos) 

length(unique(df_long$patient_study_id)) 

create_trend_variable <- function(df, cont_var, quint_var, new_trend_var) {
  
  df[[quint_var]] <- as.factor(df[[quint_var]])
  
  median_vals <- aggregate(
    df[[cont_var]],
    by = list(df[[quint_var]]),
    FUN = median,
    na.rm = TRUE
  )
  
  colnames(median_vals) <- c(quint_var, "median_value")
  
  df <- merge(df, median_vals, by = quint_var, all.x = TRUE, sort = FALSE)
  
  df[[new_trend_var]] <- df$median_value
  df$median_value <- NULL
  
  return(df)
}

pfas_vars <- c("pfoa", "pfos")

#df_long_exposure <- df_long


for (x in pfas_vars) {
  cont_var <- paste0(x, "_dx")
  quint_var <- paste0(x, "_quin")
  trend_var <- paste0(x, "_dx_t")
  
  df_long <- create_trend_variable(
    df_long,
    cont_var,
    quint_var,
    trend_var
  )
}


imp_list <- complete(imputed, action = "all")

long_imp_list <- lapply(imp_list, function(x) {
  left_join(df_long, x, by = "patient_study_id")
})

sapply(long_imp_list, function(dat) {
  length(unique(dat$patient_study_id))
})



############## adding continuous PFAS per IQR term ###############

iqr_pfoa <- IQR(long_imp_list[[1]]$pfoa_dx, na.rm = TRUE)
iqr_pfos <- IQR(long_imp_list[[1]]$pfos_dx, na.rm = TRUE)

long_imp_list <- lapply(long_imp_list, function(dat) {
  dat$pfoa_iqr <- dat$pfoa_dx / iqr_pfoa
  dat$pfos_iqr <- dat$pfos_dx / iqr_pfos
  dat
})

###################### saving long_imp_list #####################

save(
  long_imp_list,
  file = "~/R/registry/long_list_imp.Rda"
)


############################# Cox models ###############################

# load("~/R/registry/long_list_imp.Rda")


# defining exposures and outcomes 
pfas_quint <- c("pfoa_quin", "pfos_quin")
pfas_trend <- c("pfoa_dx_t", "pfos_dx_t")
pfas_iqr <- c("pfoa_iqr", "pfos_iqr")
pfas_mcl   <- c("pfoa_mcl", "pfos_mcl")

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
      )
    )
  
  dat
})

long_imp_list[[1]] |>
  distinct(patient_study_id, pfoa_mcl, pfos_mcl) |>
  summarise(
    pfoa_below = sum(pfoa_mcl == "Below MCL", na.rm = TRUE),
    pfoa_above = sum(pfoa_mcl == "At or Above MCL", na.rm = TRUE),
    pfos_below = sum(pfos_mcl == "Below MCL", na.rm = TRUE),
    pfos_above = sum(pfos_mcl == "At or Above MCL", na.rm = TRUE)
  )

# seeing how many patients ever exceed the MCL
long_imp_list[[1]] |>
  group_by(patient_study_id) |>
  summarise(
    pfoa_above = any(pfoa_mcl == "At or Above MCL", na.rm = TRUE),
    pfos_above = any(pfos_mcl == "At or Above MCL", na.rm = TRUE)
  ) |>
  summarise(
    pfoa_above = sum(pfoa_above),
    pfos_above = sum(pfos_above)
  )


################################################# CASE COUNTS AND PERSON YEARS ###########################
load("~/R/registry/long_list_imp.Rda")

dat_counts <- long_imp_list[[1]]

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
  "~/quint_counts.csv",
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
  "~/overall_counts.csv",
  row.names = FALSE
)

# by pfas
pfas <- c("pfoa", "pfos")

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
  "~/pfas_continuous_counts.csv",
  row.names = FALSE
)


############### running models #####################################

results_quint <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
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

results_trend <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
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

results_iqr <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
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
  


results_mcl <- do.call(
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


summary(long_imp_list[[1]][, c("pfoa_dx", "pfos_dx")])

final_results <- bind_rows(results_quint, results_trend, results_mcl, results_iqr) %>%
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
  "~/pfas_mortality_imputed_cox_results_main.csv",
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


results_age_dxyr_quint <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
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

results_age_dxyr_trend <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
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

results_iqr_age_dxyr <- do.call(rbind, lapply(names(outcomes), function(outcome_name) {
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



results_mcl_age_dxyr <- bind_rows(
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



final_results_age_dxyr <- bind_rows(
  results_age_dxyr_quint,
  results_age_dxyr_trend,
  results_mcl_age_dxyr,
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
  "~/pfas_mortality_imputed_cox_age_dxyr_results.csv",
  row.names = FALSE
)
