#################################################################
# R Code - IPW - Cox-model-based IPCW estimator

# Megan Shanahan
# Date: 07/17/2026

# Updated: Hari Iyer - added p-trend
# Date: 09/28/2026

# Goal: Use ipcw for cvd and prostate cancer mortality with a Cox-model-based IPCW estimator
#################################################################

##############################################################################
# COX-MODEL-BASED IPCW ANALYSES
#
# Outcomes:
#   1. Cardiovascular disease mortality
#   2. Prostate cancer mortality
#
# Exposures:
#   1. PFAS quintiles
#   2. PFAS quintile-median trend variables
#   3. PFAS MCL indicators
#
# Memory strategy:
#   - Run age + diagnosis-year models in one R session
#   - Restart R
#   - Run fully adjusted models in another R session
#   - Process one imputed dataset at a time
##############################################################################


##############################################################################
# Select which outcome-model adjustment set to run
##############################################################################

# Use "age_dxyr" for age + diagnosis-year models
# Use "full" for fully adjusted models

RUN_MODEL <- "age_dxyr"

if (!RUN_MODEL %in% c("age_dxyr", "full")) {
  stop('RUN_MODEL must be either "age_dxyr" or "full".')
}


##############################################################################
# Packages
##############################################################################

library(survival)
library(mice)
library(dplyr)
library(tibble)


##############################################################################
# Load long multiply imputed datasets
##############################################################################

## load("~/R/registry/long_list_imp.Rda")
load("/extra/home/hi97/R/registry/nj_res_hx/long_list_imp.Rda")


# Confirm that the expected object was loaded
if (!exists("long_imp_list")) {
  stop("The object long_imp_list was not found after loading the Rda file.")
}

cat(
  "Number of imputed datasets:",
  length(long_imp_list),
  "\n"
)

cat(
  "Unique participants in each imputation:\n"
)

print(
  sapply(
    long_imp_list,
    function(dat) {
      dplyr::n_distinct(dat$patient_study_id)
    }
  )
)


##############################################################################
# Create continuous PFAS per-IQR variables
##############################################################################


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





##############################################################################
# Define PFAS exposure variables
##############################################################################

pfas_quint <- c("pfna_cat")
pfas_trend <- c( "pfna_dx_t")
pfas_iqr <- c("pfna_iqr")
pfas_mcl   <- c("pfna_mcl")


##############################################################################
# Create PFAS MCL indicators
##############################################################################


long_imp_list <- lapply(
  long_imp_list,
  function(dat) {
    
    dat %>%
      mutate(
        pfna_mcl = factor(
          case_when(
            is.na(pfna_dx) ~ NA_character_,
            pfna_dx >= 13 ~ "At or Above MCL",
            pfna_dx < 13 ~ "Below MCL"
          ),
          levels = c(
            "Below MCL",
            "At or Above MCL"
          )
        )
      )
  }
)


##############################################################################
# Outcomes requiring IPCW
#
# All-cause mortality is excluded because there is no competing cause when
# all deaths are counted as outcome events.
##############################################################################

outcomes_ipcw <- c(
  cvd = "cvd_death",
  prostate = "prostate_death"
)


##############################################################################
# Covariates for denominator censoring model
##############################################################################

ipcw_covariates <- c(
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


##############################################################################
# Helper function:
# Evaluate a cumulative baseline hazard at requested follow-up times
##############################################################################

get_cumhaz_at_time <- function(
    basehaz_data,
    times
) {
  
  if ("hazard" %in% names(basehaz_data)) {
    
    hazard_variable <- "hazard"
    
  } else if ("cumhaz" %in% names(basehaz_data)) {
    
    hazard_variable <- "cumhaz"
    
  } else {
    
    stop(
      "Could not identify the cumulative-hazard column returned by basehaz()."
    )
  }
  
  # For each requested time, identify the most recent censoring-event time
  # that occurred at or before that time.
  index <- findInterval(
    times,
    basehaz_data$time
  )
  
  cumulative_hazard <- numeric(
    length(times)
  )
  
  has_prior_event_time <- index > 0
  
  cumulative_hazard[has_prior_event_time] <-
    basehaz_data[[hazard_variable]][
      index[has_prior_event_time]
    ]
  
  cumulative_hazard
}


##############################################################################
# Function to calculate Cox-model-based stabilized IPC weights
##############################################################################

calculate_ipcw_cox <- function(
    data,
    outcome,
    exposure,
    covariates,
    id_var = "patient_study_id",
    tstart_var = "tstart",
    tstop_var = "tstop",
    trim_probs = c(0.01, 0.99),
    probability_bound = 0.0001
) {
  
  ##########################################################################
  # Check required variables
  ##########################################################################
  
  required_variables <- unique(
    c(
      id_var,
      tstart_var,
      tstop_var,
      "all_death",
      "cvd_death",
      "prostate_death",
      exposure,
      covariates
    )
  )
  
  missing_variables <- setdiff(
    required_variables,
    names(data)
  )
  
  if (length(missing_variables) > 0) {
    
    stop(
      "Missing required variables: ",
      paste(
        missing_variables,
        collapse = ", "
      )
    )
  }
  
  
  ##########################################################################
  # Preserve original row order
  ##########################################################################
  
  data <- data %>%
    mutate(
      .original_row_order = row_number()
    )
  
  
  ##########################################################################
  # Define competing-cause censoring event
  ##########################################################################
  
  if (outcome == "cvd_death") {
    
    # A death from any non-CVD cause is a competing censoring event.
    data <- data %>%
      mutate(
        competing_censor = case_when(
          all_death == 1 & cvd_death == 0 ~ 1L,
          TRUE ~ 0L
        )
      )
    
  } else if (outcome == "prostate_death") {
    
    # A death from any non-prostate-cancer cause is a competing censoring event.
    data <- data %>%
      mutate(
        competing_censor = case_when(
          all_death == 1 & prostate_death == 0 ~ 1L,
          TRUE ~ 0L
        )
      )
    
  } else {
    
    stop(
      "Outcome must be cvd_death or prostate_death."
    )
  }
  
  
  ##########################################################################
  # Sort intervals chronologically within participant
  ##########################################################################
  
  data <- data %>%
    arrange(
      .data[[id_var]],
      .data[[tstart_var]],
      .data[[tstop_var]]
    )
  
  
  ##########################################################################
  # Check for invalid survival intervals
  ##########################################################################
  
  invalid_interval <- (
    is.na(data[[tstart_var]]) |
      is.na(data[[tstop_var]]) |
      data[[tstop_var]] <= data[[tstart_var]]
  )
  
  if (any(invalid_interval)) {
    
    stop(
      sum(invalid_interval),
      " rows have missing or invalid start-stop intervals."
    )
  }
  
  
  ##########################################################################
  # Check for missing censoring-model values
  ##########################################################################
  
  censoring_model_variables <- unique(
    c(
      id_var,
      tstart_var,
      tstop_var,
      "competing_censor",
      exposure,
      covariates
    )
  )
  
  missing_counts <- sapply(
    data[censoring_model_variables],
    function(x) {
      sum(is.na(x))
    }
  )
  
  if (any(missing_counts > 0)) {
    
    missing_message <- paste(
      names(missing_counts)[missing_counts > 0],
      missing_counts[missing_counts > 0],
      sep = "=",
      collapse = ", "
    )
    
    stop(
      "Missing values remain in censoring-model variables: ",
      missing_message
    )
  }
  
  
  ##########################################################################
  # Denominator censoring model
  #
  # Models the conditional hazard of a competing-cause death using:
  #   - PFAS exposure
  #   - all measured censoring-model covariates
  ##########################################################################
  
  denominator_variables <- unique(
    c(
      exposure,
      covariates
    )
  )
  
  denominator_formula <- as.formula(
    paste0(
      "Surv(",
      tstart_var,
      ", ",
      tstop_var,
      ", competing_censor) ~ ",
      paste(
        denominator_variables,
        collapse = " + "
      )
    )
  )
  
  denominator_fit <- coxph(
    formula = denominator_formula,
    data = data,
    ties = "efron",
    na.action = na.fail,
    model = FALSE,
    x = TRUE,
    y = TRUE
  )
  
  
  ##########################################################################
  # Numerator censoring model
  #
  # Exposure is retained in the numerator stabilization model.
  ##########################################################################
  
  numerator_formula <- as.formula(
    paste0(
      "Surv(",
      tstart_var,
      ", ",
      tstop_var,
      ", competing_censor) ~ ",
      exposure
    )
  )
  
  numerator_fit <- coxph(
    formula = numerator_formula,
    data = data,
    ties = "efron",
    na.action = na.fail,
    model = FALSE,
    x = TRUE,
    y = TRUE
  )
  
  
  ##########################################################################
  # Extract baseline cumulative hazards
  ##########################################################################
  
  denominator_basehaz <- basehaz(
    denominator_fit,
    centered = FALSE
  )
  
  numerator_basehaz <- basehaz(
    numerator_fit,
    centered = FALSE
  )
  
  
  ##########################################################################
  # Evaluate baseline cumulative hazard at each interval boundary
  ##########################################################################
  
  denominator_h0_start <- get_cumhaz_at_time(
    denominator_basehaz,
    data[[tstart_var]]
  )
  
  denominator_h0_stop <- get_cumhaz_at_time(
    denominator_basehaz,
    data[[tstop_var]]
  )
  
  numerator_h0_start <- get_cumhaz_at_time(
    numerator_basehaz,
    data[[tstart_var]]
  )
  
  numerator_h0_stop <- get_cumhaz_at_time(
    numerator_basehaz,
    data[[tstop_var]]
  )
  
  
  ##########################################################################
  # Obtain linear predictors
  ##########################################################################
  
  denominator_lp <- predict(
    denominator_fit,
    newdata = data,
    type = "lp",
    reference = "zero"
  )
  
  numerator_lp <- predict(
    numerator_fit,
    newdata = data,
    type = "lp",
    reference = "zero"
  )
  
  
  ##########################################################################
  # Calculate the censoring cumulative-hazard increment during each interval
  ##########################################################################
  
  denominator_hazard_increment <- pmax(
    denominator_h0_stop -
      denominator_h0_start,
    0
  )
  
  numerator_hazard_increment <- pmax(
    numerator_h0_stop -
      numerator_h0_start,
    0
  )
  
  
  ##########################################################################
  # Estimate probability of remaining uncensored during each interval
  ##########################################################################
  
  data <- data %>%
    mutate(
      p_uncensored_denom_interval = exp(
        -denominator_hazard_increment *
          exp(denominator_lp)
      ),
      
      p_uncensored_numer_interval = exp(
        -numerator_hazard_increment *
          exp(numerator_lp)
      )
    )
  
  
  ##########################################################################
  # Stop if predictions are invalid
  ##########################################################################
  
  invalid_probability <- (
    is.na(data$p_uncensored_denom_interval) |
      is.na(data$p_uncensored_numer_interval) |
      !is.finite(data$p_uncensored_denom_interval) |
      !is.finite(data$p_uncensored_numer_interval)
  )
  
  if (any(invalid_probability)) {
    
    stop(
      sum(invalid_probability),
      " invalid censoring survival probabilities were generated."
    )
  }
  
  
  ##########################################################################
  # Bound probabilities away from zero and one
  ##########################################################################
  
  data <- data %>%
    mutate(
      p_uncensored_denom_interval = pmin(
        pmax(
          p_uncensored_denom_interval,
          probability_bound
        ),
        1
      ),
      
      p_uncensored_numer_interval = pmin(
        pmax(
          p_uncensored_numer_interval,
          probability_bound
        ),
        1
      ),
      
      interval_ipcw =
        p_uncensored_numer_interval /
        p_uncensored_denom_interval
    )
  
  
  ##########################################################################
  # Construct cumulative stabilized weights
  #
  # The weight used for an interval is the cumulative weight at the START
  # of that interval. This avoids incorporating censoring information that
  # occurs during the current interval into its own analysis weight.
  ##########################################################################
  
  data <- data %>%
    group_by(
      .data[[id_var]]
    ) %>%
    mutate(
      cumulative_ipcw_end = cumprod(
        interval_ipcw
      ),
      
      ipcw = lag(
        cumulative_ipcw_end,
        default = 1
      )
    ) %>%
    ungroup()
  
  
  ##########################################################################
  # Trim cumulative weights at selected quantiles
  ##########################################################################
  
  weight_quantiles <- quantile(
    data$ipcw,
    probs = trim_probs,
    na.rm = TRUE,
    names = TRUE
  )
  
  lower_limit <- unname(
    weight_quantiles[1]
  )
  
  upper_limit <- unname(
    weight_quantiles[2]
  )
  
  data <- data %>%
    mutate(
      ipcw_trim = pmin(
        pmax(
          ipcw,
          lower_limit
        ),
        upper_limit
      )
    )
  
  
  ##########################################################################
  # Weight diagnostics
  ##########################################################################
  
  weight_diagnostics <- tibble(
    outcome = outcome,
    exposure = exposure,
    
    n_rows = nrow(data),
    
    n_patients = n_distinct(
      data[[id_var]]
    ),
    
    n_competing_censor_events = sum(
      data$competing_censor == 1,
      na.rm = TRUE
    ),
    
    untrimmed_min = min(
      data$ipcw,
      na.rm = TRUE
    ),
    
    untrimmed_mean = mean(
      data$ipcw,
      na.rm = TRUE
    ),
    
    untrimmed_median = median(
      data$ipcw,
      na.rm = TRUE
    ),
    
    untrimmed_max = max(
      data$ipcw,
      na.rm = TRUE
    ),
    
    trim_lower = lower_limit,
    trim_upper = upper_limit,
    
    trimmed_min = min(
      data$ipcw_trim,
      na.rm = TRUE
    ),
    
    trimmed_mean = mean(
      data$ipcw_trim,
      na.rm = TRUE
    ),
    
    trimmed_median = median(
      data$ipcw_trim,
      na.rm = TRUE
    ),
    
    trimmed_max = max(
      data$ipcw_trim,
      na.rm = TRUE
    ),
    
    n_untrimmed_gt_5 = sum(
      data$ipcw > 5,
      na.rm = TRUE
    ),
    
    n_untrimmed_gt_10 = sum(
      data$ipcw > 10,
      na.rm = TRUE
    ),
    
    n_untrimmed_gt_20 = sum(
      data$ipcw > 20,
      na.rm = TRUE
    ),
    
    n_untrimmed_gt_100 = sum(
      data$ipcw > 100,
      na.rm = TRUE
    )
  )
  
  
  ##########################################################################
  # Return data to original order
  ##########################################################################
  
  data <- data %>%
    arrange(
      .original_row_order
    ) %>%
    select(
      -.original_row_order
    )
  
  
  ##########################################################################
  # Remove large temporary model objects before returning
  ##########################################################################
  
  rm(
    denominator_fit,
    numerator_fit,
    denominator_basehaz,
    numerator_basehaz,
    denominator_lp,
    numerator_lp,
    denominator_h0_start,
    denominator_h0_stop,
    numerator_h0_start,
    numerator_h0_stop
  )
  
  gc(
    verbose = FALSE
  )
  
  
  ##########################################################################
  # Return weighted dataset and diagnostics
  ##########################################################################
  
  list(
    data = data,
    quantiles = weight_quantiles,
    diagnostics = weight_diagnostics
  )
}


##############################################################################
# Assign readable exposure contrast labels
##############################################################################

add_contrast_labels <- function(
    tab,
    exposure_type
) {
  
  if (exposure_type == "quintile") {
    
    tab <- tab %>%
      arrange(
        as.character(term)
      ) %>%
      mutate(
        contrast = case_when(
          row_number() == 1 ~ "Q2 vs Q1",
          row_number() == 2 ~ "Q3 vs Q1",
          row_number() == 3 ~ "Q4 vs Q1",
          row_number() == 4 ~ "Q5 vs Q1",
          TRUE ~ as.character(term)
        )
      )
    
  } else if (exposure_type == "trend") {
    
    tab <- tab %>%
      mutate(
        contrast = "Trend"
      )
    
  } else if (exposure_type == "iqr") {
    
    tab <- tab %>%
      mutate(
        contrast = "Per IQR increase"
      )
    
  } else if (exposure_type == "mcl") {
    
    tab <- tab %>%
      mutate(
        contrast = "At or Above MCL vs Below MCL"
      )
    
  } else {
    
    tab <- tab %>%
      mutate(
        contrast = as.character(term)
      )
  }
  
  tab %>%
    mutate(
      contrast = as.character(contrast)
    )
}


##############################################################################
# Run one exposure-outcome combination
#
# This function:
#   - processes one imputation at a time
#   - creates weights
#   - immediately fits the requested outcome model
#   - deletes the weighted dataset before moving to the next imputation
##############################################################################

run_pfas_cox_ipcw <- function(
    exposure,
    outcome,
    long_imp_list,
    exposure_type,
    model_type = c(
      "age_dxyr",
      "full"
    )
) {
  
  model_type <- match.arg(
    model_type
  )
  
  
  ##########################################################################
  # Define requested outcome-model formula
  ##########################################################################
  
  if (model_type == "age_dxyr") {
    
    outcome_formula <- as.formula(
      paste0(
        "Surv(tstart, tstop, ",
        outcome,
        ") ~ ",
        exposure,
        " + age_cat.f + dxyr_cat.f",
        " + cluster(patient_study_id)"
      )
    )
    
    model_label <- "age + diagnosis year"
    
  } else {
    
    outcome_formula <- as.formula(
      paste0(
        "Surv(tstart, tstop, ",
        outcome,
        ") ~ ",
        exposure,
        " + age_cat.f + dxyr_cat.f + ",
        "nhrace + married.f + insur.f + stage.f + ",
        "txsurg.f + txrad.f + psa.f + ",
        "nses_quin + no2_new + pm25_new + ndvi27_final",
        " + cluster(patient_study_id)"
      )
    )
    
    model_label <- "fully adjusted"
  }
  
  
  ##########################################################################
  # Create containers for fitted models and weight diagnostics
  ##########################################################################
  
  number_imputations <- length(
    long_imp_list
  )
  
  fit_list <- vector(
    mode = "list",
    length = number_imputations
  )
  
  diagnostics_list <- vector(
    mode = "list",
    length = number_imputations
  )
  
  
  ##########################################################################
  # Process one imputation at a time
  ##########################################################################
  
  for (
    imputation_number in seq_along(long_imp_list)
  ) {
    
    message(
      paste0(
        "\nImputation ",
        imputation_number,
        "/",
        number_imputations,
        "; outcome = ",
        outcome,
        "; exposure = ",
        exposure,
        "; model = ",
        model_type
      )
    )
    
    
    ########################################################################
    # Estimate Cox-model-based IPC weights
    ########################################################################
    
    weighted_object <- calculate_ipcw_cox(
      data = long_imp_list[[imputation_number]],
      outcome = outcome,
      exposure = exposure,
      covariates = ipcw_covariates,
      id_var = "patient_study_id",
      tstart_var = "tstart",
      tstop_var = "tstop",
      trim_probs = c(
        0.01,
        0.99
      ),
      probability_bound = 0.0001
    )
    
    weighted_data <- weighted_object$data
    
    diagnostics_list[[imputation_number]] <-
      weighted_object$diagnostics %>%
      mutate(
        imputation = imputation_number,
        model = model_label
      )
    
    
    ########################################################################
    # Fit weighted outcome model
    ########################################################################
    
    fit_list[[imputation_number]] <- coxph(
      formula = outcome_formula,
      data = weighted_data,
      weights = ipcw_trim,
      robust = TRUE,
      ties = "efron",
      na.action = na.fail,
      
      # Avoid storing the full analysis frame inside each fitted model.
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
    
    
    ########################################################################
    # Delete the large weighted dataset before the next imputation
    ########################################################################
    
    rm(
      weighted_object,
      weighted_data
    )
    
    gc(
      verbose = FALSE
    )
  }
  
  
  ##########################################################################
  # Pool outcome-model coefficients across imputations
  ##########################################################################
  
  pooled_model <- mice::pool(
    mice::as.mira(
      fit_list
    )
  )
  
  model_results <- summary(
    pooled_model
  ) %>%
    filter(
      startsWith(
        as.character(term),
        exposure
      )
    ) %>%
    mutate(
      exposure = exposure,
      outcome = outcome,
      exposure_type = exposure_type,
      model = model_label,
      
      HR = exp(
        estimate
      ),
      
      HR_95LL = exp(
        estimate -
          1.96 * std.error
      ),
      
      HR_95UL = exp(
        estimate +
          1.96 * std.error
      )
    )
  
  
  ##########################################################################
  # Add readable exposure contrast labels
  ##########################################################################
  
  model_results <- add_contrast_labels(
    tab = model_results,
    exposure_type = exposure_type
  )
  
  
  ##########################################################################
  # Combine weight diagnostics across imputations
  ##########################################################################
  
  weight_diagnostics <- bind_rows(
    diagnostics_list
  )
  
  
  ##########################################################################
  # Remove fitted models after pooling
  ##########################################################################
  
  rm(
    fit_list,
    pooled_model,
    diagnostics_list
  )
  
  gc(
    verbose = FALSE
  )
  
  
  ##########################################################################
  # Return results and diagnostics
  ##########################################################################
  
  list(
    results = model_results,
    diagnostics = weight_diagnostics
  )
}


##############################################################################
# Run all PFAS models for one adjustment level
##############################################################################

run_all_ipcw_models <- function(
    long_imp_list,
    model_type
) {
  
  exposure_sets <- list(
    quintile = pfas_quint,
    trend = pfas_trend,
    iqr = pfas_iqr,
    mcl = pfas_mcl
  )
  
  result_list <- list()
  diagnostics_list <- list()
  
  result_number <- 1L
  
  
  ##########################################################################
  # Loop through exposure types, outcomes and individual PFAS variables
  ##########################################################################
  
  for (
    exposure_type in names(exposure_sets)
  ) {
    
    exposures <- exposure_sets[[
      exposure_type
    ]]
    
    for (
      outcome_name in names(outcomes_ipcw)
    ) {
      
      outcome_variable <- outcomes_ipcw[[
        outcome_name
      ]]
      
      for (
        exposure in exposures
      ) {
        
        message(
          paste0(
            "\n============================================================\n",
            "Starting: ",
            outcome_name,
            " / ",
            exposure,
            " / ",
            exposure_type,
            " / ",
            model_type,
            "\n============================================================"
          )
        )
        
        model_output <- run_pfas_cox_ipcw(
          exposure = exposure,
          outcome = outcome_variable,
          long_imp_list = long_imp_list,
          exposure_type = exposure_type,
          model_type = model_type
        )
        
        result_list[[result_number]] <-
          model_output$results %>%
          mutate(
            outcome_label = outcome_name
          )
        
        diagnostics_list[[result_number]] <-
          model_output$diagnostics %>%
          mutate(
            outcome_label = outcome_name,
            exposure_type = exposure_type
          )
        
        rm(
          model_output
        )
        
        gc(
          verbose = FALSE
        )
        
        result_number <- result_number + 1L
      }
    }
  }
  
  
  ##########################################################################
  # Return all model results and all weight diagnostics
  ##########################################################################
  
  list(
    results = bind_rows(
      result_list
    ),
    
    diagnostics = bind_rows(
      diagnostics_list
    )
  )
}


##############################################################################
# Organize final results
##############################################################################

organize_ipcw_results <- function(
    results
) {
  
  results %>%
    select(
      outcome_label,
      outcome,
      exposure,
      exposure_type,
      contrast,
      model,
      HR,
      HR_95LL,
      HR_95UL,
      p.value
    ) %>%
    mutate(
      outcome_label = factor(
        outcome_label,
        levels = c(
          "cvd",
          "prostate"
        ),
        labels = c(
          "Cardiovascular mortality",
          "Prostate cancer mortality"
        )
      ),
      
      pfas = case_when(
        grepl(
          "^pfna",
          exposure
        ) ~ "PFNA",
        
        TRUE ~ NA_character_
      ),
      
      pfas = factor(
        pfas,
        levels = c(
          "PFNA"
        )
      ),
      
      exposure_type = factor(
        exposure_type,
        levels = c(
          "quintile",
          "trend",
          "iqr",
          "mcl"
        ),
        labels = c(
          "Quintile",
          "Trend",
          "Per IQR",
          "MCL"
        )
      ),
      
      contrast_order = case_when(
        contrast == "Q2 vs Q1" ~ 1,
        contrast == "Q3 vs Q1" ~ 2,
        contrast == "Q4 vs Q1" ~ 3,
        contrast == "Q5 vs Q1" ~ 4,
        contrast %in% c(
          "trend",
          "Trend"
        ) ~ 1,
        contrast == "Per IQR increase" ~ 1,
        contrast == "At or Above MCL vs Below MCL" ~ 1,
        TRUE ~ 99
      )
    ) %>%
    arrange(
      exposure_type,
      outcome_label,
      pfas,
      contrast_order
    ) %>%
    select(
      outcome = outcome_label,
      pfas,
      exposure_type,
      contrast,
      HR,
      HR_95LL,
      HR_95UL,
      p.value,
      model
    )
}


##############################################################################
# Run selected model adjustment level
##############################################################################

all_ipcw_output <- run_all_ipcw_models(
  long_imp_list = long_imp_list,
  model_type = RUN_MODEL
)


names(all_ipcw_output)

organized_ipcw_results <- organize_ipcw_results(
  all_ipcw_output$results
)

weight_diagnostics <- all_ipcw_output$diagnostics


##############################################################################
# Save age + diagnosis-year results
##############################################################################

if (RUN_MODEL == "age_dxyr") {
  
  ipcw_age_dxyr_results <- organized_ipcw_results %>%
    filter(
      model == "age + diagnosis year"
    ) %>%
    select(
      -model
    )
  
  print(
    ipcw_age_dxyr_results
  )
  
  write.csv(
    ipcw_age_dxyr_results,
    file = "~/pfas_cox_ipcw_age_dxyr_results_pfna.csv",
    row.names = FALSE
  )
  
  write.csv(
    weight_diagnostics,
    file = "~/pfas_cox_ipcw_age_dxyr_weight_diagnostics_pfna.csv",
    row.names = FALSE
  )
  
  save(
    ipcw_age_dxyr_results,
    weight_diagnostics,
    file = "~/pfas_cox_ipcw_age_dxyr_results_pfna.Rda"
  )
}


##############################################################################
# Save fully adjusted results
##############################################################################

if (RUN_MODEL == "full") {
  
  ipcw_fully_adjusted_results <- organized_ipcw_results %>%
    filter(
      model == "fully adjusted"
    ) %>%
    select(
      -model
    )
  
  print(
    ipcw_fully_adjusted_results
  )
  
  write.csv(
    ipcw_fully_adjusted_results,
    file = "~/pfas_cox_ipcw_fully_adjusted_results_pfna.csv",
    row.names = FALSE
  )
  
  write.csv(
    weight_diagnostics,
    file = "~/pfas_cox_ipcw_fully_adjusted_weight_diagnostics_pfna.csv",
    row.names = FALSE
  )
  
  save(
    ipcw_fully_adjusted_results,
    weight_diagnostics,
    file = "~/pfas_cox_ipcw_fully_adjusted_results_pfna.Rda"
  )
}


##############################################################################
# Display summary of weight diagnostics
##############################################################################

weight_diagnostics_summary <- weight_diagnostics %>%
  group_by(
    model,
    outcome_label,
    exposure_type,
    exposure
  ) %>%
  summarise(
    number_imputations = n(),
    
    mean_trim_lower = mean(
      trim_lower,
      na.rm = TRUE
    ),
    
    mean_trim_upper = mean(
      trim_upper,
      na.rm = TRUE
    ),
    
    mean_untrimmed_weight = mean(
      untrimmed_mean,
      na.rm = TRUE
    ),
    
    largest_untrimmed_weight = max(
      untrimmed_max,
      na.rm = TRUE
    ),
    
    mean_trimmed_weight = mean(
      trimmed_mean,
      na.rm = TRUE
    ),
    
    largest_trimmed_weight = max(
      trimmed_max,
      na.rm = TRUE
    ),
    
    total_weights_gt_10 = sum(
      n_untrimmed_gt_10,
      na.rm = TRUE
    ),
    
    total_weights_gt_100 = sum(
      n_untrimmed_gt_100,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(
  weight_diagnostics_summary
)


##############################################################################
# Save diagnostic summary
##############################################################################

if (RUN_MODEL == "age_dxyr") {
  
  write.csv(
    weight_diagnostics_summary,
    file = "~/pfas_cox_ipcw_age_dxyr_weight_diagnostics_summary_pfna.csv",
    row.names = FALSE
  )
  
} else {
  
  write.csv(
    weight_diagnostics_summary,
    file = "~/pfas_cox_ipcw_fully_adjusted_weight_diagnostics_summary_pfna.csv",
    row.names = FALSE
  )
}


##############################################################################
# Clean up
##############################################################################

rm(
  all_ipcw_output,
  organized_ipcw_results
)

gc()