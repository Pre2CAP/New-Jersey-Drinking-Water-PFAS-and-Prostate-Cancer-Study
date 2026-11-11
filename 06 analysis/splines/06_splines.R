#################################################################
# R Code - multiple imputation

# Megan Shanahan
# Date: 07/14/2026

# Goal: Using splines for non-linearity
#################################################################

# load data
load("~/R/registry/long_list_imp.Rda")


################### Define PFAS variables and labels #####################

library(survival)
library(mice)
library(miceadds)
library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(tidyr)

pfas_spline_vars <- c(
  pfoa_dx = "PFOA (ng/L)",
  pfos_dx = "PFOS (ng/L)"
#  pfna_dx = "PFNA (ng/L)"
)

outcomes <- c(
  all_cause = "all_death",
  cvd       = "cvd_death",
  prostate  = "prostate_death"
)

outcome_titles <- c(
  all_cause = "All-cause mortality",
  cvd       = "Cardiovascular mortality",
  prostate  = "Prostate cancer mortality"
)

################### Function to create one spline plot ##################
get_pfas_spline_plot <- function(
    data,
    exposure,
    outcome,
    x_label,
    outcome_title,
    spline_df = 3,
    reference = c("median", "zero")
) {
  
  reference <- match.arg(reference)
  
  # Complete observations for the exposure and survival outcome
  model_data <- data |>
    filter(
      !is.na(.data[[exposure]]),
      !is.na(.data[[outcome]]),
      !is.na(tstart),
      !is.na(tstop),
      tstop > tstart
    )
  
  # Construct the Cox model formula
  spline_formula <- as.formula(
    paste0(
      "Surv(tstart, tstop, ", outcome, ") ~ ",
      "pspline(", exposure, ", df = ", spline_df, ") + ",
      "age_cat.f + dxyr_cat.f + ",
      "nhrace + married.f + insur.f + stage.f + ",
      "txsurg.f + txrad.f + psa.f + ",
      "nses_quin + no2_new + pm25_new + ndvi27_final"
    )
  )
  
  # Fit spline Cox model
  spline_fit <- coxph(
    spline_formula,
    data = model_data,
    x = TRUE,
    model = TRUE
  )
  
  # Extract the estimated spline function
  spline_term <- termplot(
    spline_fit,
    term = 1,
    se = TRUE,
    plot = FALSE
  )[[1]]
  
  x_values  <- spline_term$x
  log_hr    <- spline_term$y
  se_values <- spline_term$se
  
  # Choose the reference concentration
  if (reference == "median") {
    reference_value <- median(
      model_data[[exposure]],
      na.rm = TRUE
    )
  } else {
    reference_value <- 0
  }
  
  # Find the plotted value closest to the desired reference
  reference_index <- which.min(
    abs(x_values - reference_value)
  )
  
  plotted_reference <- x_values[reference_index]
  
  # Center the spline so HR = 1 at the reference value
  centered_log_hr <- log_hr - log_hr[reference_index]
  
  plot_data <- tibble(
    exposure_value = x_values,
    HR = exp(centered_log_hr),
    HR_lower = exp(centered_log_hr - 1.96 * se_values),
    HR_upper = exp(centered_log_hr + 1.96 * se_values)
  )
  
  # One exposure value per participant for the rug plot
  #
  # This avoids giving participants with longer follow-up more rug marks
  # simply because they have more rows in the long dataset.
  rug_data <- model_data |>
    distinct(
      patient_study_id,
      exposure_value = .data[[exposure]]
    )
  
  ggplot(
    plot_data,
    aes(x = exposure_value, y = HR)
  ) +
    geom_ribbon(
      aes(
        ymin = HR_lower,
        ymax = HR_upper
      ),
      alpha = 0.20
    ) +
    geom_line(
      linewidth = 1
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = plotted_reference,
      linetype = "dotted"
    ) +
    geom_rug(
      data = rug_data,
      aes(x = exposure_value),
      sides = "b",
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    labs(
      title = outcome_title,
      subtitle = paste0(
        "Reference = ",
        round(plotted_reference, 2),
        " (ng/L)"
      ),
      x = x_label,
      y = "Hazard ratio"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        size = 11,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 9,
        hjust = 0.5
      ),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

################# Create plots from the first imputation #####################

plot_imp <- long_imp_list[[1]]

## All-cause mortality plots
all_cause_spline_plots <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline <- wrap_plots(
  all_cause_spline_plots,
  nrow = 1
)

print(combined_all_cause_spline)

ggsave(
  filename = "pfas_all_cause_spline.png",
  plot = combined_all_cause_spline,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots 

cvd_spline_plots <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline <- wrap_plots(
  cvd_spline_plots,
  nrow = 1
)

print(combined_cvd_spline)

ggsave(
  filename = "pfas_cvd_spline.png",
  plot = combined_cvd_spline,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer

prostate_spline_plots <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_prostate_spline <- wrap_plots(
  prostate_spline_plots,
  nrow = 1
)

print(combined_prostate_spline)

ggsave(
  filename = "pfas_prostate_spline.png",
  plot = combined_prostate_spline,
  width = 15,
  height = 4.5,
  dpi = 300
)


##################### using log transformed pfas ##############################

pfas_spline_vars_log <- c(
  pfoa_log = "Log PFOA",
  pfos_log = "Log PFOS"
)

## All-cause mortality plots (log)
all_cause_spline_plots_log <- lapply(
  names(pfas_spline_vars_log),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars_log[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline_log <- wrap_plots(
  all_cause_spline_plots_log,
  nrow = 1
)

print(combined_all_cause_spline_log)

ggsave(
  filename = "pfas_all_cause_spline_log.png",
  plot = combined_all_cause_spline_log,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots (log)

cvd_spline_plots_log <- lapply(
  names(pfas_spline_vars_log),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars_log[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline_log <- wrap_plots(
  cvd_spline_plots_log,
  nrow = 1
)

print(combined_cvd_spline_log)

ggsave(
  filename = "pfas_cvd_spline_log.png",
  plot = combined_cvd_spline_log,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer (log)

prostate_spline_plots_log <- lapply(
  names(pfas_spline_vars_log),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars_log[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_prostate_spline_log <- wrap_plots(
  prostate_spline_plots_log,
  nrow = 1
)

print(combined_prostate_spline_log)

ggsave(
  filename = "pfas_prostate_spline_log.png",
  plot = combined_prostate_spline_log,
  width = 15,
  height = 4.5,
  dpi = 300
)


########################### STRATIFIED SPLINES #############################################

plot_imp <- long_imp_list[[1]]

# less than 65
plot_imp_lt65 <- plot_imp |>
  dplyr::filter(dxage < 65)

# 65 +
plot_imp_65plus <- plot_imp |>
  dplyr::filter(dxage >= 65)

# Non-hispanic white
plot_imp_nhw <- plot_imp |>
  dplyr::filter(nhrace == "NHW")

# Non - white
plot_imp_nonwhite <- plot_imp |>
  dplyr::filter(nhrace != "NHW", !is.na(nhrace))

# Non-hispanic black 
plot_imp_nhb <- plot_imp |>
  dplyr::filter(nhrace == "NHB")

# localized 
plot_imp_local <- plot_imp |>
  dplyr::filter(stage.f == "Localized")

# regional/distant 
plot_imp_regdis <- plot_imp |>
  dplyr::filter(stage.f %in% c("Regional", "Distant"))


# age stratification
## < 65


## All-cause mortality plots 
all_cause_spline_plots_lt65 <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_lt65,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline_lt65 <- wrap_plots(
  all_cause_spline_plots_lt65,
  nrow = 1
)

print(combined_all_cause_spline_lt65)

ggsave(
  filename = "pfas_all_cause_spline_lt65.png",
  plot = combined_all_cause_spline_lt65,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots

cvd_spline_plots_lt65 <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_lt65,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline_lt65 <- wrap_plots(
  cvd_spline_plots_lt65,
  nrow = 1
)

print(combined_cvd_spline_lt65)

ggsave(
  filename = "pfas_cvd_spline_lt65.png",
  plot = combined_cvd_spline_lt65,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer 

prostate_spline_plots_lt65 <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_lt65,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_prostate_spline_lt65 <- wrap_plots(
  prostate_spline_plots_lt65,
  nrow = 1
)

print(combined_prostate_spline_lt65)

ggsave(
  filename = "pfas_prostate_spline_lt65.png",
  plot = combined_prostate_spline_lt65,
  width = 15,
  height = 4.5,
  dpi = 300
)




# 65 + 

## All-cause mortality plots 
all_cause_spline_plots_65plus <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_65plus,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline_65plus <- wrap_plots(
  all_cause_spline_plots_65plus,
  nrow = 1
)

print(combined_all_cause_spline_65plus)

ggsave(
  filename = "pfas_all_cause_spline_65plus.png",
  plot = combined_all_cause_spline_65plus,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots 

cvd_spline_plots_65plus <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_65plus,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline_65plus <- wrap_plots(
  cvd_spline_plots_65plus,
  nrow = 1
)

print(combined_cvd_spline_65plus)

ggsave(
  filename = "pfas_cvd_spline_65plus.png",
  plot = combined_cvd_spline_65plus,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer 

prostate_spline_plots_65plus <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_65plus,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_prostate_spline_65plus <- wrap_plots(
  prostate_spline_plots_65plus,
  nrow = 1
)

print(combined_prostate_spline_65plus)

ggsave(
  filename = "pfas_prostate_spline_65plus.png",
  plot = combined_prostate_spline_65plus,
  width = 15,
  height = 4.5,
  dpi = 300
)


# race stratification 

## non-hispanic white 

## All-cause mortality plots 
all_cause_spline_plots_nhw <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nhw,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline_nhw <- wrap_plots(
  all_cause_spline_plots_nhw,
  nrow = 1
)

print(combined_all_cause_spline_nhw)

ggsave(
  filename = "pfas_all_cause_spline_nhw.png",
  plot = combined_all_cause_spline_nhw,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots 

cvd_spline_plots_nhw <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nhw,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline_nhw <- wrap_plots(
  cvd_spline_plots_nhw,
  nrow = 1
)

print(combined_cvd_spline_nhw)

ggsave(
  filename = "pfas_cvd_spline_nhw.png",
  plot = combined_cvd_spline_nhw,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer 

prostate_spline_plots_nhw <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nhw,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_prostate_spline_nhw <- wrap_plots(
  prostate_spline_plots_nhw,
  nrow = 1
)

print(combined_prostate_spline_nhw)

ggsave(
  filename = "pfas_prostate_spline_nhw.png",
  plot = combined_prostate_spline_nhw,
  width = 15,
  height = 4.5,
  dpi = 300
)


## non-hispanic black 

## All-cause mortality plots 
all_cause_spline_plots_nhb <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nhb,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline_nhb <- wrap_plots(
  all_cause_spline_plots_nhb,
  nrow = 1
)

print(combined_all_cause_spline_nhb)

ggsave(
  filename = "pfas_all_cause_spline_nhb.png",
  plot = combined_all_cause_spline_nhb,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots 

cvd_spline_plots_nhb <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nhb,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline_nhb <- wrap_plots(
  cvd_spline_plots_nhb,
  nrow = 1
)

print(combined_cvd_spline_nhb)

ggsave(
  filename = "pfas_cvd_spline_nhb.png",
  plot = combined_cvd_spline_nhb,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer 

pfas_spline_vars_nhb <- pfas_spline_vars[
  names(pfas_spline_vars) != "pfna_dx"
]

prostate_spline_plots_nhb <- lapply(
  names(pfas_spline_vars_nhb),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nhb,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars_nhb[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_prostate_spline_nhb <- patchwork::wrap_plots(
  prostate_spline_plots_nhb,
  nrow = 1
)



ggsave(
  filename = "pfas_prostate_spline_nhb.png",
  plot = combined_prostate_spline_nhb,
  width = 15,
  height = 4.5,
  dpi = 300
)




# non-white

## All-cause mortality plots 
all_cause_spline_plots_nonwhite <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nonwhite,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline_nonwhite <- wrap_plots(
  all_cause_spline_plots_nonwhite,
  nrow = 1
)

print(combined_all_cause_spline_nonwhite)

ggsave(
  filename = "pfas_all_cause_spline_nonwhite.png",
  plot = combined_all_cause_spline_nonwhite,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots 

cvd_spline_plots_nonwhite <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nonwhite,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline_nonwhite <- wrap_plots(
  cvd_spline_plots_nonwhite,
  nrow = 1
)

print(combined_cvd_spline_nonwhite)

ggsave(
  filename = "pfas_cvd_spline_nonwhite.png",
  plot = combined_cvd_spline_nonwhite,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer 

prostate_spline_plots_nonwhite <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_nonwhite,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3, # changing to 2 due to sparse cells error
      reference = "zero"
    )
  }
)

combined_prostate_spline_nonwhite <- wrap_plots(
  prostate_spline_plots_nonwhite,
  nrow = 1
)

print(combined_prostate_spline_nonwhite)

ggsave(
  filename = "pfas_prostate_spline_nonwhite.png",
  plot = combined_prostate_spline_nonwhite,
  width = 15,
  height = 4.5,
  dpi = 300
)



#### stage 

#### localized 


## All-cause mortality plots 
all_cause_spline_plots_local <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_local,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline_local <- wrap_plots(
  all_cause_spline_plots_local,
  nrow = 1
)

print(combined_all_cause_spline_local)

ggsave(
  filename = "pfas_all_cause_spline_local.png",
  plot = combined_all_cause_spline_local,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots 

cvd_spline_plots_local <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_local,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline_local <- wrap_plots(
  cvd_spline_plots_local,
  nrow = 1
)

print(combined_cvd_spline_local)

ggsave(
  filename = "pfas_cvd_spline_local.png",
  plot = combined_cvd_spline_local,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer 

prostate_spline_plots_local <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_local,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3, 
      reference = "zero"
    )
  }
)

combined_prostate_spline_local <- wrap_plots(
  prostate_spline_plots_local,
  nrow = 1
)

print(combined_prostate_spline_local)

ggsave(
  filename = "pfas_prostate_spline_local.png",
  plot = combined_prostate_spline_local,
  width = 15,
  height = 4.5,
  dpi = 300
)


#### regional/distant


## All-cause mortality plots 
all_cause_spline_plots_regdis <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_regdis,
      exposure = exposure,
      outcome = "all_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "All-cause mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_all_cause_spline_regdis <- wrap_plots(
  all_cause_spline_plots_regdis,
  nrow = 1
)

print(combined_all_cause_spline_regdis)

ggsave(
  filename = "pfas_all_cause_spline_regdis.png",
  plot = combined_all_cause_spline_regdis,
  width = 15,
  height = 4.5,
  dpi = 300
)


## cvd plots 

cvd_spline_plots_regdis <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_regdis,
      exposure = exposure,
      outcome = "cvd_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Cardiovascular mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_cvd_spline_regdis <- wrap_plots(
  cvd_spline_plots_regdis,
  nrow = 1
)

print(combined_cvd_spline_regdis)

ggsave(
  filename = "pfas_cvd_spline_regdis.png",
  plot = combined_cvd_spline_regdis,
  width = 15,
  height = 4.5,
  dpi = 300
)

#### prostate cancer 

prostate_spline_plots_regdis <- lapply(
  names(pfas_spline_vars),
  function(exposure) {
    
    get_pfas_spline_plot(
      data = plot_imp_regdis,
      exposure = exposure,
      outcome = "prostate_death",
      x_label = pfas_spline_vars[[exposure]],
      outcome_title = "Prostate cancer mortality",
      spline_df = 3,
      reference = "zero"
    )
  }
)

combined_prostate_spline_regdis <- wrap_plots(
  prostate_spline_plots_regdis,
  nrow = 1
)

print(combined_prostate_spline_regdis)

ggsave(
  filename = "pfas_prostate_spline_regdis.png",
  plot = combined_prostate_spline_regdis,
  width = 15,
  height = 4.5,
  dpi = 300
)

####################################### LIKELIHOOD RATIO TEST ################################################
get_spline_lrt <- function(
    data,
    exposure,
    outcome,
    spline_df = 3
) {
  
  # Keep identical observations for both models
  model_data <- data |>
    filter(
      !is.na(.data[[exposure]]),
      !is.na(.data[[outcome]]),
      !is.na(tstart),
      !is.na(tstop),
      tstop > tstart
    )
  
  # Linear model
  linear_formula <- as.formula(
    paste0(
      "Surv(tstart, tstop, ", outcome, ") ~ ",
      exposure, " + ",
      "age_cat.f + dxyr_cat.f + ",
      "nhrace + married.f + insur.f + stage.f + ",
      "txsurg.f + txrad.f + psa.f + ",
      "nses_quin + no2_new + pm25_new + ndvi27_final"
    )
  )
  
  # Spline model
  spline_formula <- as.formula(
    paste0(
      "Surv(tstart, tstop, ", outcome, ") ~ ",
      "pspline(", exposure, ", df = ", spline_df, ") + ",
      "age_cat.f + dxyr_cat.f + ",
      "nhrace + married.f + insur.f + stage.f + ",
      "txsurg.f + txrad.f + psa.f + ",
      "nses_quin + no2_new + pm25_new + ndvi27_final"
    )
  )
  
  # Fit models
  linear_fit <- coxph(
    linear_formula,
    data = model_data,
    x = TRUE,
    model = TRUE
  )
  
  spline_fit <- coxph(
    spline_formula,
    data = model_data,
    x = TRUE,
    model = TRUE
  )
  
  # Likelihood ratio statistic
  lr_linear <- as.numeric(logLik(linear_fit))
  lr_spline <- as.numeric(logLik(spline_fit))
  
  LRT_stat <- -2 * (lr_linear - lr_spline)
  
  # Difference in model degrees of freedom
  df_diff <- attr(logLik(spline_fit), "df") -
    attr(logLik(linear_fit), "df")
  
  p_nonlinearity <- pchisq(
    LRT_stat,
    df = df_diff,
    lower.tail = FALSE
  )
  
  tibble(
    exposure = exposure,
    outcome = outcome,
    LRT = LRT_stat,
    df = df_diff,
    p_nonlinearity = p_nonlinearity
  )
}

spline_lrt_results <- expand_grid(
  exposure = names(pfas_spline_vars),
  outcome_name = names(outcomes)
) |>
  mutate(
    outcome = outcomes[outcome_name]
  ) |>
  pmap_dfr(
    function(exposure, outcome_name, outcome) {
      
      get_spline_lrt(
        data = plot_imp,
        exposure = exposure,
        outcome = outcome,
        spline_df = 3
      ) |>
        mutate(
          outcome_name = outcome_name
        )
    }
  )

spline_lrt_results


write.csv(
  spline_lrt_results,
  file =
    "~/spline_lrt_results.csv",
  row.names = FALSE
)



