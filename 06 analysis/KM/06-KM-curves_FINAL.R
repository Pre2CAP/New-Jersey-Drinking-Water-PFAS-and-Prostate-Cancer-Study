#############################################################################
# KM Curves for Paper
# Created by: Hari Iyer
# Date: 9/28/2026

# Goal: Present Kaplan-Meier curves for all-cause mortality, stratified
# by PFOA, PFOS, and PFNA categories/quintiles
# Present 
# - Log-rank p-values
# - deaths/population at risk
# - median survival times (if estimable)
############################################################################


library(tidyverse)
library(survival)
library(forcats)
library(broom)
library(ggplot2)
library(survminer)
library(patchwork)

load("/extra/home/hi97/R/registry/nj_res_hx/df_long_final.Rda") 
names(df_long_final)


# Define colors for quintiles
km_colors <- c("green", "blue", "darkgoldenrod", "orange", "red")



# Helper function to generate KM plots
make_km <- function(fit, indat, title, pval.y) {
  ggsurvplot(
    fit,
    data = indat,
    pval = TRUE,
    pval.coord = c(7.5, pval.y),
    pval.size = 4,
    risk.table = TRUE,
    risk.table.height = 0.25,
    xlab = "Years Since Diagnosis",
    ylab = "Survival Probability",
    xlim = c(0,17),
    ylim = c(0.2, 1.0),
    legend.labs = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    censor.size = 1,
    legend.title = "Quintile",
    palette = km_colors,
    title = title,
    size = 0.7,
    
    ggtheme = theme_bw() +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 11, hjust = 0.5)
    )
  )
}



# Helper function to generate KM plots
make_km_pfna <- function(fit, indat, title, pval.y) {
  ggsurvplot(
    fit,
    data = indat,
    pval = TRUE,
    pval.coord = c(7.5, pval.y),
    pval.size = 4,
    risk.table = TRUE,
    risk.table.height = 0.25,
    xlab = "Years Since Diagnosis",
    ylab = "Survival Probability",
    xlim = c(0,17),
    ylim = c(0.2, 1.0),
    legend.labs = c("<2 ng/L", "2-3.8 ng/L", "≥3.8 ng/L"),
    censor.size = 1,
    legend.title = "Categories",
    palette = km_colors,
    title = title,
    size = 0.7,
    
  ggtheme = theme_bw() +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 11, hjust = 0.5)
    )
  )
}


summary(df_long_final$tstart)
summary(df_long_final$tstop)


## KM requires one-row per person

df_km <- df_long_final %>%
  group_by(patient_study_id) %>%
  summarise(
    fup_time = (max(tstop, na.rm = TRUE))/365.24,
    all_death = max(all_death, na.rm = TRUE),
    pfoa_quin = first(pfoa_quin),
    pfos_quin = first(pfos_quin),
    pfna_cat  = first(pfna_cat),
    .groups = "drop"
  )

surv_allcause <- Surv(df_km$fup_time, df_km$all_death)


# --------- All-Cause KM Curves ---------


fit_pfoa_all <- survfit(
  surv_allcause ~ pfoa_quin,
  data = df_km
)

fit_pfos_all <- survfit(
  surv_allcause ~ pfos_quin,
  data = df_km
)

fit_pfna_all <- survfit(
  surv_allcause ~ pfna_cat,
  data = df_km
)




km_pfoa_all <- make_km(fit_pfoa_all, df_km, "PFOA (All-Cause Mortality)", 0.9)
km_pfos_all <- make_km(fit_pfos_all, df_km, "PFOS (All-Cause Mortality)", 0.9)
km_pfna_all <- make_km_pfna(fit_pfna_all, df_km,  "PFNA (All-Cause Mortality)", 0.9)



pfoa_panel <- km_pfoa_all$plot / km_pfoa_all$table +
  plot_layout(heights = c(3, 1))

pfos_panel <- km_pfos_all$plot / km_pfos_all$table +
  plot_layout(heights = c(3, 1))

pfna_panel <- km_pfna_all$plot / km_pfna_all$table +
  plot_layout(heights = c(3, 1))

km_all_combined <- pfoa_panel | pfos_panel | pfna_panel

km_all_combined

ggsave(
  filename = "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/output/km_all_panels.png",
  plot = km_all_combined,
  width = 20,  # adjust based on aspect ratio
  height = 12,
  dpi = 300
)

## Median survival times and 95% CI
# PFOA
pfoa_median <- summary(fit_pfoa_all)$table[, c(
  "median", "0.95LCL", "0.95UCL"
)] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group") %>%
  mutate(PFAS = "PFOA")

# PFOS
pfos_median <- summary(fit_pfos_all)$table[, c(
  "median", "0.95LCL", "0.95UCL"
)] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group") %>%
  mutate(PFAS = "PFOS")

# PFNA
pfna_median <- summary(fit_pfna_all)$table[, c(
  "median", "0.95LCL", "0.95UCL"
)] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group") %>%
  mutate(PFAS = "PFNA")

median_survival <- bind_rows(
  pfoa_median,
  pfos_median,
  pfna_median
) %>%
  rename(
    median_survival = median,
    lower_95CI = `0.95LCL`,
    upper_95CI = `0.95UCL`
  ) %>%
  select(PFAS, group, median_survival, lower_95CI, upper_95CI)

median_survival

