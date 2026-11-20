#################################################################
# R Code - merging all data together, creating exposure quintiles, and creating follow up variables

# Hari Iyer, Megan Shanahan
# Date: 10/29/2026

# Goal: Prepare long database with merged PFAS modeled data to geocodes
#################################################################

library(tidyverse)
library(haven)
library(dplyr)
library(stringr)
library(readr)
library(lubridate)

## check that we can just merge our new data with the original imputed datasets
load("/extra/home/hi97/R/registry/nj_res_hx/long_list_imp.Rda")
load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/df_analytic_imputed.Rda")
head(long_imp_list[[1]])
summary(long_imp_list[[1]])
summary(long_imp_list[[2]])
summary(long_imp_list[[3]])
summary(long_imp_list[[4]])
summary(long_imp_list[[5]])

head(df_analytic_imputed)
names(df_analytic_imputed)

df_lag_for_merge <- df_analytic_imputed %>%
  mutate(
    season = factor(
      quarter(quarter_start),
      levels = 1:4,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  ) %>%
  dplyr::select(patient_study_id, quarter_start, season, year, pfoa_exp_lag, pfos_exp_lag, pfna_exp_lag, pm25_new_lag,
                no2_new_lag, ndvi27_final_lag)



# creating exposure quintile variables 
iqr_pfoa <- IQR(df_lag_for_merge$pfoa_exp_lag, na.rm = TRUE)
iqr_pfos <- IQR(df_lag_for_merge$pfos_exp_lag, na.rm = TRUE)
iqr_pfna <- IQR(df_lag_for_merge$pfna_exp_lag, na.rm = TRUE)

df_lag_for_merge <- df_lag_for_merge %>%
  mutate(
    pfoa_lag_log = if_else(pfoa_exp_lag > 0, log(pfoa_exp_lag), NA_real_),
    pfos_lag_log = if_else(pfos_exp_lag > 0, log(pfos_exp_lag), NA_real_),
    pfna_lag_log = if_else(pfna_exp_lag > 0, log(pfna_exp_lag), NA_real_),
    
    pfoa_lag_quin = ntile(pfoa_exp_lag, 5),
    pfos_lag_quin = ntile(pfos_exp_lag, 5),
    pfna_lag_quin = ntile(pfna_exp_lag, 5),
    
    pfoa_exp_lag_per_iqr = pfoa_exp_lag / iqr_pfoa,
    pfos_exp_lag_per_iqr = pfos_exp_lag / iqr_pfos,
    pfna_exp_lag_per_iqr = pfna_exp_lag / iqr_pfna
  ) %>%
  group_by(pfoa_lag_quin) %>%
  mutate(pfoa_exp_lag_t = median(pfoa_exp_lag, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(pfos_lag_quin) %>%
  mutate(pfos_exp_lag_t = median(pfos_exp_lag, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(pfna_lag_quin) %>%
  mutate(pfna_exp_lag_t = median(pfna_exp_lag, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(across(ends_with("_quin"), factor))

summary(df_lag_for_merge)

## nSES
load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/expanded_nses_lag.Rda")
summary(expanded_nses_lag)

names(df_analytic_imputed)

nses_for_merge <- expanded_nses_lag %>%
  dplyr::select(patient_study_id, year, quarter_start, zscore_lag)

df_lag_for_merge <- df_lag_for_merge %>%
  left_join(nses_for_merge , by = c("patient_study_id","year","quarter_start"))  |>
  arrange(patient_study_id, year, quarter_start) |>
  group_by(patient_study_id) |>
  fill(zscore_lag, .direction = "down") |>
  ungroup() 

summary(df_lag_for_merge)
  
## annual imputation for nSES
  
yearly_medians <- df_lag_for_merge |>
  group_by(year) |>
  summarise(
    median_zscore_lag = median(zscore_lag, na.rm = TRUE),
    .groups = "drop"
  )

df_lag_for_merge <- df_lag_for_merge %>%
  left_join(yearly_medians, by = "year") |>
  mutate(
    zscore_lag_new = coalesce(zscore_lag, median_zscore_lag),
  ) |>
  select(-median_zscore_lag) 
  
summary(df_lag_for_merge)
length(unique(df_lag_for_merge$patient_study_id))

df_lag_for_merge <- df_lag_for_merge %>%
  mutate(
    nses_lag_quin = ntile(zscore_lag_new, 5)
  )

## we should be able to merge on patient_study_id, season, and year

merge_keys <- c("patient_study_id", "season", "year")


## confirm no dups on keys

long_imp_list_lagged <- lapply(long_imp_list, function(dat) {
  
  dat %>%
    select(
      -any_of(c("pfoa_exp_lag", "pfos_exp_lag", "pfna_exp_lag"))
    ) %>%
    inner_join(
      df_lag_for_merge,
      by = merge_keys,
      relationship = "many-to-one"
    )
})

head(long_imp_list_lagged[[1]])
summary(long_imp_list_lagged[[1]])
names(long_imp_list_lagged[[1]])

table(long_imp_list_lagged[[1]]$all_death) 
table(long_imp_list_lagged[[1]]$prostate_death) 
table(long_imp_list_lagged[[1]]$cvd_death) 


## use inner join to just keep the ids/rows 

## check which rows would be removed before joining

## long_imp_list[[1]] %>%
##  anti_join(
##    df_lag_for_merge,
##    by = c("patient_study_id", "season", "year")
##  )


length(unique(long_imp_list_lagged[[1]]$patient_study_id)) 

table(long_imp_list_lagged[[1]]$nses_lag_quin)

save(long_imp_list_lagged, file = "/extra/home/hi97/R/registry/nj_res_hx/long_imp_list_lagged.Rda")
