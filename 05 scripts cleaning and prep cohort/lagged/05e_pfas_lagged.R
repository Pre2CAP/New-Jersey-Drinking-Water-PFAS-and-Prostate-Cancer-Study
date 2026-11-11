#################################################################
# R Code - Processing PFAS data for lagged analyses

# Hari Iyer
# Date: 10/26/2026

# Goal: Use multiple imputation to resolve missingness
# Add lagged exposures and restrict to those with complete data
# 5-year lags - so 2011 diagnosis year onward
#################################################################


library(tidyverse)
library(mice)
library(survival)
library(forcats)
library(broom)
library(dplyr)
library(ggplot2)
library(stringr)



library(sf)

load("~/R/registry/nj_res_hx/long_list_imp.Rda")

df_pca_vars <- long_imp_list[[1]] %>%
  dplyr::select(patient_study_id, Addr_StudyID, dxyr, dx_date) %>%
  dplyr::distinct()

load("/home/hi97/R/registry/nj_res_hx/nj_pfas_merged.Rda")

df_pca_lagged <- nj_pfas_merged %>%
  st_drop_geometry() %>%
  select(-any_of("centroid_geom")) %>%
  left_join(df_pca_vars, by = c("patient_study_id","Addr_StudyID"))

rm(long_imp_list)
rm(df_pca_vars)

## exclusions - need to keep only diagnosis year >= 2011 to allow up to 5 years of pre-diagnosis addresses

df_pca_lagged <- df_pca_lagged %>%
  filter(dxyr >= 2011)

length(unique(df_pca_lagged$patient_study_id)) 

summary(df_pca_lagged$ToDate)

library(purrr)

study_start <- as.Date("2006-01-01") 
study_end <- as.Date("2024-10-31")



id_cols <- c("patient_study_id", "Addr_StudyID")

pfas_cols <- names(df_pca_lagged)[
  grepl("_exp_\\d{4}_[0-3]$", names(df_pca_lagged))
]

# PFAS values in long format
pfas_long <- df_pca_lagged %>%
  select(all_of(id_cols), all_of(pfas_cols)) %>%
  distinct() %>%
  pivot_longer(
    cols = all_of(pfas_cols),
    names_to = c("analyte", "year", "quarter"),
    names_pattern = "^(.*)_exp_(\\d{4})_([0-3])$",
    values_to = "exposure"
  ) %>%
  mutate(
    year = as.integer(year),
    quarter = as.integer(quarter),
    quarter_start = make_date(
      year = year,
      month = 1 + 3 * quarter,
      day = 1
    )
  )

# Expand each FromDate-ToDate interval into quarters
expanded <- df_pca_lagged %>%
  select(all_of(id_cols), FromDate, ToDate) %>%
  distinct() %>%
  mutate(
    FromDate = as.Date(FromDate),
    ToDate   = as.Date(ToDate),
    
    interval_start = pmax(FromDate, study_start),
    interval_end   = pmin(ToDate, study_end)
  ) %>%
  filter(
    !is.na(interval_start),
    !is.na(interval_end),
    interval_start <= interval_end
  ) %>%
  rowwise() %>%
  mutate(
    quarter_start = list(
      seq.Date(
        from = floor_date(interval_start, "quarter"),
        to   = floor_date(interval_end, "quarter"),
        by   = "3 months"
      )
    )
  ) %>%
  ungroup() %>%
  unnest(quarter_start) %>%
  mutate(
    quarter_end = quarter_start %m+% months(3) - days(1),
    row_start   = pmax(interval_start, quarter_start),
    row_end     = pmin(interval_end, quarter_end),
    year        = year(quarter_start),
    quarter     = quarter(quarter_start) - 1L
  ) %>%
  filter(row_start <= row_end)

# Attach the correct quarterly PFAS values
expanded_pfas <- expanded %>%
  left_join(
    pfas_long,
    by = c(
      "patient_study_id",
      "Addr_StudyID",
      "quarter_start",
      "year",
      "quarter"
    )
  )

head(expanded_pfas)

expanded_pfas_wide <- expanded_pfas %>%
  pivot_wider(
    names_from = analyte,
    values_from = exposure,
    names_glue = "{analyte}_exp"
  )

summary(expanded_pfas_wide$pfoa_exp)
summary(expanded_pfas_wide$pfos_exp)
summary(expanded_pfas_wide$pfna_exp)

save(expanded_pfas_wide, file="/extra/home/hi97/R/registry/nj_res_hx/PFAS study/expanded_pfas_wide.Rda")

rm(df_pca_lagged, pfas_long, expanded, expanded_pfas)

## Clean duplicate rows for quarter_start and patient_study_id


pfas_vars <- names(expanded_pfas_wide) %>%
  str_subset("_exp$")

weighted_nonmissing <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(x[keep], w[keep])
}

expanded_pfas_clean <- expanded_pfas_wide %>%
  mutate(
    row_start = as.Date(row_start),
    row_end   = as.Date(row_end),
    exposure_days = as.numeric(row_end - row_start) + 1
  ) %>%
  group_by(patient_study_id, quarter_start) %>%
  summarise(
    quarter_end = first(quarter_end),
    year = first(year),
    quarter = first(quarter),
    
    row_start = min(row_start, na.rm = TRUE),
    row_end   = max(row_end, na.rm = TRUE),
    
    across(
      all_of(pfas_vars),
      ~ weighted_nonmissing(.x, exposure_days)
    ),
    
    n_source_rows = n(),
    addresses_used = paste(
      unique(na.omit(as.character(Addr_StudyID))),
      collapse = "; "
    ),
    
    .groups = "drop"
  )


load("~/R/registry/nj_res_hx/long_list_imp.Rda")

df_pca_vars <- long_imp_list[[1]] %>%
  dplyr::select(patient_study_id, Addr_StudyID, dxyr, dx_date, fup_end) %>%
  dplyr::distinct()


expanded_pfas_pca <- expanded_pfas_wide %>%
  left_join(df_pca_vars, by = c("patient_study_id", "Addr_StudyID")) %>%
  mutate(
    dx_date = as.Date(dx_date),
    fup_end = as.Date(fup_end),
    
    study_start = pmax(row_start, dx_date),
    study_end = pmin(row_end, fup_end)
  ) %>%
  
  mutate(
    time_start_days = as.numeric(study_start - dx_date),
    time_end_days   = as.numeric(study_end - dx_date),
    time_start_years = time_start_days / 365.25,
    time_end_years   = time_end_days / 365.25
  ) %>%
  filter( !is.na(study_start),
          !is.na(study_end),
          study_start <= study_end
  ) 

lag_quarters <- 20 # five-year lag

## duplicates checking

duplicate_rows <- expanded_pfas_wide %>%
  add_count(
    patient_study_id,
    quarter_start,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1) %>%
  arrange(patient_study_id, quarter_start, row_start)


exposure_lookup <- expanded_pfas_wide %>%
dplyr::select(
  patient_study_id,
  Addr_StudyID,
  exposure_quarter = quarter_start,
  ends_with("_exp")
) %>%
  distinct() %>%
  rename_with(
    ~ paste0(.x,"_lag"),
    ends_with("_exp")
  )

expanded_pfas_lag <- expanded_pfas_pca %>%
  mutate(
    exposure_quarter = 
      quarter_start %m-% months(3 * lag_quarters)
  ) %>%
  left_join(
    exposure_lookup,
    by = c("patient_study_id", "Addr_StudyID", "exposure_quarter"
  ),
relationship = "many-to-one"
)

## clean rows - take weighted average if multiple addresses within quarter, take the value that is non-missing if one is missing,
## otherwise NA

pfas_vars <- names(expanded_pfas_lag) %>%
  str_subset("_exp_lag$")

weighted_nonmissing <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(x[keep], w[keep])
}

expanded_pfas_lag_clean <- expanded_pfas_lag %>%
  mutate(
    row_start = as.Date(row_start),
    row_end   = as.Date(row_end),
    exposure_days = as.numeric(row_end - row_start) + 1
  ) %>%
  group_by(patient_study_id, quarter_start) %>%
  summarise(
    quarter_end = first(quarter_end),
    year = first(year),
    quarter = first(quarter),
    
    row_start = min(row_start, na.rm = TRUE),
    row_end   = max(row_end, na.rm = TRUE),
    
    across(
      all_of(pfas_vars),
      ~ weighted_nonmissing(.x, exposure_days)
    ),
    
    n_source_rows = n(),
    addresses_used = paste(
      unique(na.omit(as.character(Addr_StudyID))),
      collapse = "; "
    ),
    
    .groups = "drop"
  )

summary(expanded_pfas_lag_clean)

## carry forward if missing after
expanded_pfas_lag_clean <- expanded_pfas_lag_clean %>%
  arrange(patient_study_id, quarter_start) %>%
  group_by(patient_study_id) %>%
  fill(
    pfoa_exp_lag,
    pfos_exp_lag,
    pfna_exp_lag,
    .direction = "down"
  ) %>%
  ungroup()
summary(expanded_pfas_lag_clean) 

## merge diagnosis data

names(expanded_pfas_pca)

## just keep the time variables, id, and quarter variables to merge
pfas_pca <- expanded_pfas_pca %>%
  dplyr::select(patient_study_id, Addr_StudyID, quarter_start, quarter_end, dx_date, fup_end, time_start_days,
                time_end_days, time_start_years, time_end_years)

pfas_pca_dups <- pfas_pca %>%
  add_count(
    patient_study_id,
    quarter_start,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1) %>%
  arrange(patient_study_id, quarter_start, time_start_days)
  
pfas_pca <- pfas_pca %>%
  group_by(patient_study_id, quarter_start) %>%
  summarise(quarter_end = first(quarter_end),
            
            time_start_days = min(time_start_days, na.rm = TRUE),
            time_end_days = max(time_end_days, na.rm = TRUE),
            
            time_start_years = min(time_start_years, na.rm = TRUE), 
            time_end_years = max(time_end_years, na.rm=TRUE),
            
            n_address_rows = n(),
            n_addresses = n_distinct(Addr_StudyID)
            )
    


expanded_pfas_lag_clean_f <- expanded_pfas_lag_clean %>%
  left_join(pfas_pca, by = c("patient_study_id","quarter_start", "quarter_end"))



save(expanded_pfas_wide, file="/extra/home/hi97/R/registry/nj_res_hx/PFAS study/expanded_pfas_wide.Rda")



