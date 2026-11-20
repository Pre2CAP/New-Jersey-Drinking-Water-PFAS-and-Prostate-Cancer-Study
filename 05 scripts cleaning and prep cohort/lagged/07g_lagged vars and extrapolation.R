#################################################################
# R Code - exclusions and then extrapolating air pollution and ndvi measures

# Hari Iyer
# Date: 10/28/2026

# Goal: Prepare long database with merged PFAS modeled data to geocodes
#################################################################

library(sf)
library(tidyverse)
library(haven)
library(dplyr)
library(stringr)
library(readr)
library(tidycensus)
library(lubridate)
library(tigris)
library(broom)


## load all datasets
load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/expanded_pfas_lag_clean_f.Rda") ## PFAS and study time
load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/expanded_nses_lag.Rda") ## nses
load("/extra/home/hi97/R/registry/nj_res_hx/person_qtr_exp.Rda")

names(person_qtr_exp)
head(person_qtr_exp)

## check for duplicate quarters and IDs
check_dups <- person_qtr_exp %>%
  add_count(
    patient_study_id,
    qtr_start,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1) %>%
  arrange(patient_study_id, qtr_start, interval_start) ## 144597 rows

## NDVI, PM2.5, NO2

person_qtr_exp <- person_qtr_exp %>%
  dplyr::select(Addr_StudyID, patient_study_id, year, pm25_val, no2_val, ndvi27_val, ndvi12_val, qtr_start, qtr_end, overlap_start, overlap_end)


tv_vars <- c("pm25_val", "no2_val", "ndvi27_val", "ndvi12_val")

weighted_nonmissing <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(x[keep], w[keep])
}

expanded_poll_clean <- person_qtr_exp %>%
  rename(quarter_start = qtr_start) %>%
  mutate(
    row_start = as.Date(overlap_start),
    row_end   = as.Date(overlap_end),
    exposure_days = as.numeric(row_end - row_start) + 1
  ) %>%
  group_by(patient_study_id, quarter_start) %>%
  summarise(
    quarter_end = first(qtr_end),
    year = first(year),
    quarter     = quarter(quarter_start) - 1L,
    
    row_start = min(row_start, na.rm = TRUE),
    row_end   = max(row_end, na.rm = TRUE),
    
    across(
      all_of(tv_vars),
      ~ weighted_nonmissing(.x, exposure_days)
    ),
    
    n_source_rows = n(),
    addresses_used = paste(
      unique(na.omit(as.character(Addr_StudyID))),
      collapse = "; "
    ),
    
    .groups = "drop"
  )

summary(expanded_poll_clean) 

expanded_poll_clean <- expanded_poll_clean %>%
  mutate(
    season = factor(
      quarter(quarter_start),
      levels = 1:4,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  )
  


length(unique(expanded_pfas_lag_clean_f$patient_study_id))


length(unique(expanded_nses_lag$patient_study_id)) 
length(unique(expanded_poll_clean$patient_study_id)) 

id_lookup <- expanded_pfas_lag_clean_f %>%
  distinct(patient_study_id)

exp_poll_subset <- expanded_poll_clean %>%
  semi_join(id_lookup, by = "patient_study_id") 

## need to merge

df_pca <- expanded_pfas_lag_clean_f %>%
  dplyr::select(patient_study_id, dx_date, fup_end) %>%
  distinct()


exp_poll_subset <- exp_poll_subset %>%
  left_join(df_pca, by = "patient_study_id")

## use 2010-2016 to extrapolate forward 

# just select id, year, season and air pollution
df_airpollution_yearly <- exp_poll_subset |>
  select(patient_study_id, year, pm25_val, no2_val, fup_end) |>
  group_by(patient_study_id, year, fup_end) |>
  summarise(
    pm25_val = mean(pm25_val, na.rm = TRUE),
    no2_val  = mean(no2_val, na.rm = TRUE),
    .groups = "drop"
  )

# predicting for 2017 to 2024

# fitting a safe model so it returns NA for people who do not have enough data to fit a linear trend
safe_lm <- function(df, outcome) {
  df2 <- df |>
    filter(!is.na(.data[[outcome]]), !is.na(year))
  
  if (nrow(df2) < 2 || n_distinct(df2$year) < 2) {
    return(NULL)
  }
  
  lm(reformulate("year", outcome), data = df2)
}


# Fit person-specific linear models using 2010-2016
person_models <- df_airpollution_yearly |>
  filter(year >= 2010, year <= 2016) |>
  group_by(patient_study_id) |>
  nest() |>
  mutate(
    pm25_model = map(data, ~ safe_lm(.x, "pm25_val")),
    no2_model  = map(data, ~ safe_lm(.x, "no2_val"))
  )




# Predict 2017-2024 values
future_dat <- expand_grid(
  patient_study_id = unique(person_models$patient_study_id),
  year = 2017:2024
)

future_pred <- person_models |>
  select(patient_study_id, pm25_model, no2_model) |>
  right_join(future_dat, by = "patient_study_id") |>
  mutate(
    pm25_extrap = map2_dbl(
      pm25_model, year,
      ~ if (is.null(.x)) NA_real_ else predict(.x, newdata = data.frame(year = .y))
    ),
    no2_extrap = map2_dbl(
      no2_model, year,
      ~ if (is.null(.x)) NA_real_ else predict(.x, newdata = data.frame(year = .y))
    )
    ) |>
  mutate(
    pm25_extrap = pmax(pm25_extrap, 0),
    no2_extrap  = pmax(no2_extrap, 0)
  ) |>
  select(patient_study_id, year, pm25_extrap, no2_extrap)


summary(future_pred$pm25_extrap)


# merging new air pollution with analytic df and creating new variable that combines them
df_long_analytic_airpollution <- exp_poll_subset |>
  left_join(future_pred, by=c("patient_study_id", "year")) |>
  mutate(no2_new = case_when(is.na(no2_val) ~ no2_extrap,
                             TRUE ~ no2_val),
         pm25_new = case_when(is.na(pm25_val) ~ pm25_extrap,
                              TRUE ~ pm25_val)) |>
  arrange(patient_study_id, year) |>
  group_by(patient_study_id) |>
  fill(no2_new, pm25_new, .direction = "down") |>
  ungroup()

df_long_analytic_airpollution <- df_long_analytic_airpollution %>%
  filter(year >= 2006) 

summary(df_long_analytic_airpollution$no2_new)

length(unique(df_long_analytic_airpollution$patient_study_id))

# impute the median air pollution measures for each year

yearly_medians <- df_long_analytic_airpollution |>
  filter(year >= 2017) |>
  group_by(year) |>
  summarise(
    median_pm25 = median(pm25_new, na.rm = TRUE),
    median_no2  = median(no2_new, na.rm = TRUE),
    .groups = "drop"
  )

df_long_analytic_airpollution <- df_long_analytic_airpollution |>
  left_join(yearly_medians, by = "year") |>
  mutate(
    pm25_new = coalesce(pm25_new, median_pm25),
    no2_new  = coalesce(no2_new, median_no2)
  ) |>
  select(-median_pm25, -median_no2) 


summary(df_long_analytic_airpollution$pm25_new)


summary(df_long_analytic_airpollution$no2_new)


## LAGGED PM and NO2


lag_quarters <- 20 # five-year lag
poll_lag_vars <- c(
  "pm25_new",
  "no2_new"
)

exposure_lookup <- df_long_analytic_airpollution %>%
  select(
    patient_study_id,
    exposure_quarter = quarter_start,
    all_of(poll_lag_vars)
  ) %>%
  distinct() %>%
  rename_with(
    ~ paste0(.x, "_lag"),
    all_of(poll_lag_vars)
  )

df_long_analytic_airpollution <- df_long_analytic_airpollution %>%
  mutate(
    exposure_quarter =
      quarter_start %m-% months(3 * lag_quarters)
  ) %>%
  left_join(
    exposure_lookup,
    by = c("patient_study_id", "exposure_quarter"),
    relationship = "many-to-one"
  )

summary(df_long_analytic_airpollution)

save(df_long_analytic_airpollution, file="/extra/home/hi97/R/registry/nj_res_hx/PFAS study/df_lag_analytic_airpollution.Rda")

#################### NDVI ################################################

table(exp_poll_subset$quarter)

quarter_year_ndvi <- exp_poll_subset |>
  mutate(
    quarter_factor = factor(
      quarter,
      levels = 0:3,
      labels = c("Q1","Q2","Q3","Q4")
    )
  ) |>
  filter(!is.na(ndvi27_val)) |>
  group_by(year, quarter_factor) |>
  summarize(
    mean_ndvi = mean(ndvi27_val, na.rm = TRUE),
    .groups = "drop"
  )

summary(quarter_year_ndvi)

fit_window <- function(start_year, end_year) {
  quarter_year_ndvi |>
    filter(year >= start_year,
           year <= end_year) |>
    group_by(quarter_factor) |>
    do(
      glance(lm(mean_ndvi ~ year, data = .))
    ) |>
    mutate(
      start = start_year,
      end = end_year
    ) |>
    select(
      quarter_factor,
      start,
      end,
      r.squared,
      adj.r.squared
    )
}

bind_rows(
  fit_window(2010, 2018),
  fit_window(2012, 2018),
  fit_window(2014, 2018),
  fit_window(2010, 2016),
  fit_window(2012, 2016)
)

get_slopes <- function(start_year, end_year) {
  quarter_year_ndvi |>
    filter(year >= start_year,
           year <= end_year) |>
    group_by(quarter_factor) |>
    summarize(
      slope = coef(lm(mean_ndvi ~ year))[2],
      .groups = "drop"
    ) |>
    mutate(start = start_year,
           end = end_year)
}

bind_rows(
  get_slopes(2010, 2018),
  get_slopes(2012, 2018),
  get_slopes(2012, 2016)
)

# linear extrapolation from 2012 to 2018 to predict

ndvi_model_data <- exp_poll_subset |>
  filter(
    !is.na(ndvi27_val),
    year >= 2012,
    year <= 2018
  ) |>
  group_by(patient_study_id, season, year) |>
  summarize(
    mean_ndvi = mean(ndvi27_val, na.rm = TRUE),
    .groups = "drop"
  )

summary(ndvi_model_data)

eligible_2018 <- ndvi_model_data |>
  filter(year == 2018) |>
  distinct(patient_study_id, season)

ndvi_model_data_2018 <- ndvi_model_data |>
  semi_join(
    eligible_2018,
    by = c("patient_study_id", "season")
  )



participant_qtr_models <- ndvi_model_data_2018 |>
  group_by(patient_study_id, season) |>
  filter(n_distinct(year) >= 3) |>
  nest() |>
  mutate(
    model = map(data, ~ lm(mean_ndvi ~ year, data = .x))
  )


# create predictions years
pred_years <- tibble(year = 2019:2024)

# predict NDVI
ndvi_predictions_2019_2024 <- participant_qtr_models |>
  mutate(
    pred = map(model, predict, newdata = pred_years)
  ) |>
  select(patient_study_id, season, pred) |>
  unnest(pred) |>
  group_by(patient_study_id, season) |>
  mutate(year = 2019:2024) |>
  ungroup() |>
  rename(pred_ndvi = pred)


## seeing which participants did not get an NDVI predicted 

ndvi_data_coverage <- ndvi_model_data |>
  group_by(patient_study_id, season) |>
  summarize(
    n_years = n_distinct(year),
    years_available = paste(sort(unique(year)), collapse = ", "),
    .groups = "drop"
  )

not_enough_data_participant_qtr <- ndvi_data_coverage |>
  filter(n_years < 3)

length(unique(not_enough_data_participant_qtr$patient_study_id))

eligible_missing <- not_enough_data_participant_qtr |>
  filter(patient_study_id %in% eligible_2018$patient_study_id)
length(unique(eligible_missing$patient_study_id))

# using last observation per quarter carried forward for the 244 IDs that don't have enough data

last_ndvi <- ndvi_model_data_2018 |>
  group_by(patient_study_id, season) |>
  slice_max(year, n = 1, with_ties = FALSE) |>
  select(patient_study_id, season, last_year = year,
         last_ndvi = mean_ndvi)

future_years <- tibble(year = 2019:2024)

locf_predictions <- last_ndvi |>
  tidyr::crossing(future_years) |>
  rename(pred_ndvi = last_ndvi) |>
  filter(patient_study_id %in% not_enough_data_participant_qtr$patient_study_id) |>
  select(-(c(last_year))) |>
  rename(ndvi27_new = pred_ndvi) |>
  anti_join(
    ndvi_predictions_2019_2024 |>
      distinct(patient_study_id, year, season),
    by = c("patient_study_id", "year", "season")
  )


# bind together 2 predictions
ndvi_predictions_2019_2024 <- ndvi_predictions_2019_2024 |>
  rename(ndvi27_new = pred_ndvi )


ndvi <- rbind(ndvi_predictions_2019_2024, locf_predictions)
length(unique(ndvi$patient_study_id)) 


summary(ndvi)

## truncate values <1 to 0, values >=1 to 1

ndvi <- ndvi %>%
  mutate(ndvi27_new = case_when(
    ndvi27_new > 1 ~ 1,
    ndvi27_new  <= 0 ~ 0,
    TRUE ~ ndvi27_new
  ))


# merging ndvi with df_long_analytic

exp_ndvi_subset <- exp_poll_subset %>%
  dplyr::select(patient_study_id, year, season, ndvi27_val, quarter_start, row_start, row_end, dx_date, fup_end) %>%
  left_join(ndvi, by = c("patient_study_id", "year", "season")) %>%
  mutate(ndvi27_final = case_when(is.na(ndvi27_val) ~ ndvi27_new,
                                  TRUE ~ ndvi27_val))

summary(exp_ndvi_subset$ndvi27_final)

exp_ndvi_subset <- exp_ndvi_subset %>%
  mutate(ndvi27_final = case_when(
    ndvi27_final < 0 ~ 0,
    TRUE ~ ndvi27_final
  ))

summary(exp_ndvi_subset$ndvi27_final) 

missingNDVI <- exp_ndvi_subset |>
  filter(is.na(ndvi27_final))

# impute closest year/season combination 

exp_ndvi_subset <- exp_ndvi_subset |>
  arrange(patient_study_id, season, year) |>
  group_by(patient_study_id, season) |>
  fill(ndvi27_final, .direction = "down") |>
  fill(ndvi27_final, .direction = "up") |>
  ungroup()

summary(exp_ndvi_subset$ndvi27_final)


## LAG


lag_quarters <- 20 # five-year lag
poll_lag_vars <- "ndvi27_final"

rm(exposure_lookup)

exposure_lookup <- exp_ndvi_subset %>%
  select(
    patient_study_id,
    exposure_quarter = quarter_start,
    all_of(poll_lag_vars)
  ) %>%
  distinct() %>%
  rename_with(
    ~ paste0(.x, "_lag"),
    all_of(poll_lag_vars)
  )

exp_ndvi_subset <- exp_ndvi_subset %>%
  mutate(
    exposure_quarter =
      quarter_start %m-% months(3 * lag_quarters)
  ) %>%
  left_join(
    exposure_lookup,
    by = c("patient_study_id", "exposure_quarter"),
    relationship = "many-to-one"
  )

summary(exp_ndvi_subset)


exp_ndvi_subset <- exp_ndvi_subset %>%
  dplyr::select(patient_study_id, quarter_start, year, row_start, row_end, ndvi27_final, ndvi27_final_lag)


### MERGE AP AND NDVI LAGGED VALUES

load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/df_lag_analytic_airpollution.Rda")


names(df_long_analytic_airpollution)
names(exp_ndvi_subset)


ap_dups <- df_long_analytic_airpollution %>%
  add_count(
    patient_study_id,
    quarter_start,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1) %>%
  arrange(patient_study_id, quarter_start, row_start)



df_long_analytic_airpollution <- df_long_analytic_airpollution %>%
  dplyr::select(patient_study_id, year, quarter_start, quarter_end, row_start, row_end, pm25_new_lag, no2_new_lag)

tv_vars <- c("pm25_new_lag", "no2_new_lag")

weighted_nonmissing <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(x[keep], w[keep])
}

df_ap_summary <- df_long_analytic_airpollution %>%
  mutate(
    quarter_start = as.Date(quarter_start),
    row_start = as.Date(row_start),
    row_end   = as.Date(row_end),
    exposure_days = as.numeric(row_end - row_start) + 1
  ) %>%
  group_by(patient_study_id, quarter_start) %>%
  summarise(
    quarter_end = first(quarter_end),
    year = first(year),
    quarter = quarter(first(quarter_start)) - 1L,
    
    row_start = min(row_start, na.rm = TRUE),
    row_end   = max(row_end, na.rm = TRUE),
    
    across(
      all_of(tv_vars),
      ~ weighted_nonmissing(.x, exposure_days)
    ),
    
    .groups = "drop"
  )

df_ap_summary %>%
  count(patient_study_id, quarter_start) %>%
  filter(n > 1)


summary(df_ap_summary) 

names(exp_ndvi_subset)

tv_vars <- c("ndvi27_final_lag")

weighted_nonmissing <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(x[keep], w[keep])
}

df_ndvi_summary <- exp_ndvi_subset %>%
  mutate(
    quarter_start = as.Date(quarter_start),
    row_start = as.Date(row_start),
    row_end   = as.Date(row_end),
    exposure_days = as.numeric(row_end - row_start) + 1
  ) %>%
  group_by(patient_study_id, quarter_start) %>%
  summarise(
    year = first(year),
    quarter = quarter(first(quarter_start)) - 1L,
    
    row_start = min(row_start, na.rm = TRUE),
    row_end   = max(row_end, na.rm = TRUE),
    
    across(
      all_of(tv_vars),
      ~ weighted_nonmissing(.x, exposure_days)
    ),
    
    .groups = "drop"
  )

df_ndvi_summary %>%
  count(patient_study_id, quarter_start) %>%
  filter(n > 1)

names(df_ndvi_summary)
names(df_ap_summary)

df_long_analytic_extrap <- df_ap_summary %>%
  left_join(df_ndvi_summary, by = c("patient_study_id","quarter_start","quarter","year","row_start","row_end")) %>%
  mutate(
    season = factor(
      quarter(quarter_start),
      levels = 1:4,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  )

summary(df_long_analytic_extrap)

### merge PFAS data

names(df_long_analytic_extrap)
names(expanded_pfas_lag_clean_f)

df_long_analytic_extrap <- expanded_pfas_lag_clean_f %>%
  left_join(df_long_analytic_extrap, by = c("patient_study_id","quarter_start","quarter_end","row_start","row_end","quarter","year")) 

summary(df_long_analytic_extrap)

length(unique(df_long_analytic_extrap$patient_study_id))

ndvi_miss <- df_long_analytic_extrap %>%
  filter(is.na(ndvi27_final_lag)) 
length(unique(ndvi_miss$patient_study_id)) 

ap_miss <- df_long_analytic_extrap %>%
  filter(is.na(pm25_new_lag))
length(unique(ap_miss$patient_study_id)) 

## check who is missing pfas

pfas_vars <- c("pfoa_exp_lag", "pfos_exp_lag", "pfna_exp_lag")  # replace with your PFAS variables

person_missing <- df_long_analytic_extrap %>%
  group_by(patient_study_id) %>%
  summarise(
    pfas_completely_missing =
      all(if_all(all_of(pfas_vars), is.na)),
    
    pm25_completely_missing = 
      all(is.na(pm25_new_lag)),
    
    no2_completely_missing =
      all(is.na(no2_new_lag)),
    
    ndvi_completely_missing =
      all(is.na(ndvi27_final_lag)),
    
    .groups = "drop"
  )

summary(person_missing)


df_analytic_clean <- df_long_analytic_extrap %>%
  semi_join(
    person_missing %>%
      filter(
        !pfas_completely_missing,
        !no2_completely_missing,
        !ndvi_completely_missing
      ),
    by = "patient_study_id"
  )

summary(df_analytic_clean)
length(unique(df_analytic_clean$patient_study_id))



df_analytic_clean %>%
  group_by(patient_study_id) %>%
  summarise(
    n_rows = n(),
    n_missing_pfoa = sum(is.na(pfoa_exp_lag)),
    n_observed_pfoa = sum(!is.na(pfoa_exp_lag)),
    .groups = "drop"
  ) %>%
  summarise(
    people_all_pfoa_missing = sum(n_observed_pfoa == 0),
    people_partially_missing = sum(
      n_missing_pfoa > 0 & n_observed_pfoa > 0
    ),
    people_with_observed_pfoa = sum(n_observed_pfoa > 0)
  )

# seasonal imputation

## carry forward

impute_vars <- c("pfoa_exp_lag","pfos_exp_lag","pfna_exp_lag","pm25_new_lag","no2_new_lag")

df_analytic_clean <- df_analytic_clean %>%
  arrange(patient_study_id, quarter_start) %>%
  group_by(patient_study_id) %>%
  tidyr::fill(
    all_of(impute_vars),
    .direction = "down"
  ) %>%
  ungroup()

summary(df_analytic_clean)

## impute for PFOA, PFOS, PFNA, and NDVI using seasonal medians

impute_vars <- c("pfoa_exp_lag","pfos_exp_lag","pfna_exp_lag","ndvi27_final_lag")

df_analytic_clean <- df_analytic_clean %>%
  mutate(calendar_quarter = lubridate::quarter(quarter_start))

seasonal_medians <- df_analytic_clean %>%
  group_by(calendar_quarter) %>%
  summarise(
    across(
      all_of(impute_vars),
      ~ median(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

df_analytic_clean <- df_analytic_clean %>%
  left_join(seasonal_medians, by = "calendar_quarter",
            suffix = c("", "_seasonal")) 

for (v in impute_vars) {
  seasonal_v <- paste0(v, "_seasonal")
  
  df_analytic_clean[[v]] <- dplyr::coalesce(
    df_analytic_clean[[v]],
    df_analytic_clean[[seasonal_v]]
  )
}

df_analytic_clean <- df_analytic_clean %>%
  select(-ends_with("_seasonal"))

summary(df_analytic_clean) 


## Just for pollution

covars <- c("pm25_new_lag","no2_new_lag")

person_flags <- df_analytic_clean %>%
  group_by(patient_study_id) %>%
  summarise(
    pm25_completely_missing = all(is.na(pm25_new_lag)),
    no2_completely_missing = all(is.na(no2_new_lag)),
    .groups = "drop"
  )

df_imp <- df_analytic_clean %>%
  semi_join(
    person_flags %>%
      filter(
        !no2_completely_missing,
        !pm25_completely_missing
      ),
    by = "patient_study_id"
  ) %>%
  arrange(patient_study_id, quarter_start) %>%
  group_by(patient_study_id) %>%
  tidyr::fill(
    all_of(covars),
    .direction = "down"
  ) %>%
  ungroup()


imputation_medians <- df_imp %>%
  summarise(
    across(
      all_of(covars),
      ~ median(.x, na.rm = TRUE)
    )
  ) %>%
  as.list()

df_analytic_imputed <- df_analytic_clean %>%
  mutate(
    across(
      all_of(covars),
      ~ tidyr::replace_na(
        .x,
        imputation_medians[[cur_column()]]
      )
    )
  )


df_analytic_imputed <- df_analytic_imputed %>%
  dplyr::select(-c(season))

summary(df_analytic_imputed)
length(unique(df_analytic_imputed$patient_study_id))
length(unique(df_analytic_clean$patient_study_id))
summary(df_analytic_clean)

# end of air pollution and ndvi extrapolution - df_long_analytic_extrap is df with new air pollution and ndvi measures


########################### just saving time varying variables #############################


save(df_analytic_imputed, file = "/extra/home/hi97/R/registry/nj_res_hx/PFAS study/df_analytic_imputed.Rda")

