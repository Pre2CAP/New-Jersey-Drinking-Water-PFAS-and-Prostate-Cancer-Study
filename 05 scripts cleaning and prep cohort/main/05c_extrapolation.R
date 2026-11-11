#################################################################
# R Code - exclusions and then extrapolating air pollution and ndvi measures

# Megan Shanahan
# Date: 06/02/2026

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

Sys.setenv(HOME = "/extra/home/ms3283")

load("~/R/registry/nj_res_hx/df_long_timevarying.Rda") 
load("~/R/registry/df_long_death.Rda") # 1947606

df_long_analytic <- person_qtr_exp_nses_psa |>
  filter(patient_study_id %in% df_long_death$patient_study_id) |>
  filter(year >= 2004) |>
  mutate(
    season = factor(
      quarter(qtr_start),
      levels = 1:4,
      labels = c("Q1", "Q2", "Q3", "Q4")
    ))


#######################################################
# Examining trends in PM2.5 and NO2 for extrapolation
######################################################

library(dplyr)
library(ggplot2)

annual_means <- df_long_death %>%
  group_by(year) %>%
  summarize(
    mean_pm25 = mean(pm25_val, na.rm = TRUE),
    mean_no2  = mean(no2_val, na.rm = TRUE),
    n = n()
  )

ggplot(annual_means, aes(year, mean_pm25)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_bw()

ggplot(annual_means, aes(year, mean_no2)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_bw()

## use 2010-2016 to extrapolate forward 

# just select id, year, season and air pollution
df_airpollution_yearly <- df_long_analytic |>
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


# checking to make sure everyone with air pollution data in 2016 has a prediction
airpollution_eligible_2016 <- df_airpollution_yearly |>
  filter(year==2016)

# merging new air pollution with analytic df and creating new variable that combines them
df_long_analytic_airpollution <- df_long_analytic |>
  left_join(future_pred, by=c("patient_study_id", "year")) |>
  mutate(no2_new = case_when(is.na(no2_val) ~ no2_extrap,
                             TRUE ~ no2_val),
         pm25_new = case_when(is.na(pm25_val) ~ pm25_extrap,
                              TRUE ~ pm25_val)) |>
  arrange(patient_study_id, year) |>
  group_by(patient_study_id) |>
  fill(no2_new, pm25_new, .direction = "down") |>
  ungroup()

# impute the median air pollution measures for each year for IDs completely missing air pollution data

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


#################### NDVI ################################################

quarter_year_ndvi <- df_long_analytic |>
  mutate(
    quarter_factor = factor(
      quarter(qtr_start),
      levels = 1:4,
      labels = c("Q1","Q2","Q3","Q4")
    )
  ) |>
  filter(!is.na(ndvi27_val)) |>
  group_by(year, quarter_factor) |>
  summarize(
    mean_ndvi = mean(ndvi27_val, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  quarter_year_ndvi,
  aes(year, mean_ndvi)
) +
  geom_point() +
  geom_line() +
  facet_wrap(~quarter_factor) +
  theme_bw()

library(broom)

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

# predicting code
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(broom)


ndvi_model_data <- df_long_analytic |>
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

# using last observation per quarter carried forward for the IDs that don't have enough data

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

# merging ndvi with df_long_analytic

# there are multiple lines per patient, per year, per season because some have multiple addresses
# don't rbind with before 2018, just merge original with after 2018 and then create a new combined ndvi variable

df_long_analytic_extrap <- df_long_analytic_airpollution |>
  left_join(ndvi, by=c("patient_study_id", "year", "season")) |>
  mutate(ndvi27_final = case_when(is.na(ndvi27_val) ~ ndvi27_new,
                                  TRUE ~ ndvi27_val)) 

missingNDVI <- df_long_analytic_extrap |>
  filter(is.na(ndvi27_final))

# these IDs just have random years missing NDVI - impute closest year/season combination

df_long_analytic_extrap <- df_long_analytic_extrap |>
  arrange(patient_study_id, season, year) |>
  group_by(patient_study_id, season) |>
  fill(ndvi27_final, .direction = "down") |>
  fill(ndvi27_final, .direction = "up") |>
  ungroup()


# end of air pollution and ndvi extrapolution - df_long_analytic_extrap is df with new air pollution and ndvi measures


########################### just saving time varying variables #############################

df_long_analytic_extrap <- df_long_analytic_extrap |>
  select(c(Addr_StudyID, patient_study_id, year, season, pfoa, pfos, pfna, PWSID, nses_quin, psa_perc, no2_new, pm25_new, ndvi27_final ))


save(df_long_analytic_extrap, file = "~/R/registry/df_long_analytic_extrap.Rda")

