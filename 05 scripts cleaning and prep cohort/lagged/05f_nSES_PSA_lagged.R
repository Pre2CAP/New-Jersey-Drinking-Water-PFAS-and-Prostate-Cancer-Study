#################################################################
# R Code - add nSES and county level PSA to long time varying data

# Hari Iyer
# Date: Updated 10/28/2026

# Goal: Prepare long database with merged PFAS modeled data to geocodes
#################################################################

library(tidyverse)
library(dplyr)

########################################################################## 
# Merging person_qtr_exp data with nses data
##########################################################################

load("~/R/registry/nj_res_hx/nj_res_nses_f.Rda")
load("/extra/home/hi97/R/registry/nj_res_hx/person_qtr_exp.Rda")
load("/extra/home/hi97/R/registry/nj_res_hx/PFAS study/expanded_pfas_lag_clean_f.Rda")

## First, expand the nSES file to capture quarters

nj_res_nses <- nj_res_nses_f %>%
  dplyr::select(Addr_StudyID, patient_study_id, FromDate, ToDate, mdinc, mdhval, edubach, percfamint, pohse, percwhite, 
                percblack, percforeign, unemp,
                mdinc10, mdhval10, edubach10, percfamint10, pohse10, percwhite10, 
                percblack10, percforeign10, unemp10) %>%
  mutate(
    across(
      -(c(patient_study_id, Addr_StudyID, FromDate, ToDate)),
      ~ as.numeric(.x)
    )
  )


study_start <- as.Date("2006-01-01") ## these are the date ranges we need
study_end <- as.Date("2024-10-31")


expanded_nses <- nj_res_nses %>%
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


duplicate_rows <- expanded_nses %>%
  add_count(
    patient_study_id,
    quarter_start,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1) %>%
  arrange(patient_study_id, quarter_start, row_start) 


nses_vars <- c("mdinc", "mdhval", "edubach", "percfamint", "pohse", "percwhite", 
               "percblack", "percforeign", "unemp",
               "mdinc10", "mdhval10", "edubach10", "percfamint10", "pohse10", "percwhite10", 
               "percblack10", "percforeign10", "unemp10")

weighted_nonmissing <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(x[keep], w[keep])
}

expanded_nses_clean <- expanded_nses %>%
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
      all_of(nses_vars),
      ~ weighted_nonmissing(.x, exposure_days)
    ),
    
    n_source_rows = n(),
    addresses_used = paste(
      unique(na.omit(as.character(Addr_StudyID))),
      collapse = "; "
    ),
    
    .groups = "drop"
  )

summary(expanded_nses_clean) 


## carry forward if missing after
expanded_nses_clean <- expanded_nses_clean %>%
  arrange(patient_study_id, quarter_start) %>%
  group_by(patient_study_id) %>%
  fill(
    mdinc, mdhval, edubach, percfamint, pohse, percwhite, 
    percblack, percforeign, unemp,
    mdinc10, mdhval10, edubach10, percfamint10, pohse10, percwhite10, 
    percblack10, percforeign10, unemp10,
    .direction = "down"
  ) %>%
  ungroup()
summary(expanded_nses_clean) 

## replace 2006-2009 with 2000 census, >=2010 with 2010 ACS
nses_base <- c(
  "mdinc",
  "mdhval",
  "edubach",
  "percfamint",
  "pohse",
  "percwhite",
  "percblack",
  "percforeign",
  "unemp"
)

expanded_nses_clean_tv <- expanded_nses_clean

for (v in nses_base) {
  
  v10 <- paste0(v, "10")
  
  base_value <- as.numeric(
    as.character(expanded_nses_clean_tv[[v]])
  )
  
  value_10 <- as.numeric(
    as.character(expanded_nses_clean_tv[[v10]])
  )
  
  expanded_nses_clean_tv[[v]] <- dplyr::case_when(
    expanded_nses_clean_tv$year < 2010 ~
      dplyr::coalesce(base_value, value_10),
    
    expanded_nses_clean_tv$year >= 2010 ~
      dplyr::coalesce(value_10, base_value),
    
    TRUE ~ dplyr::coalesce(base_value, value_10)
  )
}

expanded_nses_clean_tv <- expanded_nses_clean_tv %>%
  select(-all_of(paste0(nses_base, "10")))

summary(expanded_nses_clean_tv)

## calculate zscore


expanded_nses_clean_tv <- expanded_nses_clean_tv %>%
  mutate(
    across(
      -(c(patient_study_id, quarter_start, quarter_end, year, quarter, row_start, row_end, n_source_rows, addresses_used)),
      ~ as.numeric(scale(.x))
    )
  ) %>%
  mutate(
    zscore = mdinc + mdhval + edubach +
      percfamint + pohse +
      percwhite + percblack +
      percforeign + unemp
  )

summary(expanded_nses_clean_tv$zscore) 

## Apply lag

lag_quarters <- 20 # five-year lag
nses_lag_vars <- c(
  "mdinc",
  "mdhval",
  "edubach",
  "percfamint",
  "pohse",
  "percwhite",
  "percblack",
  "percforeign",
  "unemp",
  "zscore"
)

exposure_lookup <- expanded_nses_clean_tv %>%
  select(
    patient_study_id,
    exposure_quarter = quarter_start,
    all_of(nses_lag_vars)
  ) %>%
  distinct() %>%
  rename_with(
    ~ paste0(.x, "_lag"),
    all_of(nses_lag_vars)
  )

expanded_nses_lag <- expanded_nses_clean_tv %>%
  mutate(
    exposure_quarter =
      quarter_start %m-% months(3 * lag_quarters)
  ) %>%
  left_join(
    exposure_lookup,
    by = c("patient_study_id", "exposure_quarter"),
    relationship = "many-to-one"
  )

summary(expanded_nses_lag)

save(expanded_nses_lag, file="/extra/home/hi97/R/registry/nj_res_hx/PFAS study/expanded_nses_lag.Rda")
