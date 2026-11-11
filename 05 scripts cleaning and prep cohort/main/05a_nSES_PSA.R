#################################################################
# R Code - add nSES and county level PSA to long time varying data

# Megan Shanahan
# Date: Updated 7/2/2026

# Goal: Prepare long database with merged PFAS modeled data to geocodes
#################################################################

library(tidyverse)
library(dplyr)

########################################################################## 
# Merging person_qtr_exp data with nses data
##########################################################################

load("~/R/registry/nj_res_hx/nj_res_nses_f.Rda")
load("~/R/registry/data/person_qtr_exp.Rda")

summary(person_qtr_exp$pfoa)

length(unique(person_qtr_exp$patient_study_id))

class(person_qtr_exp$year)
class(nj_res_nses_f$nses_year)

person_qtr_exp_nses_tomerge <- person_qtr_exp |>
  mutate(nses_year = case_when(year < 2010 ~ 2000, 
                               year >= 2010 ~ 2010))

names(person_qtr_exp_nses_tomerge)
names(nj_res_nses_f)


length(unique(person_qtr_exp$patient_study_id)) 
length(unique(nj_res_nses_f$patient_study_id))  

missing_nses <- person_qtr_exp |>
  filter(!(patient_study_id %in% nj_res_nses_f$patient_study_id ))
length(unique(missing_nses$patient_study_id))  

missing_ids <- nj_res_nses_f |>
  filter(!(patient_study_id %in% person_qtr_exp$patient_study_id ))
length(unique(missing_ids$patient_study_id)) 

# creating nses variable and cleaning nontime varying data
nses <- nj_res_nses_f |>
  select(c(patient_study_id, nses_year, Addr_StudyID, STATEFP, COUNTYFP, mdinc, mdhval, edubach, percfamint, pohse, percwhite, percblack, percforeign, unemp))

nses <- nses %>%
  mutate(across(c(mdinc, mdhval, edubach, percfamint, pohse, percwhite, percblack, percforeign, unemp), as.numeric))

nses <- nses %>%
  mutate(
    across(
      -(c(patient_study_id, nses_year, Addr_StudyID, STATEFP, COUNTYFP)),
      ~ as.numeric(scale(.x))
    )
  )

# if any of these are missing for 2010, carry forward the values for 2000
vars_to_fill <- c(
  "mdinc", "mdhval", "edubach", "percfamint", "pohse",
  "percwhite", "percblack", "percforeign", "unemp"
)

nses <- nses %>%
  mutate(
    across(
      all_of(vars_to_fill),
      ~ ifelse(.x == ".", 0, .x)
    )
  ) %>%
  mutate(
    across(
      all_of(vars_to_fill),
      as.numeric
    )
  ) %>%
  group_by(patient_study_id, Addr_StudyID) %>%
  mutate(
    across(
      all_of(vars_to_fill),
      ~ case_when(
        nses_year == 2010 & is.na(.x) ~ .x[nses_year == 2000 & !is.na(.x)][1],
        nses_year == 2000 & is.na(.x) ~ .x[nses_year == 2010 & !is.na(.x)][1],
        TRUE ~ .x
      )
    )
  ) %>%
  ungroup()


nses <- nses %>%
  mutate(
    zscore = mdinc + mdhval + edubach +
      percfamint + pohse +
      percwhite + percblack +
      percforeign + unemp
  )

summary(nses$zscore) 

# if z score is completely missing for one of the years, use the other year
nses <- nses %>%
  group_by(patient_study_id, Addr_StudyID) %>%
  mutate(
    zscore = coalesce(
      zscore,
      zscore[nses_year == 2000][1],
      zscore[nses_year == 2010][1]
    )
  ) %>%
  ungroup()

summary(nses$zscore)


#assigning quintiles of nses
nses <- nses |>
  mutate(nses_quin = factor(ntile(zscore, 5))) |>
  select(patient_study_id, Addr_StudyID, nses_year, nses_quin, STATEFP, COUNTYFP)


person_qtr_exp_nses <- left_join(person_qtr_exp_nses_tomerge, nses, by=c("patient_study_id", "nses_year", "Addr_StudyID"))

##################################### PSA ###########################################################

psa00 <- read_csv("~/R/registry/data/nj_psa00.csv") |>
  select(patient_study_id:psa08)

length(unique(psa00$patient_study_id)) 

psa10 <- read_csv("~/R/registry/data/nj_psa10.csv") |>
  select(patient_study_id:psa12)

length(unique(psa10$patient_study_id)) 

psa <- full_join(psa00, psa10, by=c("patient_study_id", "STATE", "county"))


psa_long <- psa |>
  pivot_longer(starts_with("psa"), names_to = "PSAYear", values_to = "psa_perc") |>
  mutate(
    year = str_extract(PSAYear, "\\d{2}$")
  ) |>
  select(-(c(PSAYear))) |>
  mutate(psa_year = case_when(year=="04" ~ 2004,
                              year=="06" ~ 2006,
                              year=="08" ~ 2008,
                              year=="10" ~ 2010,
                              year=="12" ~ 2012),
         STATEFP = STATE,
         COUNTYFP = substr(county, 3, 5)) |>
  select(-c(STATE, county, year))

length(unique(psa_long$patient_study_id)) 
class(psa_long$psa_year)
class(person_qtr_exp_nses$year)
length(unique(person_qtr_exp_nses$patient_study_id)) 

person_qtr_exp_nses_psa <- person_qtr_exp_nses |>
  mutate(
    psa_year = case_when(
      year <= 2005 ~ 2004,
      year %in% c(2006, 2007) ~ 2006,
      year %in% c(2008, 2009) ~ 2008,
      year %in% c(2010, 2011) ~ 2010,
      year >= 2012 ~ 2012
    )
  ) |>
  left_join(
    psa_long,
    by = c("patient_study_id", "psa_year", "STATEFP", "COUNTYFP")
  )

psaperc_check <- psa_long |>
  filter(!(patient_study_id %in% person_qtr_exp_nses_psa$patient_study_id))


length(unique(psa_long$patient_study_id)) 

length(unique(person_qtr_exp_nses_psa$patient_study_id)) 
length(unique(psa_long$patient_study_id)) 

load("~/R/registry/df_long_final.Rda")

length(unique(df_long_final$patient_study_id)) 

missing_psa_final <- df_long_final |>
  filter(!(patient_study_id %in% psa_long$patient_study_id)) |>
  select(patient_study_id) |>
  distinct(patient_study_id, .keep_all = TRUE) 

save(missing_psa_final, file = "~/R/registry/nj_res_hx/missing_psa_final.Rda")

person_qtr_exp_nses_psa <- person_qtr_exp_nses_psa |>
  arrange(patient_study_id, overlap_start) |>
  group_by(patient_study_id) |>
  fill(psa_perc, .direction = "downup") |>
  ungroup()

person_qtr_exp_nses_psa |>
  group_by(patient_study_id) |>
  summarise(
    any_missing = any(is.na(psa_perc)),
    any_available = any(!is.na(psa_perc)),
    .groups = "drop"
  ) |>
  filter(any_missing, any_available)


save(person_qtr_exp_nses_psa, file = "~/R/registry/nj_res_hx/df_long_timevarying.Rda")