#################################################################
# R Code - merging all data together, creating exposure quintiles, and creating follow up variables

# Megan Shanahan
# Date: 06/02/2026

# Goal: Prepare long database with merged PFAS modeled data to geocodes
#################################################################

library(tidyverse)
library(haven)
library(dplyr)
library(stringr)
library(readr)
library(lubridate)

load("~/R/registry/df_long_death.Rda")
load("~/R/registry/df_long_analytic_extrap.Rda")

problem_IDs <- df_long_death |>
  filter(!(patient_study_id %in% df_long_analytic_extrap$patient_study_id))
length(unique(problem_IDs$patient_study_id))


df_long_death <- df_long_death |>
  select(-c(pfoa, nses_quin, qtr_start, qtr_end))

summary(df_long_analytic_extrap$pfoa)
class(df_long_analytic_extrap$pfoa)

# carrying pfoa down and up
df_long_analytic_extrap <- df_long_analytic_extrap |>
  arrange(patient_study_id, year, season) |>
  group_by(patient_study_id) |>
  fill(pfoa, .direction = "down") |>
  fill(pfoa, .direction = "up") |>
  fill(pfos, .direction = "down") |>
  fill(pfos, .direction = "up") |>
  fill(pfna, .direction = "down") |>
  fill(pfna, .direction = "up") |>
  fill(nses_quin, .direction = "down") |>
  fill(nses_quin, .direction = "up") |>
  ungroup()


df_long_final <- left_join(df_long_death, df_long_analytic_extrap, by = c("patient_study_id", "year", "season", "Addr_StudyID"))

df_long_final <- df_long_final |>
  mutate(new_interval_start = pmax(overlap_start, newdxdt),
         new_interval_end = pmin(overlap_end, fup_final ), 
         tstart = as.numeric(new_interval_start - newdxdt),
         tstop = as.numeric(new_interval_end - newdxdt)) 

max(df_long_final$overlap_end)
max(df_long_final$fup_final)

length(unique(df_long_final$patient_study_id)) 

df_long_final <- df_long_final |>
  filter(tstop >= tstart)

length(unique(df_long_final$patient_study_id)) 

summary(df_long_final$pfoa)
summary(df_long_final$pfos)
summary(df_long_final$pfna)


# creating exposure quintile variables 

# first log transformaing
df_long_final <- df_long_final |>
  mutate(pfoa_log = log(pfoa),
         pfos_log = log(pfos),
         pfna_log = log(pfna))

df_long_final <- df_long_final |>
  mutate(pfoa_quin = as.factor(ntile(pfoa_log, 5)),
         pfos_quin = as.factor(ntile(pfos_log, 5)),
         pfna_quin = as.factor(ntile(pfna_log, 5)))

# cleaning covariate variables 

df_long_final <- df_long_final |>
  mutate(stage.f = factor(case_when(advstage3== 1 ~ "Localized",
                                    advstage3== 2 ~ "Regional", 
                                    advstage3== 3 ~ "Distant",
                                    TRUE ~ NA),
                          levels=c("Localized", "Regional", "Distant" )),
         mets = case_when(stage.f %in% c("Regional", "Distant") ~ 1,
                          stage.f=="Localized" ~ 0, 
                          TRUE ~ NA),
         married.f = factor(case_when(married==0 ~ "No",
                                      married==1 ~ "Yes",
                                      married==2 ~ "Missing"),
                            levels=c("No", "Yes", "Missing")),
         insur.f = factor(case_when(insur==1 ~ "Uninsured",
                                    insur==0 ~ "Private",
                                    insur==2 ~ "Medicaid",
                                    insur==3 ~ "Medicare",
                                    insur==4 ~ "Other government",
                                    TRUE ~ NA),
                          levels=c("Uninsured", "Private", "Medicaid", "Medicare", "Other government" )),
         age_3cat.f = factor(case_when(dxage <55 ~ "<55",
                                      dxage >=55 & dxage <70 ~ "55-69",
                                      dxage>= 70 ~ "70+"),
                            levels=c("<55", "55-69", "70+")),
         age_cat.f = factor(case_when(age_cat==1 ~ "<50",
                                      age_cat==2 ~ "50-60",
                                      age_cat==3 ~ "60-70",
                                      age_cat==4 ~ "70-80",
                                      age_cat==5 ~ "80+"),
                            levels=c("<50", "50-60", "60-70", "70-80", "80+" )),
         dxyr_group.f = factor(case_when(dxyr_group==1 ~ "2000-2005",
                                         dxyr_group==2 ~ "2005-2010",
                                         dxyr_group==3 ~ "2010-2015"),
                               levels=c("2000-2005", "2005-2010", "2010-2015")),
         dxyr_cat.f = factor(case_when(dxyr < 2010 ~ "2006-2010",
                                       dxyr >= 2010 & dxyr < 2014 ~ "2010-2014",
                                       dxyr >= 2014 & dxyr < 2020 ~ "2014 -2017"),
                             levels=c("2006-2010", "2010-2014", "2014 -2017")),
         txchemo.f = factor(case_when(txchemo_val==1 ~ "Yes",
                                      txchemo_val==0 ~ "No",
                                      TRUE ~ NA),
                            levels=c("No", "Yes")),
         txhorm.f = factor(case_when(txhorm_val==1 ~ "Yes",
                                     txhorm_val==0 ~ "No",
                                     TRUE ~ NA),
                           levels=c("No", "Yes")),
         txrad.f = factor(case_when(txrad_val==1 ~ "Yes",
                                    txrad_val==0 ~ "No",
                                    TRUE ~ NA),
                          levels=c("No", "Yes")),
         txsurg.f = factor(case_when(txsurg_any==1 ~ "Yes",
                                     txsurg_any==0 ~ "No"),
                           levels=c("No", "Yes")),
         grade_new.f=factor(case_when(grade_new==1 ~ "<=4",
                                      grade_new==2 ~ "5-6",
                                      grade_new==3 ~ "7+",
                                      TRUE ~ NA),
                            levels=c("<=4", "5-6", "7+" )),
         gleason.f=factor(case_when(glscss6 %in% c("002", "003", "004", "005", "006") ~ "<7",
                                    glscss6 %in% c("007") ~ "7",
                                    glscss6 %in% c("008", "009", "010") ~ "8-10",
                                    TRUE ~ NA),
                          levels=c("<7", "7", "8-10"  )),
         grade_combo.f = factor(case_when(grade_new.f %in% c("<=4", "5-6") | gleason.f %in% c("<7") ~ "<7",
                                          grade_new.f %in% c("7+") | gleason.f %in% c("7", "8-10" ) ~ "7+"),
                                levels=c("<7", "7+" )),
         afr_amer.f = factor(case_when(afr_amer==1 ~ "Black",
                                afr_amer==0 ~ "White",
                                afr_amer==2 ~ "Missing"),
                             levels=c("Black", "White", "Missing"  )),
         hisp= case_when(!is.na(nhia) & nhia > 0 ~ 1,
                         TRUE ~ 0),
         hisp.f = case_when(hisp==1 ~ "Hispanic",
                            hisp==0 ~ "Non-Hispanic"),
         nhrace = factor(case_when(hisp.f=="Non-Hispanic" & afr_amer.f=="White" ~ "NHW",
                                   hisp.f=="Non-Hispanic" & afr_amer.f=="Black" ~ "NHB",
                                   hisp.f=="Hispanic" ~ "Hispanic",
                                   TRUE ~ "NHAAPI"),
                         levels=c("NHW", "NHB", "NHAAPI", "Hispanic" )),
         psa_div = psa_lab/10,
         psa.f = factor(case_when(psa_div < 10 ~ "<10",
                           psa_div >= 10 & psa_div < 20 ~ "10-20",
                           psa_div >= 20 ~ "20+",
                           TRUE ~ NA),
                        levels=c("<10", "10-20", "20+" ))
         )


summary(df_long_final$psa_perc)
psaperc_missing <- df_long_final |>
  filter(is.na(psa_perc)) |>
  dplyr::select(c(patient_study_id, psa_perc, year)) |>
  distinct(patient_study_id, .keep_all = TRUE)
 save(psaperc_missing, file = "~/R/registry/psaperc_missing.Rda")

summary(df_long_final$nses_quin) 

df_long_final <- df_long_final |>
  mutate(psa_perc_cat = ntile(psa_perc, 4),
         psa_perc_cat.f = case_when(psa_perc_cat==1 ~ "Q1",
                                    psa_perc_cat==2 ~ "Q2",
                                    psa_perc_cat==3 ~ "Q3",
                                    psa_perc_cat==4 ~ "Q4",
                                    TRUE ~ "Missing"))

table(df_long_final$psa_perc_cat.f, useNA="always")
summary(df_long_final$nses_quin)

save(df_long_final, file = "~/R/registry/df_long_final.Rda")
