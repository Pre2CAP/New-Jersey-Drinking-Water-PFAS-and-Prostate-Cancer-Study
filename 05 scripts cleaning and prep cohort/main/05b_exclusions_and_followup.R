#################################################################
# R Code - exclusions and define follow up
## registry data

# Megan Shanahan
# Date: 6/25/2026

######################################################################

library(tidyverse)
library(dplyr)
library(haven)
library(khsverse)
library(rifttable)
library(rsvg)
library(DiagrammeRsvg)

# load data with death data and follow up variables
load("~/R/registry/data/njscr_final_1jul26.Rda")

nj_final <- nj_final |>
  rename(patient_study_id = study_id) |>
  mutate(patient_study_id = as.numeric(patient_study_id))

## read in NJ reshx data
nj_res <- read_sas("~/R/registry/nj_res_hx/reshist_ext.sas7bdat") 
length(unique(nj_res$patient_study_id)) 

nj_check <- nj_res %>%
  mutate(
    AddrSt = trimws(as.character(AddrSt))
  ) %>%
  group_by(patient_study_id) %>%
  summarise(
    n_rows = n(),
    all_nj = all(AddrSt == "NJ", na.rm = FALSE),
    any_non_nj = any(AddrSt != "NJ" | is.na(AddrSt)),
    states_seen = paste(sort(unique(AddrSt)), collapse = ", "),
    .groups = "drop"
  )

all_nj <- nj_check |>
  filter(all_nj == TRUE)

# loading in long data for other follow up variable 
load("~/R/registry/nj_res_hx/df_long_timevarying.Rda")
long <- person_qtr_exp_nses_psa |>
  select(Addr_StudyID, patient_study_id, dx_date, fup_end, overlap_start, overlap_end, year, pfoa, nses_quin, qtr_start, qtr_end) |>
  mutate(
    season = factor(
      quarter(qtr_start),
      levels = 1:4,
      labels = c("Q1", "Q2", "Q3", "Q4")
    ))


long_exclusions_fup <- left_join(nj_final, long, by="patient_study_id")
length(unique(long_exclusions_fup$patient_study_id))

##################################### SURV ACTIVE MOS VARIABLE ####################################

# creating follow up variable from surv month active follow up

long_exclusions_fup <- long_exclusions_fup |>
  mutate(
    end_fup_date = as.Date(newdxdt %m+% months(surv_mos_active_fup))
  )


long_exclusions_fup <- long_exclusions_fup |>
  mutate(fup_final = case_when(is.na(end_fup_date) ~ newenddt, 
                               TRUE ~ end_fup_date))

summary(long_exclusions_fup$fup_final)


############################## creating follow up time variable #########################################

df_long_fup_check <- long_exclusions_fup |>
  mutate(new_interval_start = pmax(overlap_start, newdxdt),
         new_interval_end = pmin(overlap_end, fup_final ), 
         tstart = as.numeric(new_interval_start - newdxdt),
         tstop = as.numeric(new_interval_end - newdxdt)) 

df_long_start_check <- df_long_fup_check |>
  select(c(patient_study_id, Addr_StudyID, newdxdt, overlap_start, tstart))


####################### CARRYING FORWARD/BACKWARD PFAS ###########################
# for participants with missing PFAS either at diagnosis or over follow up, carry forward or backwards any other PFAS values

df_long_fup_check_pfoa <- df_long_fup_check |>
  arrange(patient_study_id, season, year) |>
  group_by(patient_study_id) |>
  fill(pfoa, .direction = "down") |>
  fill(pfoa, .direction = "up") |>
  ungroup()


###################################### EXCLUSIONS #####################################

# creating a version that is one line per person for exclusion flow chart
df_exclusions_onerow <- df_long_fup_check_pfoa |>
  mutate(
    diagnosis_interval_distance = case_when(
      is.na(newdxdt) |
        is.na(overlap_start) |
        is.na(overlap_end) ~ Inf,
      
      overlap_start <= newdxdt &
        newdxdt <= overlap_end ~ 0,
      
      newdxdt < overlap_start ~
        as.numeric(overlap_start - newdxdt),
      
      newdxdt > overlap_end ~
        as.numeric(newdxdt - overlap_end)
    )
  ) |>
  group_by(patient_study_id) |>
  slice_min(
    order_by = diagnosis_interval_distance,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  mutate(
    diagnosis_interval_distance =
      if_else(is.infinite(diagnosis_interval_distance),
              NA_real_,
              diagnosis_interval_distance)
  )


length(df_exclusions_onerow$patient_study_id) 

# exclusion critera: continuous NJ resident, not missing diagnosis date, diagnosis after jan 1 2006, diagnosis not on autopsy, available PFAS data, 


# identifying patients missing address at diagnosis
# filter to just first row (at diagnosis)
df_diagnosis <- df_long_fup_check_pfoa |>
  filter(
    overlap_start <= newdxdt,
    newdxdt <= overlap_end
  ) %>%
  distinct(patient_study_id, .keep_all = TRUE)

length(unique(df_diagnosis$patient_study_id)) 

no_address_at_diagnosis <- df_long_fup_check_pfoa |>
  filter(!(patient_study_id %in% df_diagnosis$patient_study_id))

length(unique(no_address_at_diagnosis$patient_study_id))


# creating variable for difference between diagnosis date and newenddt
df_exclusions_onerow <- df_exclusions_onerow |>
  mutate(diag_death_dif = as.numeric(newenddt - newdxdt))

## missing pfoa
missing_pfoa <- df_exclusions_onerow |>
  filter(patient_study_id %in% all_nj$patient_study_id) |>
  filter(newdxdt >= as.Date("2006-01-01")) |>
  filter(diag_death_dif > 1) |>
  filter(is.na(pfoa))



# missing pfoa over follow up
missing_pfoa_long <- long |>
  filter(patient_study_id %in% missing_pfoa$patient_study_id)


missing_pfoa_partial <- missing_pfoa_long |>
  filter(!(is.na(pfoa)))

# negative follow up check
address_before <- df_exclusions_onerow |>
  filter(patient_study_id %in% all_nj$patient_study_id) |>
  filter(newdxdt >= as.Date("2006-01-01")) |>
  filter(diag_death_dif > 1) |>
  filter(!(is.na(pfoa))) |>
  filter((patient_study_id %in% neg_fup$patient_study_id)) |>
  select(patient_study_id, Addr_StudyID, dx_date, overlap_start, overlap_end, pfoa, fup_final)
  

exclusion_criteria <- bind_rows(
  tibble(
    left = "Initial registry population",
    right = "No continuous NJ residency",
    filter = list(expr(patient_study_id %in% all_nj$patient_study_id))
  ), 
  tibble(
    left = "NJ only population",
    right = "Diagnosis before 2006",
    filter = list(expr(newdxdt >= as.Date("2006-01-01")))
  ), 
  tibble(
    left = "Diagnosed after Jan 1 2006",
    right = "Diagnosed on autopsy/no follow up time",
    filter = list(expr(diag_death_dif > 1 & (is.na(surv_mos_active_fup) | surv_mos_active_fup != 0)))
  ),
  tibble(
    left = "Diagnosed after 2006, not on autopsy",
    right = "Missing PFAS at diagnosis",
    filter = list(expr((!is.na(pfoa)))) 
  ),
  tibble(
    left = "Population with PFAS data",
    right = "Address history only available before diagnosis",
    filter = list(expr((!(patient_study_id %in% neg_fup$patient_study_id)))) 
  ),
  tibble(
    left = "Population after address history exclusions",
    right = "Address not available at diagnosis",
    filter = list(expr((!(patient_study_id %in% no_address_at_diagnosis$patient_study_id)))) 
  ),
  tibble(
    left = "Population with address at diagnosis",
    right = "Missing nSES",
    filter = list(expr(!(is.na(nses_quin)))) 
  ),
  tibble(
    left = "Analytical population",
    right = "",
    filter = list(expr(TRUE))
  )
)

#make exclusions for the flowchart
# Print criteria
exclusion_criteria

# Make exclusions
exclusions_flowchart <- make_exclusions(
  criteria = exclusion_criteria,
  data = df_exclusions_onerow
)

# Pull put the analytical population
flowchart <- exclusions_flowchart |>
  pluck("included") |>
  last()

# Draw and save a flowchart (and display it *outside* the code chunk)
exclusions_flowchart |>
  exclusion_flowchart(width = 5) |>
  export_svg() |>
  charToRaw() |>
  rsvg_svg(file = "study_pop_exclusion.svg")

#make exclusions for the analytical dataset

df_long_after_exclusions <- long_exclusions_fup |>
  filter(patient_study_id %in% flowchart$patient_study_id)

length(unique(df_long_after_exclusions$patient_study_id))

# carry forward and backward PFOA in this dataset

df_long_after_exclusions <- df_long_after_exclusions |>
  arrange(patient_study_id, season, year) |>
  group_by(patient_study_id) |>
  fill(pfoa, .direction = "down") |>
  fill(pfoa, .direction = "up") |>
  ungroup()

################### DEATH EVENT ##################################
df_long_death_first <- df_long_after_exclusions  %>%
  mutate(
    all_death = case_when(dth==1 & overlap_start <= fup_final & fup_final <= overlap_end ~ 1,
                          TRUE ~ 0),
    cvd_death = case_when(dth==1 & overlap_start <= fup_final &  fup_final <= overlap_end & dthcause_label == "Diseases of the circulatory system" ~ 1,
                          TRUE ~ 0),
    prostate_death  = case_when(dth==1 & overlap_start <= fup_final & fup_final <= overlap_end & dthcause_label == "Prostate Cancer" ~ 1, 
                                TRUE ~ 0)
  ) %>%
  arrange(patient_study_id, overlap_start) %>%
  group_by(patient_study_id) %>%
  filter(cumsum(lag(all_death, default = 0)) == 0) %>%  # drops rows after death
  ungroup() 

length(unique(df_long_death_first$patient_study_id)) 




######### Creating follow up time variables #################
df_long_death <- df_long_death_first |>
  mutate(new_interval_start = pmax(overlap_start, newdxdt),
         new_interval_end = pmin(overlap_end + 1, fup_final +1 ),
         tstart = as.numeric(new_interval_start - newdxdt),
         tstop = as.numeric(new_interval_end - newdxdt)) |>
  filter(tstop > tstart)


####################### saving data #####################################

save(df_long_death, file = "~/R/registry/df_long_death.Rda")
