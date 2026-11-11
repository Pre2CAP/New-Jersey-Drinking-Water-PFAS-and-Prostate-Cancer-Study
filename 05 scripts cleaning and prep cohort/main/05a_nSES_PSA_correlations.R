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

nj_res_nses_f <- nj_res_nses_f %>%
  mutate(
    STATEFP = case_when(
      nses_year == 2000 ~ STATEFP00,
      nses_year == 2010 ~ STATEFP10
    ),
    COUNTYFP = case_when(
      nses_year == 2000 ~ COUNTYFP00,
      nses_year == 2010 ~ COUNTYFP10
    )
  )

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

nses_df <- nses %>%
  arrange(patient_study_id, nses_year) %>%
  group_by(patient_study_id) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  dplyr::select(patient_study_id, zscore)

## merge with df_table

load("/extra/home/hi97/R/registry/nj_res_hx/df_table1_diagnosis.Rda")
names(df_table1_diagnosis)

df_tab1_corr <- df_table1_diagnosis %>%
  dplyr::select(patient_study_id, pfoa, pfoa_log, pfos, pfos_log, pfna, pfna_log, 
                ndvi27_final, pm25_new, no2_new) %>%
  left_join(nses_df, by = "patient_study_id")


### variable labels so they have better names
var_labels <- c(
  pfoa_log     = "PFOA",
  pfos_log     = "PFOS",
  pfna_log     = "PFNA",
  ndvi27_final = "NDVI\n(270)",
  pm25_new     = "PM₂.5",
  no2_new      = "NO₂",
  zscore       = "nSES"
)
# use zscore for nses so it's continuous
vars <- c("pfoa_log", "pfos_log", "pfna_log", "ndvi27_final", "pm25_new", "no2_new", "zscore")

corr_matrix <- cor(df_tab1_corr[, vars], use = "complete.obs")

ggcorrplot(corr_matrix, 
           method = "square",
           type = "lower",        
           lab = TRUE,
           lab_size = 3,
           colors = c("blue", "white", "red"), 
           show.legend = TRUE,    # Remove the legend
           ggtheme = theme_minimal()) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # Rotate labels
        panel.grid = element_blank()) +  # Remove grid lines
  scale_x_discrete(labels = var_labels) +  # Rename X-axis labels
  scale_y_discrete(labels = var_labels)    # Rename Y-axis labels


colnames(corr_matrix)  # x-axis order
rownames(corr_matrix)  # y-axis order