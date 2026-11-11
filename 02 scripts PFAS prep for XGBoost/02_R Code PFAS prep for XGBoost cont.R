###################################################################
# R Code: Prepping PFAS data for XGBoost
# Goal: Get a matrix with labels for each CWS in NJ over time
# Structure: We will have CWS, Season-year, and label
# Labels: 0 = non-detect, 1=detect, < MCL, 2 = detect, > MCL
# Created by: Hari Iyer
# Date: 11/21/2025

# Note: sourced from
# C:/Users/hi97/Box/Grants/K award/NJDEP/Time series/Final Scripts/R Code_pfas_2018_2025_v6.R
###################################################################


library(tidyverse)
library(stringr)
library(sf)
library(tmap)


## install.packages("shiny")
## install.packages(c("shinyjs", "kableExtra"))
library(shinyjs)
library(kableExtra)
library(shiny)
library(cols4all)
## cols4all::c4a_gui()  # Opens an interactive color picker

library(tmap)
## install.packages("tidyquant")
library(tidyquant)
library(nlme)
library(lmtest)
library(lubridate)

library(broom)
## install.packages("broom.mixed")
library(broom.mixed)
library(tseries)

library(forecast)

####################################################
# Get all data cleaned from multiple sources
####################################################

pfas_update <- read.csv("C:/Users/hi97/Box/Grants/K award/NJDEP/pfas_nj_pwsid_2018_2025.csv", stringsAsFactors = F, header = T)
names(pfas_update)

pfas_update_clean <- pfas_update %>%
  rename(SYSTEM.NAME = System.Name,
         SAMPLE.COLLECTION.DATE = Sample.Date,
         PARAMETER.NAME = Analyte,
         LABORATORY.NAME = Lab.Name..NJ.Lab.ID..) %>%
  dplyr::select(PWSID, SYSTEM.NAME, SAMPLE.COLLECTION.DATE,
                PARAMETER.NAME, LABORATORY.NAME, Analytical.Method, Result) %>%
  filter(PARAMETER.NAME %in% c("PERFLUOROCTANOIC ACID (PFOA)","PERFLUOROCTANE SULFONIC ACID (PFOS)","PERFLUORONONANOIC ACID (PFNA)")) 

table(pfas_update_clean$Analytical.Method)


pfas_update_clean$QUALIFIER <- ifelse(substr(pfas_update_clean$Result,1,1) == "<","<",NA)
pfas_update_clean$UNIT <- sapply(str_split(pfas_update_clean$Result, " "), `[`, 2)
pfas_update_clean$RESULT <- str_extract_all(pfas_update_clean$Result, "\\d+\\.?\\d*")
pfas_update_clean$RESULT <- as.numeric(unlist(pfas_update_clean$RESULT))



## Only PFOA, PFOS, PFNA




table(pfas_update_clean$QUALIFIER)

table(pfas_update_clean$QUALIFIER, pfas_update_clean$PARAMETER.NAME)
table(pfas_update_clean$PARAMETER.NAME)


## Update result so it's only ng/L

pfas_update_clean$RESULT_NG <- ifelse(pfas_update_clean$UNIT == "NG/L",pfas_update_clean$RESULT,
                                      ifelse(pfas_update_clean$UNIT == "UG/L", pfas_update_clean$RESULT * 1000,NA))

pfas_update_clean$RESULT_NG[pfas_update_clean$RESULT_NG > 100]<-2/sqrt(2)


pfas_update_clean$RESULT_NG <- ifelse(
  !is.na(pfas_update_clean$QUALIFIER) & pfas_update_clean$QUALIFIER == "<",
  pfas_update_clean$RESULT_NG / sqrt(2),
  pfas_update_clean$RESULT_NG
)
table(pfas_update_clean$UNIT) 

## restrict to PFOA, PFOS, PFNA

pfas_update_clean <- pfas_update_clean %>%
  filter(PARAMETER.NAME %in% c("PERFLUOROCTANOIC ACID (PFOA)","PERFLUOROCTANE SULFONIC ACID (PFOS)","PERFLUORONONANOIC ACID (PFNA)")) 


length(unique(pfas_update_clean$PWSID)) 

table(pfas_update_clean$PARAMETER.NAME)
table(pfas_update_clean$QUALIFIER)

## get # of PWSID
length(unique(pfas_update_clean$PWSID)) 


pfas_og <- read.csv("C:/Users/hi97/Box/Grants/K award/NJDEP/2018.07.26.Query for OPRA 233491_PFAS.csv",
                    stringsAsFactors = F, header = T) 



pfas_og_sm <- pfas_og %>%
  dplyr::select(PWSID, SYSTEM.NAME, SAMPLE.COLLECTION.DATE, SAMPLE.TYPE, PARAMETER.NAME, ANALYTICAL.METHOD, LABORATORY.NAME, QUALIFIER, UNIT,
                RESULT) %>%
  filter(PARAMETER.NAME %in% c("Perfluoro octanoic acid - PFOA","Perfluoro octanesulfonic acid - PFOS","Perfluoro-n-nonanoic acid - PFNA")) %>%
  filter(! SAMPLE.TYPE %in% c("Raw", "Raw (Combined wells)")) %>%
  filter(!QUALIFIER %in% c("BE", "BR", "FA", "J", "J(FA)", "JB", "R", "S6")) ## 3091



pfas_og_sm$RESULT_NG <- pfas_og_sm$RESULT
## apply LOD
pfas_og_sm$RESULT_NG <- ifelse(
  !is.na(pfas_og_sm$QUALIFIER) & pfas_og_sm$QUALIFIER == "<",
  pfas_og_sm$RESULT_NG / sqrt(2),
  pfas_og_sm$RESULT_NG
) 


table(pfas_og_sm$QUALIFIER, pfas_og_sm$PARAMETER.NAME)
table(pfas_og_sm$PARAMETER.NAME)

pfas_og_sm <- pfas_og_sm %>%
  filter(PARAMETER.NAME %in% c("Perfluoro octanoic acid - PFOA","Perfluoro octanesulfonic acid - PFOS","Perfluoro-n-nonanoic acid - PFNA"))


table(pfas_og_sm$PWSID)



length(unique(pfas_og_sm$PWSID)) 


pfas_update_clean <- pfas_update_clean %>%
  dplyr::select(-c(Result,Analytical.Method))

pfas_og_sm <- pfas_og_sm %>%
  dplyr::select(-c(ANALYTICAL.METHOD, SAMPLE.TYPE)) %>%
  mutate(QUALIFIER = if_else(QUALIFIER == "<", "<", NA_character_))


pfas_full <- rbind(pfas_og_sm, pfas_update_clean) 
ncol(pfas_og_sm)
names(pfas_og_sm)
ncol(pfas_update_clean)
names(pfas_update_clean)
length(unique(pfas_full$PWSID)) 

#############################################################
# Create matrix
# Variables
# PWS/CWS
# Season
# Year
# Category (0=ND, 1=<MCL, 2=>MCL)
#############################################################

pfas_xg <- pfas_full %>%
  mutate(pfas_label = case_when(
    PARAMETER.NAME %in% c("Perfluoro octanoic acid - PFOA","PERFLUOROCTANOIC ACID (PFOA)") & QUALIFIER == "<"  ~ "ND",
    PARAMETER.NAME %in% c("Perfluoro octanoic acid - PFOA","PERFLUOROCTANOIC ACID (PFOA)") & is.na(QUALIFIER) & RESULT_NG < 14  ~ "<MCL",
    PARAMETER.NAME %in% c("Perfluoro octanoic acid - PFOA","PERFLUOROCTANOIC ACID (PFOA)") & is.na(QUALIFIER) & RESULT_NG >= 14  ~ ">=MCL",
    PARAMETER.NAME %in% c("Perfluoro octanesulfonic acid - PFOS","PERFLUOROCTANE SULFONIC ACID (PFOS)") & QUALIFIER == "<"  ~ "ND",
    PARAMETER.NAME %in% c("Perfluoro octanesulfonic acid - PFOS","PERFLUOROCTANE SULFONIC ACID (PFOS)") & is.na(QUALIFIER) & RESULT_NG < 13  ~ "<MCL",
    PARAMETER.NAME %in% c("Perfluoro octanesulfonic acid - PFOS","PERFLUOROCTANE SULFONIC ACID (PFOS)") & is.na(QUALIFIER) & RESULT_NG >= 13  ~ ">=MCL",
    PARAMETER.NAME %in% c("Perfluoro-n-nonanoic acid - PFNA","PERFLUORONONANOIC ACID (PFNA)") & QUALIFIER == "<"  ~ "ND",
    PARAMETER.NAME %in% c("Perfluoro-n-nonanoic acid - PFNA","PERFLUORONONANOIC ACID (PFNA)") & is.na(QUALIFIER) & RESULT_NG < 13  ~ "<MCL",
    PARAMETER.NAME %in% c("Perfluoro-n-nonanoic acid - PFNA","PERFLUORONONANOIC ACID (PFNA)") & is.na(QUALIFIER) & RESULT_NG >= 13  ~ ">=MCL",
    is.na(PARAMETER.NAME) ~ "Missing",
    TRUE ~ "Missing"
  ),
  
  # numeric label for XGBoost: 0 = ND, 1 = <MCL, 2 = >=MCL
  pfas_label_num = case_when(
    pfas_label == "ND"    ~ 0L,
    pfas_label == "<MCL"  ~ 1L,
    pfas_label == ">=MCL" ~ 2L,
    TRUE                  ~ NA_integer_
  ),
  season = case_when(
    month(as.Date(SAMPLE.COLLECTION.DATE, format="%m/%d/%Y")) %in% 1:3 ~ 0L,
    month(as.Date(SAMPLE.COLLECTION.DATE, format="%m/%d/%Y")) %in% 4:6 ~ 1L,
    month(as.Date(SAMPLE.COLLECTION.DATE, format="%m/%d/%Y")) %in% 7:9 ~ 2L,
    month(as.Date(SAMPLE.COLLECTION.DATE, format="%m/%d/%Y")) %in% 10:12 ~ 3L,
  TRUE                                              ~ NA_integer_
  ),
  year = year(as.Date(SAMPLE.COLLECTION.DATE, format="%m/%d/%Y")),
  pfas_name = case_when(
    PARAMETER.NAME == "Perfluoro octanoic acid - PFOA" ~ "PFOA",
    PARAMETER.NAME == "PERFLUOROCTANOIC ACID (PFOA)" ~ "PFOA",
    PARAMETER.NAME == "Perfluoro octanesulfonic acid - PFOS" ~ "PFOS",
    PARAMETER.NAME == "PERFLUOROCTANE SULFONIC ACID (PFOS)" ~ "PFOS",
    PARAMETER.NAME == "Perfluoro-n-nonanoic acid - PFNA" ~ "PFNA",
    PARAMETER.NAME == "PERFLUORONONANOIC ACID (PFNA)" ~ "PFNA",
    is.na(PARAMETER.NAME) ~ "Missing",
    TRUE ~ "Missing"
  )
) %>%
dplyr::select(PWSID, pfas_name, pfas_label, pfas_label_num, season, year, RESULT_NG)



###############################################################################
# UCMR3
###############################################################################

######################
## read in UCMR3 ##
######################

ucmr3.state <- read.delim("C:/Users/hi97/Box/PFAS Grant/UCMR/ucmr3/ucmr3-occurrence-data-by-state/UCMR3_All_MA_WY.txt", sep="\t", dec = ".", header=T, fileEncoding="latin1")
ucmr3.nj <- ucmr3.state %>%
  filter(State == "NJ" & Contaminant %in% c("PFOA","PFOS","PFNA")) %>%
  mutate(RESULT_NG = as.numeric(AnalyticalResultValue)*1000) %>%
  mutate(pfas_label = case_when(
    Contaminant == "PFOA" & AnalyticalResultsSign == "<"  ~ "ND",
    Contaminant == "PFOA" & AnalyticalResultsSign == "=" & RESULT_NG < 20  ~ "<20ng/L",
    Contaminant == "PFOA" & AnalyticalResultsSign == "=" & RESULT_NG >= 20  ~ ">=20ng/L",
    Contaminant == "PFOS" & AnalyticalResultsSign == "<"  ~ "ND",
    Contaminant == "PFOS" & AnalyticalResultsSign == "=" & RESULT_NG < 40  ~ "<40ng/L",
    Contaminant == "PFOS" & AnalyticalResultsSign == "=" & RESULT_NG >= 40  ~ ">=40ng/L",
    Contaminant == "PFNA" & AnalyticalResultsSign == "<"  ~ "ND",
    Contaminant == "PFNA" & AnalyticalResultsSign == "=" & RESULT_NG < 20  ~ "<40ng/L",
    Contaminant == "PFNA" & AnalyticalResultsSign == "=" & RESULT_NG >= 20  ~ ">=40ng/L",
    is.na(Contaminant) ~ "Missing",
    TRUE ~ "Missing"
  ),
  
  # numeric label for XGBoost: 0 = ND, 1 = <MCL, 2 = >=MCL
  pfas_label_num = case_when(
    pfas_label == "ND"    ~ 0L,
    pfas_label %in% c("<20ng/L","<40ng/L")  ~ 1L,
    pfas_label %in% c(">=20ng/L",">=40ng/L") ~ 2L,
    TRUE                  ~ NA_integer_
  ),
  season = case_when(
    month(as.Date(CollectionDate, format="%m/%d/%Y")) %in% 1:3 ~ 0L,
    month(as.Date(CollectionDate, format="%m/%d/%Y")) %in% 4:6 ~ 1L,
    month(as.Date(CollectionDate, format="%m/%d/%Y")) %in% 7:9 ~ 2L,
    month(as.Date(CollectionDate, format="%m/%d/%Y")) %in% 10:12 ~ 3L,
    TRUE                                              ~ NA_integer_
  ),
  year = year(as.Date(CollectionDate, format="%m/%d/%Y")),
  pfas_name = case_when(
    Contaminant == "PFOA" ~ "PFOA",
    Contaminant == "PFOS" ~ "PFOS",
    Contaminant == "PFNA" ~ "PFNA",
    is.na(Contaminant) ~ "Missing",
    TRUE ~ "Missing"
  )
  ) %>%
  dplyr::select(PWSID, pfas_name, pfas_label, pfas_label_num, season, year, RESULT_NG)
rm(ucmr3.state)
length(unique(ucmr3.nj$PWSID)) 

###############################################################################
# UCMR5
###############################################################################


ucmr5.state <- read.delim("C:/Users/hi97/Box/PFAS Grant/UCMR/ucmr5/ucmr5-occurrence-data-by-state/UCMR5_All_MA_WY.txt", sep="\t", dec = ".", header=T, fileEncoding="latin1")

table(ucmr5.state$Contaminant)

ucmr5.nj <- ucmr5.state %>%
  filter(State == "NJ" & Contaminant %in% c("PFOA","PFOS","PFNA")) %>%
  mutate(RESULT_NG = as.numeric(AnalyticalResultValue)*1000) %>%
  mutate(pfas_label = case_when(
    Contaminant == "PFOA" & AnalyticalResultsSign == "<"  ~ "ND",
    Contaminant == "PFOA" & AnalyticalResultsSign == "=" & RESULT_NG < 14  ~ "<MCL",
    Contaminant == "PFOA" & AnalyticalResultsSign == "=" & RESULT_NG >= 14  ~ ">=MCL",
    Contaminant == "PFOS" & AnalyticalResultsSign == "<"  ~ "ND",
    Contaminant == "PFOS" & AnalyticalResultsSign == "=" & RESULT_NG < 13  ~ "<MCL",
    Contaminant == "PFOS" & AnalyticalResultsSign == "=" & RESULT_NG >= 13  ~ ">=MCL",
    Contaminant == "PFNA" & AnalyticalResultsSign == "<"  ~ "ND",
    Contaminant == "PFNA" & AnalyticalResultsSign == "=" & RESULT_NG < 13  ~ "<MCL",
    Contaminant == "PFNA" & AnalyticalResultsSign == "=" & RESULT_NG >= 13  ~ ">=MCL",
    is.na(Contaminant) ~ "Missing",
    TRUE ~ "Missing"
  ),
  
  # numeric label for XGBoost: 0 = ND, 1 = <MCL, 2 = >=MCL
  pfas_label_num = case_when(
    pfas_label == "ND"    ~ 0L,
    pfas_label == "<MCL"  ~ 1L,
    pfas_label == ">=MCL" ~ 2L,
    TRUE                  ~ NA_integer_
  ),
  season = case_when(
    month(as.Date(CollectionDate, format="%m/%d/%Y")) %in% 1:3 ~ 0L,
    month(as.Date(CollectionDate, format="%m/%d/%Y")) %in% 4:6 ~ 1L,
    month(as.Date(CollectionDate, format="%m/%d/%Y")) %in% 7:9 ~ 2L,
    month(as.Date(CollectionDate, format="%m/%d/%Y")) %in% 10:12 ~ 3L,
    TRUE                                              ~ NA_integer_
  ),
  year = year(as.Date(CollectionDate, format="%m/%d/%Y")),
  pfas_name = case_when(
    Contaminant == "PFOA" ~ "PFOA",
    Contaminant == "PFOS" ~ "PFOS",
    Contaminant == "PFNA" ~ "PFNA",
    is.na(Contaminant) ~ "Missing",
    TRUE ~ "Missing"
  )
  ) %>%
  dplyr::select(PWSID, pfas_name, pfas_label, pfas_label_num, season, year, RESULT_NG)
rm(ucmr5.state)
length(unique(ucmr5.nj$PWSID))

table(ucmr5.nj$pfas_label)

############################################################################
# Combine all data
############################################################################

pfas_full <- rbind(pfas_xg, ucmr3.nj, ucmr5.nj) 

summary(pfas_full$RESULT_NG) 

pfas_full_miss <- pfas_full %>%
  filter(is.na(RESULT_NG))

table()

save(pfas_full, file = "C:/Users/hi97/Box/PFAS Grant/XGBoost/data/pfas_xg.Rda")



dups <- pfas_full %>%
  filter(duplicated(.) | duplicated(., fromLast = TRUE)) %>% 
  arrange(across(everything()))

dups_counts <- dups %>%
  count(across(everything())) %>%
  filter(n > 1) 

table(dups_counts$n)

pfas_unique <- pfas_full %>% 
  distinct() 

length(unique(pfas_full$PWSID)) 
