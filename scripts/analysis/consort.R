library(visR)
library(tidyverse)

mpg = read_csv('./data/MPG_DATA_COMPILED.csv')

quants = read_csv('./data/quant_data.csv')

pool = read_csv('./data/quant_allPooled.csv')

data = readRDS("./data/processed_quant_data.rds")

data %>%
  select(contains('HIST')) %>%
  filter_all(., any_vars(!is.na(.)))

data %>%
  # filter(CANC == 'cer') %>%
  filter_at(vars(HIST_LP, HIST_RP, HIST_PALN), any_vars(!is.na(.))) %>%
  select(-c(PATIENT_ID, TRIAL_STATUS, CANC),
         # -contains('HIST'),
         -contains('PT'),
         -contains('SA'),
         -contains('LA')) %>% 
  filter_all(., any_vars(!is.na(.)))

# endo
pool %>%
  filter(CANC == 'endo') %>%
  sapply(function(x) sum(!is.na(x))) 

pool %>%
  filter(CANC == 'endo') %>%
  filter(region == 'lp' | region == 'rp') %>% 
  sapply(function(x) sum(!is.na(x))) 
  

pool %>%
  filter(CANC == 'endo') %>%
  filter(region == 'paln') %>% 
  sapply(function(x) sum(!is.na(x))) 


pool %>%
  filter(CANC == 'endo') %>%
  filter(!is.na(FDG_SUV)) 

# cer
pool %>%
  filter(CANC == 'cer') %>%
  sapply(function(x) sum(!is.na(x))) 

pool %>%
  filter(CANC == 'cer') %>%
  filter(region == 'lp' | region == 'rp') %>% 
  sapply(function(x) sum(!is.na(x))) 


pool %>%
  filter(CANC == 'cer') %>%
  filter(region == 'paln') %>% 
  sapply(function(x) sum(!is.na(x))) 