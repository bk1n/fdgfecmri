library(visR)
library(tidyverse)

data = readRDS("./data/processed_quant_data.rds")

# nodal status confirmed by histology
data %>%
  select(contains('HIST')) %>%
  filter_all(., any_vars(!is.na(.)))

## lp hist
data %>%
   filter(!is.na(HIST_LP))
## rp hist
data %>%
  filter(!is.na(HIST_RP))
## paln hist
data %>%
  filter(!is.na(HIST_PALN))

# had quant imaging data
cols = data %>%
  select(contains('ADC') | contains('FDG') | contains('FEC'),
         -contains('PT'),
         -contains('SA'),
         -contains('LA')) %>%
  colnames
## any canc
data %>%
  select(contains('ADC') | contains('FDG') | contains('FEC'),
         -contains('PT'),
         -contains('SA'),
         -contains('LA')) %>% 
  filter_all(., any_vars(!is.na(.)))
## endo
data %>%
  filter_at(cols, any_vars(!is.na(.))) %>%
  filter(CANC == 'endo')
## cer
data %>%
  filter_at(cols, any_vars(!is.na(.))) %>%
  filter(CANC == 'cer')

# had fdg
cols = data %>%
  select(contains('FDG'),
         -contains('PT'),
         -contains('SA'),
         -contains('LA')) %>%
  colnames
## all canc
data %>%
  filter_at(cols, any_vars(!is.na(.)))
## endo 
data %>%
  filter_at(cols, any_vars(!is.na(.))) %>%
  filter(CANC == 'endo')
## cer
data %>%
  filter_at(cols, any_vars(!is.na(.))) %>%
  filter(CANC == 'cer')

# had fec
cols = data %>%
  select(contains('FEC'),
         -contains('PT'),
         -contains('SA'),
         -contains('LA')) %>%
  colnames
## all canc
data %>%
  filter_at(cols, any_vars(!is.na(.)))
## endo 
data %>%
  filter_at(cols, any_vars(!is.na(.))) %>%
  filter(CANC == 'endo')
## cer
data %>%
  filter_at(cols, any_vars(!is.na(.))) %>%
  filter(CANC == 'cer')

# had adc
cols = data %>%
  select(contains('ADC'),
         -contains('PT'),
         -contains('SA'),
         -contains('LA')) %>%
  colnames
## all canc
data %>%
  filter_at(cols, any_vars(!is.na(.)))
## endo 
data %>%
  filter_at(cols, any_vars(!is.na(.))) %>%
  filter(CANC == 'endo')
## cer
data %>%
  filter_at(cols, any_vars(!is.na(.))) %>%
  filter(CANC == 'cer')



## NODES ####
pool = read_csv('./data/quant_allPooled.csv')

# fdg
## endo
### all
pool %>%
  filter(CANC == 'endo') %>%
  select('FDG_SUV',
         'region') %>%
  filter(!is.na(FDG_SUV)) 
### pelvic
pool %>%
  filter(CANC == 'endo') %>%
  select('FDG_SUV',
         'region') %>%
  filter(!is.na(FDG_SUV)) %>%
  filter(region == 'lp' | region == 'rp')
### pa
pool %>%
  filter(CANC == 'endo') %>%
  select('FDG_SUV',
         'region') %>%
  filter(!is.na(FDG_SUV)) %>%
  filter(region == 'paln')

# fdg
## cer
### all
pool %>%
  filter(CANC == 'cer') %>%
  select('FDG_SUV',
         'region') %>%
  filter(!is.na(FDG_SUV)) 
### pelvic
pool %>%
  filter(CANC == 'cer') %>%
  select('FDG_SUV',
         'region') %>%
  filter(!is.na(FDG_SUV)) %>%
  filter(region == 'lp' | region == 'rp')
### pa
pool %>%
  filter(CANC == 'cer') %>%
  select('FDG_SUV',
         'region') %>%
  filter(!is.na(FDG_SUV)) %>%
  filter(region == 'paln')

# fec
## endo
### all
pool %>%
  filter(CANC == 'endo') %>%
  select('FEC_SUV',
         'region') %>%
  filter(!is.na(FEC_SUV)) 
### pelvic
pool %>%
  filter(CANC == 'endo') %>%
  select('FEC_SUV',
         'region') %>%
  filter(!is.na(FEC_SUV)) %>%
  filter(region == 'lp' | region == 'rp')
### pa
pool %>%
  filter(CANC == 'endo') %>%
  select('FEC_SUV',
         'region') %>%
  filter(!is.na(FEC_SUV)) %>%
  filter(region == 'paln')

# fec
## cer
### all
pool %>%
  filter(CANC == 'cer') %>%
  select('FEC_SUV',
         'region') %>%
  filter(!is.na(FEC_SUV)) 
### pelvic
pool %>%
  filter(CANC == 'cer') %>%
  select('FEC_SUV',
         'region') %>%
  filter(!is.na(FEC_SUV)) %>%
  filter(region == 'lp' | region == 'rp')
### pa
pool %>%
  filter(CANC == 'cer') %>%
  select('FEC_SUV',
         'region') %>%
  filter(!is.na(FEC_SUV)) %>%
  filter(region == 'paln')

# mri
## endo
### all
pool %>%
  filter(CANC == 'endo') %>%
  select('ADC',
         'region') %>%
  filter(!is.na(ADC)) 
### pelvic
pool %>%
  filter(CANC == 'endo') %>%
  select('ADC',
         'region') %>%
  filter(!is.na(ADC)) %>%
  filter(region == 'lp' | region == 'rp')
### pa
pool %>%
  filter(CANC == 'endo') %>%
  select('ADC',
         'region') %>%
  filter(!is.na(ADC)) %>%
  filter(region == 'paln')

# mri
## cer
### all
pool %>%
  filter(CANC == 'cer') %>%
  select('ADC',
         'region') %>%
  filter(!is.na(ADC)) 
### pelvic
pool %>%
  filter(CANC == 'cer') %>%
  select('ADC',
         'region') %>%
  filter(!is.na(ADC)) %>%
  filter(region == 'lp' | region == 'rp')
### pa
pool %>%
  filter(CANC == 'cer') %>%
  select('ADC',
         'region') %>%
  filter(!is.na(ADC)) %>%
  filter(region == 'paln')
