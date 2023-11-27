library(tidyverse)

vis = read_csv('./data/MAPPING_expert_assessment_scores.csv')

quants = read_csv('./data/quant_allPooled.csv')

lp  = vis %>%
  select(`Patient ID`, contains('LP_')) %>%
  mutate(region = 'lp') %>%
  rename_with(.fn = function(col) gsub('LP_', '', col))

rp = vis %>%
  select(`Patient ID`, contains('RP_')) %>%
  mutate(region = 'rp') %>%
  rename_with(.fn = function(col) gsub('RP_', '', col))

paln = vis %>%
  select(`Patient ID`, contains('PALN_')) %>%
  mutate(region = 'paln') %>%
  rename_with(.fn = function(col) gsub('PALN_', '', col))

vis = rbind(lp, rp, paln) %>%
  rename(VIS_MORPH = 'MORPH',
         VIS_MRI = 'MRI',
         VIS_FDG = 'FDG',
         VIS_FEC = 'FEC',
         PATIENT_ID = 'Patient ID')

write.csv(vis, './data/visual_data.csv', row.names = F)
