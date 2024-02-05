library(tidyverse)

data = as_tibble(read.csv("./data/fdgpet_data.csv", header = T))
vis = read_csv('./data/MAPPING_expert_assessment_scores.csv')
clin_data = as_tibble(read.csv("./data/clinical_data.csv", header = T))


# clean data ----
data_clean  = data %>%
  #get means of all CR for FDG
  mutate(FDG_SUV_RP_CR = rowMeans(select(.,FDG_SUV_RP_1_CR, FDG_SUV_RP_2_CR), na.rm = T),
         FDG_SUV_LP_CR = rowMeans(select(.,FDG_SUV_LP_1_CR, FDG_SUV_LP_2_CR), na.rm = T),
         FDG_SUV_PALN_CR = rowMeans(select(.,FDG_SUV_PALN_1_CR, FDG_SUV_PALN_2_CR), na.rm = T)) %>%
  select(-c(FDG_SUV_LP_1_CR:FDG_SUV_PALN_2_CR)) %>% 
  
  #get means of all CR for FEC
  mutate(FEC_SUV_RP_CR = rowMeans(select(.,FEC_SUV_RP_1_CR, FEC_SUV_RP_2_CR), na.rm = T),
         FEC_SUV_LP_CR = rowMeans(select(.,FEC_SUV_LP_1_CR, FEC_SUV_LP_2_CR), na.rm = T),
         FEC_SUV_PALN_CR = rowMeans(select(.,FEC_SUV_PALN_1_CR, FEC_SUV_PALN_2_CR), na.rm = T)) %>%
  select(-c(FEC_SUV_LP_1_CR:FEC_SUV_PALN_2_CR)) %>% 
  
  #get means of all CR for ADC
  mutate(ADC_PT_CR = rowMeans(select(., ADC_PT_1_CR, ADC_PT_2_CR), na.rm = T),
         ADC_RP_CR = rowMeans(select(., ADC_RP_1_CR, ADC_RP_2_CR), na.rm = T),
         ADC_LP_CR = rowMeans(select(., ADC_LP_1_CR, ADC_LP_2_CR), na.rm = T),
         ADC_PALN_CR = rowMeans(select(., ADC_PALN_1_CR, ADC_PALN_2_CR), na.rm = T)) %>%
  select(-c(ADC_PT_1_CR:ADC_PALN_2_CR)) %>%
  
  #FDG - if CR is NA, take Tara's read; else take CR
  mutate(FDG_SUV_RP = case_when(is.na(FDG_SUV_RP_CR) ~ FDG_SUV_RP,
                                TRUE ~ FDG_SUV_RP_CR),
         FDG_SUV_LP = case_when(is.na(FDG_SUV_LP_CR) ~ FDG_SUV_LP,
                                TRUE ~ FDG_SUV_LP_CR),
         FDG_SUV_PALN = case_when(is.na(FDG_SUV_PALN_CR) ~ FDG_SUV_PALN,
                                  TRUE ~ FDG_SUV_PALN_CR)) %>%
  select(-c(FDG_SUV_RP_CR:FDG_SUV_PALN_CR)) %>%
  
  #FEC - if CR is NA, take Tara's read; else take CR
  mutate(FEC_SUV_RP = case_when(is.na(FEC_SUV_RP_CR) ~ FEC_SUV_RP,
                                TRUE ~ FEC_SUV_RP_CR),
         FEC_SUV_LP = case_when(is.na(FEC_SUV_LP_CR) ~ FEC_SUV_LP,
                                TRUE ~ FEC_SUV_LP_CR),
         FEC_SUV_PALN = case_when(is.na(FEC_SUV_PALN_CR) ~ FEC_SUV_PALN,
                                  TRUE ~ FEC_SUV_PALN_CR)) %>%
  select(-c(FEC_SUV_RP_CR:FEC_SUV_PALN_CR)) %>%
  
  rename(ADC_PT = ADC_PT_CR, ADC_RP = ADC_RP_CR, ADC_LP = ADC_LP_CR, ADC_PALN = ADC_PALN_CR) %>%
  relocate(ADC_PT:ADC_PALN, .after = FEC_SUV_PALN) %>%
  
  #remove ""
  mutate(across(FDG_SUV_PT:FDG_PALN_SIZE, ~ ifelse(.x == "", NA, .x))) %>%
  
  #split LN dimensions for FEC and FDG
  mutate(across(FEC_RP_SIZE:FDG_PALN_SIZE, ~ str_trim(.x, side = "both"))) %>% 
  separate(FEC_RP_SIZE, into = c("FEC_RP_LA", "FEC_RP_SA"), sep = " x ") %>%
  separate(FEC_LP_SIZE, into = c("FEC_LP_LA", "FEC_LP_SA"), sep = " x ") %>%
  separate(FEC_PALN_SIZE, into = c("FEC_PALN_LA", "FEC_PALN_SA"), sep = " x ") %>%
  separate(FDG_RP_SIZE, into = c("FDG_RP_LA", "FDG_RP_SA"), sep = " x ") %>%
  separate(FDG_LP_SIZE, into = c("FDG_LP_LA", "FDG_LP_SA"), sep = " x ") %>%
  separate(FDG_PALN_SIZE, into = c("FDG_PALN_LA", "FDG_PALN_SA"), sep = " x ") %>%
  
  mutate(across(FDG_SUV_PT:FDG_PALN_SA, ~ as.numeric(.x))) %>%
  mutate(across(HIST_LP:HIST_PALN, ~ as.factor(.x))) %>%
  
  mutate(TRIAL_STATUS = ifelse(TRIAL_STATUS == "WITHDRAWN", NA, TRIAL_STATUS),
         TRIAL_STATUS = ifelse(TRIAL_STATUS == "", NA, TRIAL_STATUS)) %>%
  
  filter(!is.na(CANC_ENDO) & !is.na(CANC_CER)) %>%
  
  mutate(CANC = as.factor(if_else(CANC_ENDO == 1, "endo", "cer"))) %>%
  select(-c(CANC_ENDO, CANC_CER))

# regional data ----
# convert data from wide format into longer (with a col for regions)
#lp
lp = data_clean %>%
  select(PATIENT_ID, TRIAL_STATUS, CANC, FDG_SUV_PT, FEC_SUV_PT, ADC_PT, FDG_SUV_LP, FEC_SUV_LP, ADC_LP, FEC_LP_LA, FEC_LP_SA, FDG_LP_LA, FDG_LP_SA, HIST_LP) %>% 
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_LP)) %>%
  select(PATIENT_ID, CANC:HIST_LP) %>%
  rename(FDG_SUV = FDG_SUV_LP,
         FEC_SUV = FEC_SUV_LP,
         ADC = ADC_LP,
         FDG_LA = FDG_LP_LA,
         FDG_SA = FDG_LP_SA,
         FEC_LA = FEC_LP_LA,
         FEC_SA = FEC_LP_SA,
         HIST = HIST_LP) %>%
  mutate(region = "lp")

#rp
rp = data_clean %>%
  select(PATIENT_ID, TRIAL_STATUS, CANC, FDG_SUV_PT, FEC_SUV_PT, ADC_PT, FDG_SUV_RP, FEC_SUV_RP, ADC_RP, FEC_RP_LA, FEC_RP_SA, FDG_RP_LA, FDG_RP_SA, HIST_RP) %>% 
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_RP)) %>%
  select(PATIENT_ID, CANC:HIST_RP) %>%
  rename(FDG_SUV = FDG_SUV_RP,
         FEC_SUV = FEC_SUV_RP,
         ADC = ADC_RP,
         FDG_LA = FDG_RP_LA,
         FDG_SA = FDG_RP_SA,
         FEC_LA = FEC_RP_LA,
         FEC_SA = FEC_RP_SA,
         HIST = HIST_RP) %>%
  mutate(region = "rp")

#paln
paln = data_clean %>%
  select(PATIENT_ID, TRIAL_STATUS, CANC, FDG_SUV_PT, FEC_SUV_PT, ADC_PT, FDG_SUV_PALN, FEC_SUV_PALN, ADC_PALN, FEC_PALN_LA, FEC_PALN_SA, FDG_PALN_LA, FDG_PALN_SA, HIST_PALN) %>% 
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_PALN)) %>%
  select(PATIENT_ID, CANC:HIST_PALN) %>%
  rename(FDG_SUV = FDG_SUV_PALN,
         FEC_SUV = FEC_SUV_PALN,
         ADC = ADC_PALN,
         FDG_LA = FDG_PALN_LA,
         FDG_SA = FDG_PALN_SA,
         FEC_LA = FEC_PALN_LA,
         FEC_SA = FEC_PALN_SA,
         HIST = HIST_PALN) %>%
  mutate(region = "paln")

# combine
all = rbind(lp, rp, paln)

# calc quant measures
all = all %>%
  mutate(FDG_NTR = FDG_SUV / FDG_SUV_PT,
         FEC_NTR = FEC_SUV / FEC_SUV_PT,
         FDG_STAR = FDG_SUV / ADC,
         FEC_STAR = FEC_SUV / ADC,
         FDG_SNSA = FDG_SUV / FDG_SA,
         FEC_SNSA = FEC_SUV / FEC_SA,
         ADC_NTR = ADC / ADC_PT)

# combine w visual measures ----
# convert vis to long format
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

# combine
all = all %>%
  left_join(vis, by = c('PATIENT_ID', 'region'))

# combine w clinical data ----
# clean clinical data
clin_data = clin_data %>%
  #collapse multiple ENDO histology columns into single 
  mutate(ENDO_HIST = case_when(ENDO_HIST_ENDO == T ~ 'endo',
                               ENDO_HIST_SERCC == T ~ 'sercc',
                               ENDO_HIST_MM == T ~ 'mm',
                               ENDO_HIST_OTHER == T ~ 'other',
                               ENDO_HIST_MIXED == T ~ 'mixed')) %>%
  select(-c(ENDO_HIST_ENDO:ENDO_HIST_MIXED)) %>%
  
  mutate(across(CER_HIST:ENDO_HIST, ~ ifelse(.x == "", NA, .x))) %>%
  
  mutate(CER_LVSI = ifelse(CER_LVSI == 1, 'Yes', CER_LVSI)) %>%
  mutate(BMI = FDG_PW / (FDG_PH / 100)**2) %>%
  select(-TRIAL_STATUS)

# join clin data w all
all = all %>%
  left_join(clin_data, by = 'PATIENT_ID')

# save to csv ----
write.csv(all, './data/processed_data.csv', row.names = F)
