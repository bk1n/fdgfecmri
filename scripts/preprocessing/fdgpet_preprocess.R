library(pacman)
p_load(tidyverse)

data = as_tibble(read.csv("./data/fdgpet_data.csv", header = T))

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


saveRDS(data_clean, "./data/processed_quant_data.rds")


