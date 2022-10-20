library(pacman)
p_load(tidyverse)

fdg_data = as_tibble(read.csv("./data/fdgpet_data.csv", header = T))

fdg_data_clean  = fdg_data %>%
  
  mutate(RP_SUV_CR = rowMeans(select(.,SUV_RP_1_CR, SUV_RP_2_CR), na.rm = T),
         LP_SUV_CR = rowMeans(select(.,SUV_LP_1_CR, SUV_LP_2_CR), na.rm = T),
         PALN_SUV_CR = rowMeans(select(.,SUV_PALN_1_CR, SUV_PALN_2_CR), na.rm = T)) %>%
  select(-c(SUV_LP_1_CR:SUV_PALN_2_CR)) %>% 
  
  mutate(RP_SUV = case_when(is.na(RP_SUV_CR) ~ RP_SUV,
                            TRUE ~ RP_SUV_CR),
         LP_SUV = case_when(is.na(LP_SUV_CR) ~ LP_SUV,
                            TRUE ~ LP_SUV_CR),
         PALN_SUV = case_when(is.na(PALN_SUV_CR) ~ PALN_SUV,
                            TRUE ~ PALN_SUV_CR)) %>%
  select(-c(RP_SUV_CR:PALN_SUV_CR)) %>%
  
  mutate(RP_ADC = rowMeans(select(., ADC_RP_1, ADC_RP_2), na.rm = T),
         LP_ADC = rowMeans(select(., ADC_LP_1, ADC_LP_2), na.rm = T),
         PALN_ADC = rowMeans(select(., ADC_PALN_1, ADC_PALN_2), na.rm = T)) %>%
  select(-c(ADC_LP_1:ADC_PALN_2)) %>%
  
  mutate(across(PT_SUV:PALN_SIZE, ~ ifelse(.x == "", NA, .x))) %>%
  mutate(across(RP_SIZE:PALN_SIZE, ~ str_trim(.x, side = "both"))) %>% 
  separate(RP_SIZE, into = c("RP_LA", "RP_SA"), sep = " x ") %>%
  separate(LP_SIZE, into = c("LP_LA", "LP_SA"), sep = " x ") %>%
  separate(PALN_SIZE, into = c("PALN_LA", "PALN_SA"), sep = " x ") %>%
  
  mutate(across(PT_SUV:PALN_SA, ~ as.numeric(.x))) %>%
  mutate(across(HIST_LP:HIST_PALN, ~ as.factor(.x))) %>%
  
  mutate(TRIAL_STATUS = ifelse(TRIAL_STATUS == "WITHDRAWN", NA, TRIAL_STATUS),
         TRIAL_STATUS = ifelse(TRIAL_STATUS == "", NA, TRIAL_STATUS)) %>%
  
  filter(!is.na(CANC_ENDO) & !is.na(CANC_CER)) %>%

  mutate(CANC = as.factor(if_else(CANC_ENDO == 1, "endo", "cer"))) %>%
  select(-c(CANC_ENDO, CANC_CER))


saveRDS(fdg_data_clean, "./data/fdgpet_data.rds")

ggplot(fdg_data_clean,
       aes(x = CANC,
           y = RP_SUV)) +
  geom_boxplot()
