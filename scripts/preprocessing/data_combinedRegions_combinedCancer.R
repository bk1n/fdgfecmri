library(pacman)
p_load(tidyverse, infer)

fdg = readRDS("./data/fdgpet_data.rds")

fdg_lp = fdg %>%
  select(PATIENT_ID:CANC_ENDO, PT_SUV, LP_SUV, LP_ADC, LP_LA, LP_SA, HIST_LP) %>% 
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_LP)) %>%
  select(PT_SUV:HIST_LP) %>%
  rename(SUV = LP_SUV,
         ADC = LP_ADC,
         LA = LP_LA,
         SA = LP_SA,
         HIST = HIST_LP)

fdg_rp = fdg %>%
  select(PATIENT_ID:CANC_ENDO, PT_SUV, RP_SUV, RP_ADC, RP_LA, RP_SA, HIST_RP) %>%
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_RP)) %>%
  select(PT_SUV:HIST_RP) %>%
  rename(SUV = RP_SUV,
         ADC = RP_ADC,
         LA = RP_LA,
         SA = RP_SA,
         HIST = HIST_RP)

fdg_paln = fdg %>%
  select(PATIENT_ID:CANC_ENDO, PT_SUV, PALN_SUV, PALN_ADC, PALN_LA, PALN_SA, HIST_PALN) %>%
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_PALN)) %>%
  select(PT_SUV:HIST_PALN) %>%
  rename(SUV = PALN_SUV,
         ADC = PALN_ADC,
         LA = PALN_LA,
         SA = PALN_SA,
         HIST = HIST_PALN)

fdg_all = rbind(fdg_lp, fdg_rp, fdg_paln)

fdg_all = fdg_all %>%
  mutate(NTR = SUV / PT_SUV,
         STAR = SUV / ADC)

saveRDS(fdg_all, "./data/fdgpet_allPooled.rds")
