library(pacman)
p_load(tidyverse, infer)

fdg = readRDS("./data/fdgpet_data.rds")

fdg_lp = fdg %>%
  select(TRIAL_STATUS, CANC, PT_SUV, LP_SUV, LP_ADC, LP_LA, LP_SA, HIST_LP) %>% 
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_LP)) %>%
  select(CANC:HIST_LP) %>%
  rename(SUV = LP_SUV,
         ADC = LP_ADC,
         LA = LP_LA,
         SA = LP_SA,
         HIST = HIST_LP) %>%
  mutate(region = "lp")

fdg_rp = fdg %>%
  select(TRIAL_STATUS, CANC, PT_SUV, RP_SUV, RP_ADC, RP_LA, RP_SA, HIST_RP) %>%
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_RP)) %>%
  select(CANC:HIST_RP) %>%
  rename(SUV = RP_SUV,
         ADC = RP_ADC,
         LA = RP_LA,
         SA = RP_SA,
         HIST = HIST_RP) %>%
  mutate(region = "rp")

fdg_paln = fdg %>%
  select(TRIAL_STATUS, CANC, PT_SUV, PALN_SUV, PALN_ADC, PALN_LA, PALN_SA, HIST_PALN) %>%
  filter(!is.na(TRIAL_STATUS),
         !is.na(HIST_PALN)) %>%
  select(CANC:HIST_PALN) %>%
  rename(SUV = PALN_SUV,
         ADC = PALN_ADC,
         LA = PALN_LA,
         SA = PALN_SA,
         HIST = HIST_PALN) %>%
  mutate(region = "paln")

fdg_all = rbind(fdg_lp, fdg_rp, fdg_paln)

fdg_all = fdg_all %>%
  mutate(NTR = SUV / PT_SUV,
         STAR = SUV / ADC,
         SNSA = SUV / SA)

saveRDS(fdg_all, "./data/fdgpet_allPooled.rds")
