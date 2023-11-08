library(pacman)
p_load(tidyverse, infer)

data = readRDS("./data/processed_quant_data.rds")

#chunk data into regions
lp = data %>%
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

rp = data %>%
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

paln = data %>%
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

all = rbind(lp, rp, paln)

all = all %>%
  mutate(FDG_NTR = FDG_SUV / FDG_SUV_PT,
         FEC_NTR = FEC_SUV / FEC_SUV_PT,
         FDG_STAR = FDG_SUV / ADC,
         FEC_STAR = FEC_SUV / ADC,
         FDG_SNSA = FDG_SUV / FDG_SA,
         FEC_SNSA = FEC_SUV / FEC_SA,
         ADC_NTR = ADC / ADC_PT)

saveRDS(all, "./data/quant_allPooled.rds")
write.csv(all, './data/quant_allPooled.csv', row.names = F)
