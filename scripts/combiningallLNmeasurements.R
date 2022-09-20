library(pacman)
p_load(tidyverse, infer)

fdg = readRDS("./data/fdgpet_data.rds")

fdg_lp = fdg %>%
  select(PATIENT_ID:CANC_ENDO, LP_SUV, LP_LA, LP_SA, HIST_LP) %>%
  filter(!is.na(TRIAL_STATUS),
         !is.na(LP_SUV),
         !is.na(HIST_LP)) %>%
  select(LP_SUV:HIST_LP) %>%
  rename(SUV = LP_SUV,
         LA = LP_LA,
         SA = LP_SA,
         HIST = HIST_LP)

fdg_rp = fdg %>%
  select(PATIENT_ID:CANC_ENDO, RP_SUV, RP_LA, RP_SA, HIST_RP) %>%
  filter(!is.na(TRIAL_STATUS),
         !is.na(RP_SUV),
         !is.na(HIST_RP)) %>%
  select(RP_SUV:HIST_RP) %>%
  rename(SUV = RP_SUV,
         LA = RP_LA,
         SA = RP_SA,
         HIST = HIST_RP)

fdg_paln = fdg %>%
  select(PATIENT_ID:CANC_ENDO, PALN_SUV, PALN_LA, PALN_SA, HIST_PALN) %>%
  filter(!is.na(TRIAL_STATUS),
         !is.na(PALN_SUV),
         !is.na(HIST_PALN)) %>%
  select(PALN_SUV:HIST_PALN) %>%
  rename(SUV = PALN_SUV,
         LA = PALN_LA,
         SA = PALN_SA,
         HIST = HIST_PALN)

fdg_all = rbind(fdg_lp, fdg_rp, fdg_paln)

print(paste("Number of samples (pooled endo and cervical, all locations):", nrow(fdg_all)))

ggplot(fdg_all,
       aes(x = HIST,
           y = SUV)) +
  geom_boxplot() +
  geom_point()

t_test(fdg_all, SUV ~ HIST)
