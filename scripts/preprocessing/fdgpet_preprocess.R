library(pacman)
p_load(tidyverse)

fdg_data = as_tibble(read.csv("./data/fdgpet_data.csv", header = T))

fdg_data_clean = fdg_data %>%
  mutate(across(PT_SUV:PALN_SIZE, ~ ifelse(.x == "", NA, .x))) %>%
  mutate(across(RP_SIZE:PALN_SIZE, ~ str_trim(.x, side = "both"))) %>% 
  separate(RP_SIZE, into = c("RP_LA", "RP_SA"), sep = " x ") %>%
  separate(LP_SIZE, into = c("LP_LA", "LP_SA"), sep = " x ") %>%
  separate(PALN_SIZE, into = c("PALN_LA", "PALN_SA"), sep = " x ") %>%
  mutate(across(PT_SUV:PALN_SA, ~ as.numeric(.x))) %>%
  mutate(TRIAL_STATUS = ifelse(TRIAL_STATUS == "WITHDRAWN", NA, TRIAL_STATUS),
         TRIAL_STATUS = ifelse(TRIAL_STATUS == "", NA, TRIAL_STATUS)) %>%
  mutate(across(HIST_LP:HIST_PALN, ~ as.factor(.x)))

saveRDS(fdg_data_clean, "./data/fdgpet_data.rds")
