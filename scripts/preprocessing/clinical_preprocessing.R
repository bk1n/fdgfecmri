library(pacman)
p_load(tidyverse)

data = as_tibble(read.csv("./data/clinical_data.csv", header = T))

data_clean = data %>%
  #collapse multiple ENDO histology columns into single 
  mutate(ENDO_HIST = case_when(ENDO_HIST_ENDO == T ~ 'endo',
                               ENDO_HIST_SERCC == T ~ 'sercc',
                               ENDO_HIST_MM == T ~ 'mm',
                               ENDO_HIST_OTHER == T ~ 'other',
                               ENDO_HIST_MIXED == T ~ 'mixed')) %>%
  select(-c(ENDO_HIST_ENDO:ENDO_HIST_MIXED)) %>%
  
  mutate(across(CER_HIST:ENDO_HIST, ~ ifelse(.x == "", NA, .x))) %>%
  
  mutate(CER_LVSI = ifelse(CER_LVSI == 1, 'Yes', CER_LVSI)) %>%
  mutate(BMI = FDG_PW / (FDG_PH / 100)**2)

saveRDS(data_clean, "./data/processed_clin_data.rds")
write.csv(data_clean, './data/processed_clin_data.csv', row.names = F)

glimpse(data_clean)
table(data_clean$CER_HIST)
table(data_clean$ENDO_HIST)

table(data_clean$ENDO_LVSI)
table(data_clean$CER_LVSI)

table(data_clean$CER_DIFF)
table(data_clean$ENDO_DIFF)

table(data_clean$CER_FIGO)
table(data_clean$ENDO_FIGO)
