library(tidyverse)
library(ggpubr)
data = readRDS('./data/quant_allPooled.rds')
clin = read_csv('./data/processed_clin_data.csv')

data = data %>%
  inner_join(clin, by = 'PATIENT_ID') %>%
  mutate(PATIENT_ID = factor(PATIENT_ID))

labels = c('ADC' = 'ADCmean',
             'FDG_SUV' = 'FDG SUVmax',
             'FEC_SUV' = 'FEC SUVmax')

get_text = function(cor_res) {
  p = paste0('R = ', signif(cor_res$estimate, 2), ', p = ', signif(cor_res$p.value, 2))
  return(p)
}

# correlate nodal FDG, FEC, ADC with Age cervical ####
d = data %>%
  filter(CANC == 'cer') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, DEMO_AGE) 
d_plt = d %>% 
  pivot_longer(cols = c('FDG_SUV', 'FEC_SUV', 'ADC'))
cor_res = list(
  FDG_SUV = cor.test(~ FDG_SUV + DEMO_AGE, d),
  FEC_SUV = cor.test(~ FEC_SUV + DEMO_AGE, d),
  ADC = cor.test(~ ADC + DEMO_AGE, d)
)

txt = lapply(cor_res, get_text)
txt = data.frame(name = names(txt), 
                 value = unlist(txt))
txt = d_plt %>% 
  group_by(name) %>% 
  summarise(x = min(value, na.rm  = T)) %>%
  inner_join(txt, by = 'name')

g = ggplot(d_plt,
       aes(x = value,
           y = DEMO_AGE)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  geom_text(data = txt,
            aes(x = x,
                y = max(d_plt$DEMO_AGE, na.rm = T) + 5,
                label = value),
            hjust = 0) + 
  facet_wrap(~ name, scales = 'free_x', labeller = as_labeller(labels)) +
  ylab('Patient Age (yrs)') +
  xlab('') +
  theme_classic() 
ggsave('./figures/clinicaldata_figures/corr_age_cer.png',
       g,
       width = 6,
       height = 4,
       units = 'in',
       dpi = 300)

# correlate nodal FDG,FEC, ADC with Age endo  ####
d = data %>%
  filter(CANC == 'endo') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, DEMO_AGE) 
d_plt = d %>% 
  pivot_longer(cols = c('FDG_SUV', 'FEC_SUV', 'ADC'))
cor_res = list(
  FDG_SUV = cor.test(~ FDG_SUV + DEMO_AGE, d),
  FEC_SUV = cor.test(~ FEC_SUV + DEMO_AGE, d),
  ADC = cor.test(~ ADC + DEMO_AGE, d)
)

txt = lapply(cor_res, get_text)
txt = data.frame(name = names(txt), 
                 value = unlist(txt))
txt = d_plt %>% 
  group_by(name) %>% 
  summarise(x = min(value, na.rm  = T)) %>%
  inner_join(txt, by = 'name')

g = ggplot(d_plt,
           aes(x = value,
               y = DEMO_AGE)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  geom_text(data = txt,
            aes(x = x,
                y = max(d_plt$DEMO_AGE, na.rm = T) + 5,
                label = value),
            hjust = 0) + 
  facet_wrap(~ name, scales = 'free_x', labeller = as_labeller(labels)) +
  ylab('Patient Age (yrs)') +
  xlab('') +
  theme_classic() 
ggsave('./figures/clinicaldata_figures/corr_age_endo.png',
       g,
       width = 6,
       height = 4,
       units = 'in',
       dpi = 300)
# correlate BMI with FDG, FEC, ADC cervical ####
d = data %>%
  filter(CANC == 'cer') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, BMI) 
d_plt = d %>% 
  pivot_longer(cols = c('FDG_SUV', 'FEC_SUV', 'ADC'))
cor_res = list(
  FDG_SUV = cor.test(~ FDG_SUV + BMI, d),
  FEC_SUV = cor.test(~ FEC_SUV + BMI, d),
  ADC = cor.test(~ ADC + BMI, d)
)

txt = lapply(cor_res, get_text)
txt = data.frame(name = names(txt), 
                 value = unlist(txt))
txt = d_plt %>% 
  group_by(name) %>% 
  summarise(x = min(value, na.rm  = T)) %>%
  inner_join(txt, by = 'name')

g = ggplot(d_plt,
           aes(x = value,
               y = BMI)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  geom_text(data = txt,
            aes(x = x,
                y = max(d_plt$BMI, na.rm = T) + 5,
                label = value),
            hjust = 0) + 
  facet_wrap(~ name, scales = 'free_x', labeller = as_labeller(labels)) +
  ylab(expression(BMI~(kg/m^2))) +
  xlab('') +
  theme_classic() 
ggsave('./figures/clinicaldata_figures/corr_bmi_cer.png',
       g,
       width = 6,
       height = 4,
       units = 'in',
       dpi = 300)

# correlate BMI with FDG, FEC, ADC endo ####
d = data %>%
  filter(CANC == 'endo') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, BMI) 
d_plt = d %>% 
  pivot_longer(cols = c('FDG_SUV', 'FEC_SUV', 'ADC'))
cor_res = list(
  FDG_SUV = cor.test(~ FDG_SUV + BMI, d),
  FEC_SUV = cor.test(~ FEC_SUV + BMI, d),
  ADC = cor.test(~ ADC + BMI, d)
)

txt = lapply(cor_res, get_text)
txt = data.frame(name = names(txt), 
                 value = unlist(txt))
txt = d_plt %>% 
  group_by(name) %>% 
  summarise(x = min(value, na.rm  = T)) %>%
  inner_join(txt, by = 'name')

g = ggplot(d_plt,
           aes(x = value,
               y = BMI)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  geom_text(data = txt,
            aes(x = x,
                y = max(d_plt$BMI, na.rm = T) + 5,
                label = value),
            hjust = 0) + 
  facet_wrap(~ name, scales = 'free_x', labeller = as_labeller(labels)) +
  ylab(expression(BMI~(kg/m^2))) +
  xlab('') +
  theme_classic() 
ggsave('./figures/clinicaldata_figures/corr_bmi_endo.png',
       g,
       width = 6,
       height = 4,
       units = 'in',
       dpi = 300)
# boxplot nodal FDG, FEC, ADC with FIGO cervical ####
d = data %>%
  filter(CANC == 'cer') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, CER_FIGO, HIST) 
g1 = ggplot(d,
            aes(x = CER_FIGO,
                y = FDG_SUV,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('FDG SUVmax')
g2 = ggplot(d,
            aes(x = CER_FIGO,
                y = FEC_SUV,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('FEC SUVmax')
g3 = ggplot(d,
            aes(x = CER_FIGO,
                y = ADC,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(name = 'Histology', labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('ADCmean')
g = ggarrange(g1,g2,g3, 
          nrow = 1, 
          common.legend = T) 
ggsave('./figures/clinicaldata_figures/boxplot_figo_cer.png',
       g,
       width = 8,
       height = 4,
       units = 'in',
       dpi = 300,
       bg = 'white')

# boxplot nodal FDG, FEC, ADC with FIGO endometrial ####
d = data %>%
  filter(CANC == 'endo') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, ENDO_FIGO, HIST) %>%
  filter(ENDO_FIGO != 'Unknown')
g1 = ggplot(d,
            aes(x = ENDO_FIGO,
                y = FDG_SUV,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('FDG SUVmax')
g2 = ggplot(d,
            aes(x = ENDO_FIGO,
                y = FEC_SUV,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('FEC SUVmax')
g3 = ggplot(d,
            aes(x = ENDO_FIGO,
                y = ADC,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(name = 'Histology', labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('ADCmean')
g = ggarrange(g1,g2,g3, 
              nrow = 1, 
              common.legend = T) 
ggsave('./figures/clinicaldata_figures/boxplot_figo_endo.png',
       g,
       width = 8,
       height = 4,
       units = 'in',
       dpi = 300,
       bg = 'white')

# boxplot nodal FDG, FEC, ADC with cervical histology ####
d = data %>%
  filter(CANC == 'cer') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, CER_HIST, HIST) %>%
  filter(!is.na(CER_HIST),
         CER_HIST != 'Unspecified')
g1 = ggplot(d,
            aes(x = CER_HIST,
                y = FDG_SUV,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('FDG SUVmax') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g2 = ggplot(d,
            aes(x = CER_HIST,
                y = FEC_SUV,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('FEC SUVmax') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g3 = ggplot(d,
            aes(x = CER_HIST,
                y = ADC,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(name = 'Histology', labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('ADCmean') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g = ggarrange(g1,g2,g3, 
              nrow = 1, 
              common.legend = T) 
ggsave('./figures/clinicaldata_figures/boxplot_histology_cer.png',
       g,
       width = 8,
       height = 4,
       units = 'in',
       dpi = 300,
       bg = 'white')
# boxplot nodal FDG, FEC, ADC with endo histology ####
d = data %>%
  filter(CANC == 'endo') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, ENDO_HIST, HIST) %>%
  filter(!is.na(ENDO_HIST))
g1 = ggplot(d,
            aes(x = ENDO_HIST,
                y = FDG_SUV,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('FDG SUVmax')
g2 = ggplot(d,
            aes(x = ENDO_HIST,
                y = FEC_SUV,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('FEC SUVmax')
g3 = ggplot(d,
            aes(x = ENDO_HIST,
                y = ADC,
                color = HIST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1),
             alpha = .5) +
  theme_classic() +
  scale_color_discrete(name = 'Histology', labels = c('Benign', 'Malignant')) +
  xlab('') +
  ylab('ADCmean')
g = ggarrange(g1,g2,g3, 
              nrow = 1, 
              common.legend = T) 
ggsave('./figures/clinicaldata_figures/boxplot_histology_endo.png',
       g,
       width = 8,
       height = 4,
       units = 'in',
       dpi = 300,
       bg = 'white')
