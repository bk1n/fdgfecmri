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
get_p = function(figo_grp) {
  #fdg
  x = figo_grp[figo_grp$HIST == 1,]$FDG_SUV
  y = figo_grp[figo_grp$HIST == 0,]$FDG_SUV
  fdg = tryCatch({
    t.test(x,y)$p.value
  }, error = function(e){
    fdg = NA
  })
  #fec
  x = figo_grp[figo_grp$HIST == 1,]$FEC_SUV
  y = figo_grp[figo_grp$HIST == 0,]$FEC_SUV
  fec = tryCatch({
    t.test(x,y)$p.value
  }, error = function(e){
    fec = NA
  })
  #adc
  x = figo_grp[figo_grp$HIST == 1,]$ADC
  y = figo_grp[figo_grp$HIST == 0,]$ADC
  adc = tryCatch({
    t.test(x,y)$p.value
  }, error = function(e){
    adc = NA
  })
  return(list(
    fdg = fdg,
    fec = fec,
    adc = adc)
  )
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

p_val = unlist(lapply(cor_res, function(c) c$p.value))
p_adj = p.adjust(p_val, method = 'BH')
for(i in 1:3){
  cor_res[[i]]$p.value = p_adj[i]
}

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
  select(PATIENT_ID, HIST, FDG_SUV, FEC_SUV, ADC, BMI) 
d_plt = d %>% 
  pivot_longer(cols = c('FDG_SUV', 'FEC_SUV', 'ADC'))
cor_res = list(
  FDG_SUV_benign = cor.test(~ FDG_SUV + BMI, d %>% filter(HIST == 0)),
  FEC_SUV_benign = cor.test(~ FEC_SUV + BMI, d %>% filter(HIST == 0)),
  ADC_benign = cor.test(~ ADC + BMI, d %>% filter(HIST == 0)),
  FDG_SUV_mal = cor.test(~ FDG_SUV + BMI, d %>% filter(HIST == 1)),
  FEC_SUV_mal = cor.test(~ FEC_SUV + BMI, d %>% filter(HIST == 1)),
  ADC_mal = cor.test(~ ADC + BMI, d %>% filter(HIST == 1))
)
cor_res = lapply(cor_res, function(c) {
  c$p.value = p.adjust(c$p.value, method = 'bonferroni', n = 6)
  return(c)
})
  
txt = lapply(cor_res, get_text)
txt = data.frame(name = names(txt), 
                 value = unlist(txt))
txt = txt %>%
  mutate(HIST = if_else(grepl('_benign', name), 0, 1)) %>%
  mutate(name = gsub('_[^_]+$', '', name)) %>%
  mutate(HIST = factor(HIST)) %>%
  rename(cor_txt = value)

txt = d_plt %>% 
  group_by(name, HIST) %>% 
  inner_join(txt, by = c('name' = 'name', 'HIST' = 'HIST')) %>%
  group_by(name) %>%
  mutate(x = min(value, na.rm = T),
         y = max(d_plt$BMI, na.rm = T) + 5,
         y = if_else(HIST == 0, y * 1.03, y)) %>%
  distinct(name, HIST, cor_txt, x, y)

g = ggplot(d_plt,
           aes(x = value,
               y = BMI)) +
  geom_point(aes(color = HIST),
             alpha = .7) +
  geom_smooth(aes(color = HIST),
              method = 'lm') +
  geom_text(data = txt,
            aes(x = x,
                y = y,
                label = cor_txt,
                color = HIST),
            hjust = 0) + 
  facet_wrap(~ name, scales = 'free_x', labeller = as_labeller(labels)) +
  ylab(expression(BMI~(kg/m^2))) +
  xlab('') +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant'))
ggsave('./figures/clinicaldata_figures/corr_bmi_cer.png',
       g,
       width = 8,
       height = 5,
       units = 'in',
       dpi = 300)

# correlate BMI with FDG, FEC, ADC endo ####
d = data %>%
  filter(CANC == 'endo') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, BMI, HIST)
d_plt = d %>% 
  pivot_longer(cols = c('FDG_SUV', 'FEC_SUV', 'ADC')) 
cor_res = list(
  FDG_SUV_benign = cor.test(~ FDG_SUV + BMI, d %>% filter(HIST == 0)),
  FEC_SUV_benign = cor.test(~ FEC_SUV + BMI, d %>% filter(HIST == 0)),
  ADC_benign = cor.test(~ ADC + BMI, d %>% filter(HIST == 0)),
  FDG_SUV_mal = cor.test(~ FDG_SUV + BMI, d %>% filter(HIST == 1)),
  FEC_SUV_mal = cor.test(~ FEC_SUV + BMI, d %>% filter(HIST == 1)),
  ADC_mal = cor.test(~ ADC + BMI, d %>% filter(HIST == 1))
)
cor_res = lapply(cor_res, function(c) {
  c$p.value = p.adjust(c$p.value, method = 'bonferroni', n = 6)
  return(c)
})

txt = lapply(cor_res, get_text)
txt = data.frame(name = names(txt), 
                 value = unlist(txt))
txt = txt %>%
  mutate(HIST = if_else(grepl('_benign', name), 0, 1)) %>%
  mutate(name = gsub('_[^_]+$', '', name)) %>%
  mutate(HIST = factor(HIST)) %>%
  rename(cor_txt = value)

txt = d_plt %>% 
  group_by(name, HIST) %>% 
  inner_join(txt, by = c('name' = 'name', 'HIST' = 'HIST')) %>%
  group_by(name) %>%
  mutate(x = min(value, na.rm = T),
         y = max(d_plt$BMI, na.rm = T) + 5,
         y = if_else(HIST == 0, y * 1.03, y)) %>%
  distinct(name, HIST, cor_txt, x, y)

g = ggplot(d_plt,
           aes(x = value,
               y = BMI)) +
  geom_point(aes(color = HIST),
             alpha = .7) +
  geom_smooth(aes(color = HIST),
              method = 'lm') +
  geom_text(data = txt,
            aes(x = x,
                y = y,
                label = cor_txt,
                color = HIST),
            hjust = 0) + 
  facet_wrap(~ name, scales = 'free_x', labeller = as_labeller(labels)) +
  ylab(expression(BMI~(kg/m^2))) +
  xlab('') +
  theme_classic() +
  scale_color_discrete(labels = c('Benign', 'Malignant'))
ggsave('./figures/clinicaldata_figures/corr_bmi_endo.png',
       g,
       width = 8,
       height = 5,
       units = 'in',
       dpi = 300)
# boxplot nodal FDG, FEC, ADC with FIGO cervical ####
d = data %>%
  filter(CANC == 'cer') %>%
  select(PATIENT_ID, FDG_SUV, FEC_SUV, ADC, CER_FIGO, HIST) 
# get pvals
figo = group_split(d, CER_FIGO)
p_cer = lapply(figo, get_p)
names(p_cer) = sort(unique(d$CER_FIGO))

cer_pval = do.call(rbind, p_cer) %>%
  as.data.frame %>%
  rownames_to_column('figo') %>%
  mutate(across(fdg:adc, unlist)) %>%
  pivot_longer(cols = fdg:adc) %>%
  mutate(value = p.adjust(value, method = 'fdr'))

pvals = signif(cer_pval[cer_pval$name == 'fdg',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('FDG SUVmax') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$FDG_SUV, na.rm = T) + 3,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

pvals = signif(cer_pval[cer_pval$name == 'fec',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('FEC SUVmax') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$FEC_SUV, na.rm = T) + 1,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

pvals = signif(cer_pval[cer_pval$name == 'adc',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('ADCmean') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$ADC, na.rm = T) + 70,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

g = ggarrange(g1, g2, g3, 
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
# get pvals
figo = group_split(d, ENDO_FIGO)
p_endo = lapply(figo, get_p)
names(p_endo) = sort(unique(d$ENDO_FIGO))

endo_pval = do.call(rbind, p_endo) %>%
  as.data.frame %>%
  rownames_to_column('figo') %>%
  mutate(across(fdg:adc, unlist)) %>%
  pivot_longer(cols = fdg:adc) %>%
  mutate(value = p.adjust(value, method = 'fdr'))

pvals = signif(endo_pval[endo_pval$name == 'fdg',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('FDG SUVmax') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$FDG_SUV, na.rm = T) + 3,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

pvals = signif(endo_pval[endo_pval$name == 'fec',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('FEC SUVmax') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$FEC_SUV, na.rm = T) + 1,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

pvals = signif(endo_pval[endo_pval$name == 'adc',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('ADCmean') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$ADC, na.rm = T) + 40,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')
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

# get pvals
figo = group_split(d, CER_HIST)
p_cer = lapply(figo, get_p)
names(p_cer) = sort(unique(d$CER_HIST))

cer_pval = do.call(rbind, p_cer) %>%
  as.data.frame %>%
  rownames_to_column('figo') %>%
  mutate(across(fdg:adc, unlist)) %>%
  pivot_longer(cols = fdg:adc) %>%
  mutate(value = p.adjust(value, method = 'fdr'))

pvals = signif(cer_pval[cer_pval$name == 'fdg',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$FDG_SUV, na.rm = T) + 3,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

pvals = signif(cer_pval[cer_pval$name == 'fec',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$FEC_SUV, na.rm = T) + 1,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

pvals = signif(cer_pval[cer_pval$name == 'adc',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$ADC, na.rm = T) + 40,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

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
# get pvals
figo = group_split(d, ENDO_HIST)
p_cer = lapply(figo, get_p)
names(p_cer) = sort(unique(d$ENDO_HIST))

cer_pval = do.call(rbind, p_cer) %>%
  as.data.frame %>%
  rownames_to_column('figo') %>%
  mutate(across(fdg:adc, unlist)) %>%
  pivot_longer(cols = fdg:adc) %>%
  mutate(value = p.adjust(value, method = 'fdr'))

pvals = signif(cer_pval[cer_pval$name == 'fdg',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('FDG SUVmax') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$FDG_SUV, na.rm = T) + 3,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

pvals = signif(cer_pval[cer_pval$name == 'fec',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('FEC SUVmax') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$FEC_SUV, na.rm = T) + 1,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

pvals = signif(cer_pval[cer_pval$name == 'adc',]$value, 2)
x = seq(1, length(pvals))
x_min = x - .2
x_max = x + .2
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
  ylab('ADCmean') +
  geom_signif(
    color = 'black',
    annotations = pvals,
    xmin = x_min,
    xmax = x_max,
    y_position = max(d$ADC, na.rm = T) + 40,
    angle = 45,
    hjust = 0
  ) +
  coord_cartesian(clip = 'off')

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

