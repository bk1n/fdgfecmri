library(pacman)
p_load(tidyverse, cutpointr)

df = readRDS('./data/quant_allPooled.rds')

df = df %>%
  select(!contains('PT'))

df_cer = df %>% filter(CANC == 'cer')
df_endo = df %>% filter(CANC == 'endo')

get_opt_cut = function(measure){
  print(paste('Getting cutpoint for:', measure))
  x = pull(df_endo[,measure])
  class = df_endo$HIST
  if(measure == 'ADC' | measure == 'ADC_NTR'){
    dir = '<='
  } else {
    dir = '>='
  }
  opt_cut = cutpointr(df_endo,
          x = x,
          class = class,
          na.rm = T,
          method = maximize_metric,
          metric = F1_score,
          pos_class = 1,
          neg_class = 0, 
          direction = dir,
          boot_runs = 2000,
          boot_stratify = T)
  
  summary.opt_cut = c(measure, summary(opt_cut))
  
  return(summary.opt_cut)
}

cols_endo = df_endo %>%
  select(-c(CANC, HIST, region)) %>%
  colnames

res.opt_cut = lapply(cols_endo, get_opt_cut)

get_boot_data = function(summ){
  bd = summ$cutpointr[[1]]$boot[[1]]
  bd$name = summ[[1]]
  return(bd)
}

plt_df = do.call(rbind, lapply(res.opt_cut, get_boot_data))
logit_suv_res = read.csv('./figures/tables/lr.model_optCutOff_fullRes_SUV.csv')
logit_sa_res = read.csv('./figures/tables/lr.model_optCutOff_fullRes_SA.csv')

plt_logit_suv = logit_suv_res %>%
  select(optimal_SUVcutoff, auc, sens, spec) %>%
  rename(optimal_cutpoint = optimal_SUVcutoff)
plt_logit_suv$name = 'LR_SUV'

plt_logit_sa = logit_sa_res %>%
  select(optimal_SUVcutoff, auc, sens, spec) %>%
  rename(optimal_cutpoint = optimal_SUVcutoff)
plt_logit_sa$name = 'LR_SA'

plt_df_full = plt_df %>%
  select(optimal_cutpoint, AUC_oob, sensitivity_oob, specificity_oob, name) %>%
  rename(auc = AUC_oob,
         sens = sensitivity_oob,
         spec = specificity_oob) %>%
  rbind(., plt_logit_suv, plt_logit_sa) %>%
  mutate()

dp_df = plt_df_full %>%
  group_by(name) %>%
  summarise(across(everything(), ~mean(.x, na.rm = T))) %>%
  arrange(sens)

write.csv(plt_df_full, './figures/tables/full_roc_results_plt_df_full.csv', row.names = F)
write.csv(dp_df, './figures/tables/optimalCP_byMeasure.csv', row.names = F)

#### RUN FROM HERE

ylabs = c('FDG_SUV_PT' = 'FDG SUVmax',
          'FEC_SUV_PT' = 'FEC SUVmax',
          'ADC_PT' = 'PT ADCmean',
          'FDG_SUV' = 'FDG SUVmax',
          'FDG_SA' = 'FDG SA (mm)',
          'FDG_LA' = 'FDG LA (mm)',
          'FDG_NTR' = 'FDG SUVmax NTR',
          'FDG_STAR' = 'FDG STAR',
          'FDG_SNSA' = 'FDG SNSA',
          'FEC_SUV' = 'FEC SUVmax',
          'FEC_SA' = 'FEC SA (mm)',
          'FEC_LA' = 'FEC LA (mm)',
          'FEC_NTR' = 'FEC SUVmax NTR',
          'FEC_STAR' = 'FEC STAR',
          'FEC_SNSA' = 'FEC SNSA',
          'ADC' = 'ADCmean',
          'ADC_NTR' = 'ADCmean NTR',
          'LR_SUV' = 'LR (SUV)',
          'LR_SA' = 'LR (SA)')

plt_df_full = read.csv('./figures/tables/full_roc_results_plt_df_full.csv') %>%
  select(-1) 
dp_df = read.csv('./figures/tables/optimalCP_byMeasure.csv')

dp_df = dp_df %>%
  arrange(name) %>%
  mutate(name = ylabs[name]) %>%
  mutate(name = factor(name, levels = name)) %>%
  mutate(measure_type = case_when(grepl('FDG', name) ~ 'FDG-PET/CT',
                                  grepl('FEC', name) ~ 'FEC-PET/CT',
                                  grepl('ADC', name) ~ 'DW-MRI',
                                  grepl('LR', name) ~ 'LR'))

plt_df_full = plt_df_full %>%
  mutate(measure_type = case_when(grepl('FDG', name) ~ 'FDG-PET/CT',
                                  grepl('FEC', name) ~ 'FEC-PET/CT',
                                  grepl('ADC', name) ~ 'DW-MRI',
                                  grepl('LR', name) ~ 'LR')) %>%
  mutate(name = ylabs[name]) 

g_opt = ggplot(dp_df,
       aes(x = name,
           y = '   ')) +
  geom_text(aes(label = signif(optimal_cutpoint, digits = 3)), angle = 45)  +
  theme_classic() +
  xlab('') +
  ylab('Optimal \nCutpoint') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~ measure_type, space = 'free', scales = 'free_x')

g_auc = ggplot(plt_df_full,
       aes(x = name,
           y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = dp_df,
             size = 3,
             shape = 3,
             color = 'red') +
  theme_classic() +
  xlab('') +
  ylab('Area Under Curve (AUC)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x = element_blank()) +
  facet_grid(~ measure_type, space = 'free', scales = 'free_x')
ggsave('./figures/AUC_byMeasure.png',
       g_auc,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300) 

g_sens = ggplot(plt_df_full,
       aes(x = name,
           y = sens)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = dp_df,
             size = 3,
             shape = 3,
             color = 'red') +
  theme_classic() +
  xlab('') +
  ylab('Sensitivity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x = element_blank()) +
  facet_grid(~ measure_type, space = 'free', scales = 'free_x')

ggsave('./figures/Sens_byMeasure.png',
       g_sens,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

g_spec = ggplot(plt_df_full,
       aes(x = name,
           y = spec)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = dp_df,
             size = 3,
             shape = 3,
             color = 'red') +
  theme_classic() +
  xlab('Measure') +
  ylab('Specificity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~ measure_type, space = 'free', scales = 'free_x')

g_spec
ggsave('./figures/Spec_byMeasure.png',
       g_spec,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

library(ggpubr)
g = ggarrange(g_auc, 
              g_sens, 
              g_spec, nrow = 3, heights = c(0.8, 0.8, 1))
g

ggsave('./figures/AllPerformance_byMeasure.png',
       g,
       width = 5,
       height = 9,
       units = 'in',
       dpi = 300)


# plot ROC curves based on full dataset
library(pROC)

get_proc = function(measure){
  proc = cutpointr::roc(df_endo, !!measure, HIST, pos_class = 1, neg_class = 0)
  proc$name = measure
  return(proc)
}

roc_plt_df = do.call(rbind, lapply(cols_endo, get_proc))
lr_suv = read.csv('./figures/tables/lr_roc_results_SUV.csv') %>% mutate(name = 'LR (SUV)')
lr_sa = read.csv('./figures/tables/lr_roc_results_SA.csv') %>% mutate(name = 'LR (SA)')
roc_plt_df = rbind(roc_plt_df, lr_suv, lr_sa)

roc_plt = roc_plt_df %>%
  mutate(measure_type = case_when(grepl('FDG', name) ~ 'FDG-PET/CT',
                                  grepl('FEC', name) ~ 'FEC-PET/CT',
                                  grepl('ADC', name) ~ 'DW-MRI',
                                  name == 'LR' ~ 'LR')) %>%
  mutate(name = gsub('FDG_', '', name),
         name = gsub('FEC_', '', name),
         name = gsub('ADC_', '', name)) %>%
  mutate(name = case_when(name == 'SUV' ~'SUVmax', 
                          name == 'SA' ~ 'SA (mm)',
                          name == 'LA' ~ 'LA (mm)',
                          name == 'ADC' ~ 'ADCmean',
                          .default = name)) %>%
  mutate(measure_type = factor(measure_type, levels = c('FDG-PET/CT', 'FEC-PET/CT', 'DW-MRI', 'LR'))) %>%
  mutate(name = factor(name, levels = c('SUVmax', 'STAR', 'SNSA', 'NTR', 'SA (mm)', 'LA (mm)', 'ADCmean', 'LR (SA)', 'LR (SUV)')))

g_roc = ggplot(roc_plt,
       aes(x = fpr,
           y = 1 - fnr,
           color = name)) +
  geom_line(alpha = 0.7,
            linewidth = 1) +
  geom_line(inherit.aes = F,
            data = data.frame(x = seq(0,1), y = seq(0,1)),
            aes(x=x,
                y=y),
            linetype = 'dashed',
            color = 'red',
            alpha = 0.5) + 
  theme_classic() +
  xlab('1 - Specificity') +
  ylab('Sensitivity') +
  scale_color_brewer(name = 'Quantitative \nMeasure',
                     palette = 'Set1',
                     direction = -1) +
  scale_linetype_manual(name = 'Measure Type',
                        values = c('solid', 'twodash', 'dotdash', 'dashed')) +
  facet_wrap(~ measure_type)
g_roc

ggsave('./figures/ROC_byMeasure.png',
       g_roc,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

