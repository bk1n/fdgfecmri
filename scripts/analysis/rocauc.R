library(pacman)
p_load(tidyverse, cutpointr)

df = readRDS('./data/quant_allPooled.rds')

df = df %>%
  select(!contains('PT'))

df_cer = df %>% filter(CANC == 'cer')
df_endo = df %>% filter(CANC == 'endo')

get_opt_cut = function(measure){
  print(paste('Getting cutpoint for:', measure))
  x = pull(df[,measure])
  class = df$HIST
  opt_cut = cutpointr(df_endo,
          x = x,
          class = class,
          na.rm = T,
          method = maximize_metric,
          metric = F1_score,
          pos_class = 1,
          neg_class = 0, 
          direction = '>=',
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
logit_res = read.csv('./figures/tables/lr.model_optCutOff_fullRes_SUV.csv')

plt_logit = logit_res %>%
  select(optimal_SUVcutoff, auc, sens, spec) %>%
  rename(optimal_cutpoint = optimal_SUVcutoff)
plt_logit$name = 'LR'

plt_df_full = plt_df %>%
  select(optimal_cutpoint, AUC_oob, sensitivity_oob, specificity_oob, name) %>%
  rename(auc = AUC_oob,
         sens = sensitivity_oob,
         spec = specificity_oob) %>%
  rbind(., plt_logit)

dp_df = plt_df_full %>%
  group_by(name) %>%
  summarise(across(everything(), ~mean(.x, na.rm = T))) %>%
  arrange(sens)

write.csv(dp_df, './figures/optimalCP_byMeasure.csv', row.names = F)

dp_df = dp_df %>%
  arrange(name) %>%
  mutate(name = factor(name, levels = name))

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g_spec
ggsave('./figures/Spec_byMeasure.png',
       g_spec,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

library(ggpubr)
g = ggarrange(g_auc, 
              g_sens , 
              g_spec, nrow = 3)

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
lr = read.csv('./figures/tables/lr_roc_results_SUV.csv')
roc_plt_df = rbind(roc_plt_df, lr)

roc_plt = roc_plt_df %>%
  mutate(measure_type = case_when(grepl('FDG', name) ~ 'FDG',
                                  grepl('FEC', name) ~ 'FEC',
                                  grepl('ADC', name) ~ 'ADC',
                                  name == 'LR' ~ 'LR')) %>%
  mutate(name = gsub('FDG_', '', name),
         name = gsub('FEC_', '', name),
         name = gsub('ADC_', '', name)) %>%
  mutate(measure_type = factor(measure_type, levels = c('FDG', 'FEC', 'ADC', 'LR'))) %>%
  mutate(name = factor(name, levels = c('SUV', 'STAR', 'SNSA', 'NTR', 'SA', 'LA', 'ADC', 'LR')))

g_roc = ggplot(roc_plt,
       aes(x = fpr,
           y = 1 - fnr,
           color = name)) +
  geom_line(alpha = 0.8) +
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
  scale_color_discrete(name = 'Quantitative \nMeasure') +
  scale_linetype_manual(name = 'Measure Type',
                        values = c('solid', 'twodash', 'dotdash', 'dashed')) +
  facet_wrap(~ measure_type)

ggsave('./figures/ROC_byMeasure.png',
       g_roc,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

