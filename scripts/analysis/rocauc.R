library(pacman)
p_load(tidyverse, cutpointr)

df = readRDS('./data/fdgpet_allPooled.rds') 

df_cer = df %>% filter(CANC == 'cer')
df_endo = df %>% filter(CANC == 'endo')

get_opt_cut = function(measure){
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
logit_res = read.csv('./figures/tables/lr.model_fullRes.csv')

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
  ylab('Area Under Curve (AUC)')
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
  ylab('Sensitivity')
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
  ylab('Specificity')
ggsave('./figures/Spec_byMeasure.png',
       g_spec,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

library(ggpubr)
g = ggarrange(g_auc, g_sens, g_spec, nrow = 3)

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
lr = read.csv('./figures/tables/lr_roc_results.csv')
roc_plt_df = rbind(roc_plt_df, lr)

g_roc = ggplot(roc_plt_df,
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
  xlab('Specificity') +
  ylab('Specificity') +
  scale_color_discrete(name = 'Measure')

ggsave('./figures/ROC_byMeasure.png',
       g_roc,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

