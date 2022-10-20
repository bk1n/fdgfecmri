library(pacman)
p_load(tidyverse, infer, ggsignif)

fdg = readRDS('./data/fdgpet_allPooled.rds')

#SUV in CER vs ENDO
g_suv = ggplot(fdg,
       aes(x = CANC,
           y = SUV)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  geom_signif(comparisons = list(c('cer', 'endo')),
              annotations = c(paste0('p = ', signif(t_test(fdg, SUV ~ CANC)$p_value,3), '***'))) +
  scale_x_discrete(label = c('Cervical', 'Endometrial')) +
  ylab('SUVmax') +
  xlab('') +
  theme_classic() 

#ADC in CER vs ENDO
g_adc = ggplot(fdg,
       aes(x = CANC,
           y = ADC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  geom_signif(comparisons = list(c('cer', 'endo')),
              annotations = c(paste('p =', signif(t_test(fdg, ADC ~ CANC)$p_value,3)))) +
  scale_x_discrete(label = c('Cervical', 'Endometrial')) +
  ylab('ADCmean') +
  xlab('') +
  theme_classic() 

# can't group together SUVmax in CER + ENDO!

g = ggarrange(g_suv, g_adc, 
          ncol = 2,
          nrow = 1)

ggsave('./figures/SUV_ADC_byCancer.png',
       g,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)
