library(pacman)
p_load(tidyverse, infer, ggsignif, ggpubr)

data = read_csv('./data/quant_allPooled.csv')

#FDG_SUV in CER vs ENDO
g_fdg_suv = ggplot(data,
       aes(x = CANC,
           y = FDG_SUV)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  geom_signif(comparisons = list(c('cer', 'endo')),
              annotations = c(paste0('p = ', signif(wilcox.test(FDG_SUV ~ CANC, data=data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value,3), '*'))) +
  scale_x_discrete(label = c('Cervical', 'Endometrial')) +
  ylab('FDG SUVmax') +
  xlab('') +
  theme_classic() 

#FEC SUV in CER vs ENDO
g_fec_suv = ggplot(data,
                   aes(x = CANC,
                       y = FEC_SUV)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  geom_signif(comparisons = list(c('cer', 'endo')),
              annotations = c(paste0('p = ', signif(wilcox.test(FEC_SUV ~ CANC, data=data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value,3), ''))) +
  scale_x_discrete(label = c('Cervical', 'Endometrial')) +
  ylab('FEC SUVmax') +
  xlab('') +
  theme_classic() 

#ADC in CER vs ENDO
g_adc = ggplot(data,
       aes(x = CANC,
           y = ADC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  geom_signif(comparisons = list(c('cer', 'endo')),
              annotations = c(paste('p =', signif(wilcox.test(ADC ~ CANC, data=data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value,3)))) +
  scale_x_discrete(label = c('Cervical', 'Endometrial')) +
  ylab('ADCmean') +
  xlab('') +
  theme_classic() 

# can't group together SUVmax in CER + ENDO!

g = ggarrange(g_fdg_suv, g_fec_suv, g_adc, 
          ncol = 3,
          nrow = 1)

ggsave('./figures/FDG_FEC_ADC_byCancer.png',
       g,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)
