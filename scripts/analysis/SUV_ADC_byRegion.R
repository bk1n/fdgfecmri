library(pacman)
p_load(tidyverse, infer, ggpubr)

data = readRDS('./data/quant_allPooled.rds')

# FDG SUV in LP, RP, PALN in CER
res.aov = aov(FDG_SUV ~ region, data)
summary(res.aov)
TukeyHSD(res.aov)

g_fdg_suv = ggplot(data,
       aes(x = region,
           y = FDG_SUV)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  scale_x_discrete(label = c('Left Pelvis', 'Para-aortic', 'Right Pelvis')) +
  geom_signif(
    comparisons = list(c('paln', 'lp'), 
                       c('rp', 'lp'),
                       c('rp', 'paln')),
    annotations = c(paste('p =', signif(TukeyHSD(res.aov)$region[1:3,'p adj'], digits = 2))),
    y_position = c(30,38,34)
  ) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  ylab('FDG SUVmax') +
  xlab('') 

# FEC SUV in LP, RP, PALN in CER
res.aov = aov(FEC_SUV ~ region, data)
summary(res.aov)
TukeyHSD(res.aov)

g_fec_suv = ggplot(data,
                   aes(x = region,
                       y = FEC_SUV)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  scale_x_discrete(label = c('Left Pelvis', 'Para-aortic', 'Right Pelvis')) +
  geom_signif(
    comparisons = list(c('paln', 'lp'), 
                       c('rp', 'lp'),
                       c('rp', 'paln')),
    annotations = c(paste('p =', signif(TukeyHSD(res.aov)$region[1:3,'p adj'], digits = 2))),
    y_position = c(9,11,10)
    ) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size = 10)) +
  ylab('FEC SUVmax') +
  xlab('') 

# ADC in LP, RP, PALN
res.aov = aov(ADC ~ region, data)
summary(res.aov)
TukeyHSD(res.aov)

g_adc = ggplot(data,
       aes(x = region,
           y = ADC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  scale_x_discrete(label = c('Left Pelvis', 'Para-aortic', 'Right Pelvis')) +
  geom_signif(
    comparisons = list(c('paln', 'lp'), 
                       c('rp', 'lp'),
                       c('rp', 'paln')),
    annotations = c(paste('p =', signif(TukeyHSD(res.aov)$region[1:3,'p adj'], digits = 2))),
    y_position = c(2000,2400,2200)
  ) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  ylab('ADCmean') +
  xlab('') 

g = ggarrange(g_fdg_suv, g_fec_suv, g_adc,
          nrow = 1,
          ncol = 3)
g

ggsave('./figures/FDG_FEC_ADC_byRegion.png',
       g,
       height = 5,
       width = 5,
       units = 'in',
       dpi = 300)
