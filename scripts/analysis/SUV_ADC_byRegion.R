library(pacman)
p_load(tidyverse, infer)

fdg = readRDS('./data/fdgpet_allPooled.rds')

# SUV in LP, RP, PALN in CER
res.aov = aov(SUV ~ region, fdg)
summary(res.aov)
TukeyHSD(res.aov)

g_suv = ggplot(fdg,
       aes(x = region,
           y = SUV)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  scale_x_discrete(label = c('Left Pelvis', 'Para-aortic', 'Right Pelvis')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab('SUVmax') +
  xlab('') 

# ADC in LP, RP, PALN
res.aov = aov(ADC ~ region, fdg)
summary(res.aov)
TukeyHSD(res.aov)

g_adc = ggplot(fdg,
       aes(x = region,
           y = ADC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  scale_x_discrete(label = c('Left Pelvis', 'Para-aortic', 'Right Pelvis')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab('ADCmean') +
  xlab('') 

g = ggarrange(g_suv, g_adc,
          nrow = 1,
          ncol = 2)

ggsave('./figures/SUV_ADC_byRegion.png',
       g,
       height = 5,
       width = 5,
       units = 'in',
       dpi = 300)
