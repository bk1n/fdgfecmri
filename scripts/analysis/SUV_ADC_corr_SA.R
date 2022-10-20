library(pacman)
p_load(tidyverse, infer)

fdg = readRDS('./data/fdgpet_allPooled.rds')

g_suv = ggplot(fdg,
       aes(x = SA,
           y = SUV)) +
  geom_point(position = position_jitter(width = 0.3)) +
  geom_smooth(method = 'lm',
              formula = y~x,
              se = T) +
  stat_cor(method = 'pearson') +
  ylab('SUVmax') +
  xlab('Short axis diameter (mm)') +
  theme_classic() 
  
g_adc = ggplot(fdg,
       aes(x = SA,
           y = ADC)) +
  geom_point() +
  geom_smooth(method = 'lm',
              formula = y~x,
              se = T) +
  stat_cor(method = 'pearson') +
  xlim(c(0,20)) +
  ylab('ADCmean') +
  xlab('Short axis diameter (mm)') +
  theme_classic() 

g = ggarrange(g_suv, g_adc,
          nrow = 1,
          ncol = 2)

ggsave('./figures/SUV_ADC_corr_SA.png',
       g,
       height = 5,
       width = 5,
       units = 'in',
       dpi = 300)
