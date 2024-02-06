#### FDG, FEC, ADC by Cancer
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
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

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
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

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
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# can't group together SUVmax in CER + ENDO!

g_byCancer = ggarrange(g_fdg_suv, g_fec_suv, g_adc, 
          ncol = 3,
          nrow = 1)

ggsave('./figures/FDG_FEC_ADC_byCancer.png',
       g,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

#### FDG, FEC, ADC by Region ####
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

g_byRegion = ggarrange(g_fdg_suv, g_fec_suv, g_adc,
              nrow = 1,
              ncol = 3)

ggsave('./figures/FDG_FEC_ADC_byRegion.png',
       g,
       height = 5,
       width = 5,
       units = 'in',
       dpi = 300)

#### by Region, by Cancer ####
g = ggarrange(g_byRegion, g_byCancer,
          nrow = 1,
          labels = c('A', 'B'))
ggsave('./figures/FDG_FEC_ADC_byRegion_byCancer_combined.png',
       g, 
       height = 5,
       width = 10,
       units = 'in',
       dpi = 300)
