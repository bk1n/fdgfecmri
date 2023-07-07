library(pacman)
p_load(tidyverse, infer, Hmisc, ComplexHeatmap)

fdg =  read_csv('./data/quant_allPooled.csv')

d = fdg %>%
  filter(CANC == 'endo') %>%
  select(-c(CANC, region)) %>%
  select(-contains(c('STAR', 'NTR', 'SNSA'))) %>%
  select(-HIST)

pairwise_corr = rcorr(as.matrix(d))
sig = pairwise_corr$P
sig[is.na(sig)] = 0
sig = p.adjust(sig, method='fdr') %>% matrix(ncol=length(colnames(pairwise_corr$P)))

png('./figures/pairwise_corr_raw.png',
    width = 6,
    height = 5.5,
    units = 'in',
    res = 300)
Heatmap(pairwise_corr$r, col = c('white', 'red'),
        heatmap_legend_param = list(title = 'R'),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(sig[i, j] < 0.001) {
            grid.text("***", x, y)
          } else if(sig[i, j] < 0.01) {
            grid.text("**", x, y) 
          } else if(sig[i, j] < 0.05) {
            grid.text('*', x, y)
          }}
        )
dev.off()

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
