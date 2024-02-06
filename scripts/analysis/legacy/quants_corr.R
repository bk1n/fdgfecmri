library(pacman)
p_load(tidyverse, infer, Hmisc, ComplexHeatmap)
library(circlize)
library(ggpubr)

fdg =  read_csv('./data/quant_allPooled.csv')

d = fdg %>%
  filter(CANC == 'endo') %>%
  select(-c(CANC, region)) %>%
  select(-contains(c('STAR', 'NTR', 'SNSA', 'SA', 'LA'))) %>%
  select(-HIST)

pairwise_corr = rcorr(as.matrix(d))
sig = pairwise_corr$P
sig[is.na(sig)] = 0
sig = p.adjust(sig, method='fdr') %>% matrix(ncol=length(colnames(pairwise_corr$P)))

ylabs = c('FDG_SUV_PT' = 'FDG SUVmax (PT)',
          'FEC_SUV_PT' = 'FEC SUVmax (PT)',
          'ADC_PT' = 'ADCmean (PT)',
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
          'ADC_NTR' = 'ADCmean NTR')

colnames(pairwise_corr$r) = ylabs[as.character(colnames(pairwise_corr$r))]
rownames(pairwise_corr$r) = ylabs[as.character(rownames(pairwise_corr$r))]

col_fun = colorRamp2(c(min(pairwise_corr$r), 0, max(pairwise_corr$r)), c("blue", "white", "red"))

png('./figures/pairwise_corr_raw.png',
    width = 6,
    height = 5.5,
    units = 'in',
    res = 300)
Heatmap(pairwise_corr$r, col = col_fun,
        heatmap_legend_param = list(title = 'R'),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i != j) {
            if(sig[i, j] < 0.001) {
              text = paste0(signif(pairwise_corr$r[i,j], 2), '***')
              grid.text(text, x, y, gp = gpar(fontsize = 10))
            } else if(sig[i, j] < 0.01) {
              text = paste0(signif(pairwise_corr$r[i,j], 2), '**')
              grid.text(text, x, y, gp = gpar(fontsize = 10)) 
            } else if(sig[i, j] < 0.05) {
              text = paste0(signif(pairwise_corr$r[i,j], 2), '*')
              grid.text(text, x, y, gp = gpar(fontsize = 10)) 
            }} else {
              grid.rect(x, y, w, h, gp = gpar(fill = 'grey', col = NA))
            }}
        )
dev.off()

# g_suv = ggplot(fdg,
#        aes(x = SA,
#            y = SUV)) +
#   geom_point(position = position_jitter(width = 0.3)) +
#   geom_smooth(method = 'lm',
#               formula = y~x,
#               se = T) +
#   stat_cor(method = 'pearson') +
#   ylab('SUVmax') +
#   xlab('Short axis diameter (mm)') +
#   theme_classic() 
#   
# g_adc = ggplot(fdg,
#        aes(x = SA,
#            y = ADC)) +
#   geom_point() +
#   geom_smooth(method = 'lm',
#               formula = y~x,
#               se = T) +
#   stat_cor(method = 'pearson') +
#   xlim(c(0,20)) +
#   ylab('ADCmean') +
#   xlab('Short axis diameter (mm)') +
#   theme_classic() 
# 
# g = ggarrange(g_suv, g_adc,
#           nrow = 1,
#           ncol = 2)
# 
# ggsave('./figures/SUV_ADC_corr_SA.png',
#        g,
#        height = 5,
#        width = 5,
#        units = 'in',
#        dpi = 300)
