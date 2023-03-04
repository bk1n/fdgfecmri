library(pacman)
p_load(tidyverse, infer, ggpubr, ggsignif, pROC, cutpointr)

data = readRDS('./data/quant_allPooled.rds')

fdg_features = c('FDG_SUV', 'FDG_SA', 'FDG_LA', 'FDG_NTR', 'FDG_STAR', 'FDG_SNSA')
fec_features = c('FEC_SUV', 'FEC_SA', 'FEC_LA', 'FEC_NTR', 'FEC_STAR', 'FEC_SNSA')
adc_features = c('ADC', 'ADC_NTR')
pt_features = c('FDG_SUV_PT', 'FEC_SUV_PT', 'ADC_PT')
canc = c('cer', 'endo')


plot_myCustomGG = function(feature){
  plt_data = data %>%
    filter(CANC == canc)
  
  ylabs = c('FDG_SUV_PT' = expression(FDG~SUV[max]),
            'FEC_SUV_PT' = expression(FEC~SUV[max]),
            'ADC_PT' = expression(ADC[mean]),
            'FDG_SUV' = expression(FDG~SUV[max]),
            'FDG_SA' = expression(FDG~SA~(mm)),
            'FDG_LA' = expression(FDG~LA~(mm)),
            'FDG_NTR' = expression(FDG~NTR),
            'FDG_STAR' = expression(FDG~STAR),
            'FDG_SNSA' = expression(FDG~SNSA),
            'FEC_SUV' = expression(FEC~SUV[max]),
            'FEC_SA' = expression(FEC~SA~(mm)),
            'FEC_LA' = expression(FEC~LA~(mm)),
            'FEC_NTR' = expression(FEC~NTR),
            'FEC_STAR' = expression(FEC~STAR),
            'FEC_SNSA' = expression(FEC~SNSA),
            'ADC' = expression(ADC[mean]),
            'ADC_NTR' = expression(ADC[mean]~NTR))
  
  print(feature)

  p = signif(wilcox.test(pull(plt_data[,feature]) ~ plt_data$HIST)$p.value,2)
  if(p < 0.001){
    p = paste0('p =', p, '***')
  } else if(p >= 0.001 & p < 0.01){
    p = paste0('p = ', p, '**')
  } else if(p >= 0.01 & p < 0.05){
    p = paste0('p = ', p, '*')
  } else{
    p = paste('p = ', p)
  }
  
  g = ggplot(plt_data,
         aes(x = HIST,
             y = !!rlang::sym(feature))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.1)) +
    geom_signif(comparisons = list(c('0','1')),
                annotations = c(p)) +
    ylab(ylabs[feature]) +
    xlab('') +
    scale_x_discrete(label = c('B', 'M')) +
    theme_classic()

  return(g)
}

canc = 'cer'
gg_cer_pt = lapply(pt_features, plot_myCustomGG)
gg_cer_fdg = lapply(fdg_features, plot_myCustomGG)
gg_cer_fec = lapply(fec_features, plot_myCustomGG)
gg_cer_adc = lapply(adc_features, plot_myCustomGG) 

canc = 'endo'
gg_endo_pt = lapply(pt_features, plot_myCustomGG)
gg_endo_fdg = lapply(fdg_features, plot_myCustomGG)
gg_endo_fec = lapply(fec_features, plot_myCustomGG)
gg_endo_adc = lapply(adc_features, plot_myCustomGG) 

#PT in ENDO, CER
g = ggarrange(plotlist = c(gg_endo_pt, gg_cer_pt), 
              ncol = length(pt_features), 
              nrow = 2,
              labels = c('A', rep('', length(pt_features) - 1), 'B')) 
ggsave('./figures/quants_pt_byCancer.png',
       g,
       width = 6,
       height = 8,
       units = 'in',
       dpi = 300,
       bg = 'white')

#FDG, FEC in ENDO
g = ggarrange(plotlist = c(gg_endo_fdg, gg_endo_fec), 
              ncol = length(fdg_features), 
              nrow = 2,
              labels = c('A', rep('', length(fdg_features) - 1), 'B')) 
ggsave('./figures/fdg_fec_inEndo.png',
       g,
       width = 12,
       height =8,
       units = 'in',
       dpi = 300,
       bg = 'white')

# FDG, FEC in CER
g = ggarrange(plotlist = c(gg_cer_fdg, gg_cer_fec), 
              ncol = length(fdg_features), 
              nrow = 2,
              labels = c('A', rep('', length(fdg_features) - 1), 'B')) 
ggsave('./figures/fdg_fec_inCer.png',
       g,
       width = 12,
       height = 8,
       units = 'in',
       dpi = 300,
       bg = 'white')

#ADC in ENDO,CER
g = ggarrange(plotlist = c(gg_endo_adc, gg_cer_adc), 
              ncol = length(adc_features), 
              nrow = 2,
              labels = c('A', rep('', length(adc_features) - 1), 'B')) 
ggsave('./figures/ADC_byCancer.png',
       g,
       width = 4,
       height = 8,
       units = 'in',
       dpi = 300,
       bg = 'white')

