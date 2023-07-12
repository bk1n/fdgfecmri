library(pacman)
p_load(tidyverse, infer, ggpubr, ggsignif, pROC, cutpointr)

data = readRDS('./data/quant_allPooled.rds')

quant_row_1 = c('FDG_SUV', 'FDG_SA', 'FDG_LA', 'FDG_NTR', 'FDG_STAR', 'FDG_SNSA', 'ADC')
quant_row_2 = c('FEC_SUV', 'FEC_SA', 'FEC_LA', 'FEC_NTR', 'FEC_STAR', 'FEC_SNSA', 'ADC_NTR')

opt_cut = read_csv('./figures/tables/optimalCP_byMeasure.csv')

pt_features = c('FDG_SUV_PT', 'FEC_SUV_PT', 'ADC_PT')
canc = c('cer', 'endo')

plot_myCustomGG = function(feature, show_opt_cut = F){
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
  
  ylabs = c('FDG_SUV_PT' = 'FDG SUVmax',
            'FEC_SUV_PT' = 'FEC SUVmax',
            'ADC_PT' = 'ADCmean',
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
    geom_point(position = position_jitter(width = 0.1), alpha = .7) +
    geom_signif(comparisons = list(c('0','1')),
                annotations = c(p)) +
    ylab(ylabs[feature]) +
    xlab('') +
    scale_x_discrete(label = c('B', 'M')) +
    theme_classic()
  
  if(show_opt_cut){
    oc = opt_cut %>%
      filter(name == feature) %>%
      pull(optimal_cutpoint)
    g = g + 
      geom_hline(yintercept = oc,
                 color = 'red',
                 linetype = 'dashed',
                 linewidth = 1)
  }

  return(g)
}

canc = 'cer'
gg_cer_pt = lapply(pt_features, plot_myCustomGG)
gg_cer_r1 = lapply(quant_row_1, plot_myCustomGG)
gg_cer_r2 = lapply(quant_row_2, plot_myCustomGG)
# gg_cer_adc = lapply(adc_features, plot_myCustomGG) 

canc = 'endo'
gg_endo_pt = lapply(pt_features, plot_myCustomGG)
gg_endo_r1 = lapply(quant_row_1, plot_myCustomGG, T)
gg_endo_r2 = lapply(quant_row_2, plot_myCustomGG, T)
# gg_endo_adc = lapply(adc_features, plot_myCustomGG)

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

#FDG, FEC, ADC in ENDO
g = ggarrange(plotlist = c(gg_endo_r1, gg_endo_r2), 
              ncol = length(gg_endo_r1), 
              nrow = 2,
              labels = c('A', rep('', length(gg_endo_r1) - 2), 'C', 'B')) 
ggsave('./figures/fdg_fec_inEndo.png',
       g,
       width = 13,
       height =8,
       units = 'in',
       dpi = 300,
       bg = 'white')

# FDG, FEC, ADC in CER
g = ggarrange(plotlist = c(gg_cer_r1, gg_cer_r2),
              ncol = length(gg_cer_r1),
              nrow = 2,
              labels = c('A', rep('', length(gg_cer_r1) - 2), 'C', 'B'))
ggsave('./figures/fdg_fec_inCer.png',
       g,
       width = 13,
       height = 8,
       units = 'in',
       dpi = 300,
       bg = 'white')

#ADC in ENDO,CER
# g = ggarrange(plotlist = c(gg_endo_adc, gg_cer_adc), 
#               ncol = length(adc_features), 
#               nrow = 2,
#               labels = c('A', rep('', length(adc_features) - 1), 'B')) 
# ggsave('./figures/ADC_byCancer.png',
#        g,
#        width = 4,
#        height = 8,
#        units = 'in',
#        dpi = 300,
#        bg = 'white')

