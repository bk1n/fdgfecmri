library(pacman)
p_load(tidyverse, infer, ggpubr, ggsignif, pROC, cutpointr)

fdg = readRDS('./data/fdgpet_allPooled.rds')

fdg = fdg %>%
  mutate(SA_scaled = ((SA - min(SA, na.rm = T))/(max(SA, na.rm = T) - min(SA, na.rm = T)))+1) %>%
  mutate(SA_log = log(SA_scaled)) %>%
  mutate(S = SUV * SA_scaled)

fdg %>%
  filter(CANC == 'endo') %>%
  arrange(desc(SUV))

features = c('SUV', 'ADC', 'NTR', 'STAR', 'SNSA', 'SA', 'LA')
canc = c('cer', 'endo')

plot_myCustomGG = function(feature){
  plt_data = fdg %>%
    filter(CANC == canc)
  
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
    xlab('') +
    scale_x_discrete(label = c('B', 'M')) +
    theme_classic()

  return(g)
}

canc = 'cer'
gg_cer = lapply(features, plot_myCustomGG)

canc = 'endo'
gg_endo = lapply(features, plot_myCustomGG)

g = ggarrange(plotlist = c(gg_cer, gg_endo), 
              ncol = length(features), 
              nrow = 2,
              labels = c('A', rep('', length(features) - 1), 'B')) 
g = annotate_figure(g,
                bottom = text_grob("Lymph node histology")) 
g

ggsave('./figures/quants_byCancer.png',
       g,
       width = 15,
       height =8,
       units = 'in',
       dpi = 300,
       bg = 'white')
