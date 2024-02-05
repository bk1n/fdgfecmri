library(tidyverse)
library(cutpointr)

# filter data ----
data = read.csv('./data/processed_data.csv')

data = data %>%
  select(FDG_SUV, FDG_STAR, FDG_NTR,
         FEC_SUV, FEC_STAR, FEC_NTR,
         ADC, ADC_NTR,
         VIS_FDG, VIS_FEC, VIS_MRI,
         HIST, region, CANC) 

# opt cut func ----
get_opt_cut = function(data, measure) {
  # filter df to patients w both visual and quants
  vis_measure = if(grepl('FDG', measure)) 'VIS_FDG' else if(grepl('FEC', measure)) 'VIS_FEC' else 'VIS_MRI'
    
  df = data %>%
    dplyr::filter(!is.na(!!rlang::sym(measure)) & !is.na(!!rlang::sym(vis_measure)))
  
  print(paste('Getting cutpoint for:', measure))
  
  x = df[[measure]]
  class = df$HIST
  if(grepl('ADC', measure)){
    dir = '<='
  } else {
    dir = '>='
  }
  opt_cut = cutpointr(x = x,
                      class = class,
                      na.rm = T,
                      method = maximize_metric,
                      metric = F1_score,
                      pos_class = 1,
                      neg_class = 0, 
                      direction = dir,
                      boot_runs = 20,
                      boot_stratify = T)
  summary.opt_cut = summary(opt_cut)
  
  boot_data = summary.opt_cut$cutpointr[[1]]$boot[[1]]
  boot_data$name = measure

  return(list(opt_cut = opt_cut,
              summary = summary.opt_cut,
              boot_data = boot_data,
              df = df))
}

# run opt cut ----
run_cols = c('FDG_SUV', 'FDG_STAR', 'FDG_NTR',
             'FEC_SUV', 'FEC_STAR', 'FEC_NTR',
             'ADC', 'ADC_NTR')

cutpoint_res = lapply(run_cols, get_opt_cut, data = data)
names(cutpoint_res) = run_cols

plt_df = do.call(rbind, lapply(cutpoint_res, function(r) r$boot_data))

## save opt cut res ----

## plt opt cut res boxplots ----

# refit optimal cutpoints & calculate diagnostic performance ----
## quants
## visual
## mcnemars test
## save results ----

