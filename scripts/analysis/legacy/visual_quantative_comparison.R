library(tidyverse)

vis = read_csv('./data/visual_data.csv')
quants = read_csv('./data/quant_allPooled.csv')

q = quants %>%
  filter(CANC == 'endo') %>%
  select(ADC, HIST) %>%
  filter(!is.na(ADC)) %>%
  mutate(ADC_ = if_else(ADC >= cuts$optimal_cutpoint[cuts$name == 'ADC'], 1, 0)) %>%
  mutate(across(everything(), function(x) factor(x, levels = c(1,0))))

df = quants %>%
  left_join(vis, by = c('PATIENT_ID', 'region')) %>%
  mutate(across(contains('VIS'), ~ if_else(.x >= 5, 1, 0)))

# mutate dataframe to x measure + visual corresponding 
d = df %>%
  filter(CANC == 'endo') %>%
  mutate(FDG_SUV_ = if_else(FDG_SUV >= cuts$optimal_cutpoint[cuts$name == 'FDG_SUV'], 1, 0),
         FEC_SUV_ = if_else(FEC_SUV >= cuts$optimal_cutpoint[cuts$name == 'FEC_SUV'], 1, 0),
         ADC_ = if_else(ADC <= cuts$optimal_cutpoint[cuts$name == 'ADC'], 1, 0)) %>%
  select(FDG_SUV_, FEC_SUV_, ADC_, VIS_FDG, VIS_FEC, VIS_MRI, HIST) %>%
  mutate(across(everything(), function(x) factor(x, levels = c(1,0))))

# select q measure, visual diagnosis, actual hist
fdg = d %>%
  select(FDG_SUV_, VIS_FDG, HIST) %>%
  filter(!is.na(FDG_SUV_) & !is.na(VIS_FDG))

fec = d %>%
  select(FEC_SUV_, VIS_FEC, HIST) %>%
  filter(!is.na(FEC_SUV_) & !is.na(VIS_FEC))

adc = d %>%
  select(ADC_, VIS_MRI, HIST) %>%
  filter(!is.na(ADC_) & !is.na(VIS_MRI))

# calculate sensitivity + specificity of both
sens_spec = function(cm) {
  tp = cm[1,1]
  fn = cm[2,1]
  fp = cm[1,2]
  tn = cm[2,2]
  
  sens = (tp / (tp + fn))
  spec = (tn / (fp + tn))
  return(c(sens, spec))
}

# fdg
sens_spec(table(fdg$FDG_SUV_, fdg$HIST))
sens_spec(table(fdg$VIS_FDG, fdg$HIST))

# fec
sens_spec(table(fec$FEC_SUV_, fec$HIST))
sens_spec(table(fec$VIS_FEC, fec$HIST))

# mri
sens_spec(table(adc$ADC_, adc$HIST))
sens_spec(table(adc$VIS_MRI, adc$HIST))

# McNemar's test
# need to find + import cut-offs for each measure
cuts = read_csv('./figures/optimalCP_byMeasure.csv')

# for each measure (FDG_SUV, FEC_SUV, ADC) do mcnemars test for positive and negative 
# cm returns response columns (visual measure)

convert_cm_to_df = function(measure_name, cm, pval){
  return(
    data.frame('v_correct_q_correct' = cm[1,1],
               'v_correct_q_incorrect' = cm[2,1],
               'v_incorrect_q_correct' = cm[1,2],
               'v_incorrect_q_incorrect' = cm[2,2],
               pval = pval,
               row.names = measure_name)
  )
}

mcnemar_test = function(measure_name, visual_measure_name, cutpoint, df) {
  d = df %>%
    select(HIST, all_of(c(measure_name, visual_measure_name))) %>%
    mutate(cut = if_else(!!rlang::sym(measure_name) > cutpoint, 1, 0))
  
  hist_pos = d %>%
    filter(HIST == 1)
  hist_neg = d %>%
    filter(HIST == 0)
  
  pos_cm = table(
    factor(pull(hist_pos, 'cut') == hist_pos$HIST, levels = c(T,F)),
    factor(pull(hist_pos, visual_measure_name) == hist_pos$HIST, levels = c(T,F))
  )
  pos_mcnemar = mcnemar.test(pos_cm)
  pos_df = convert_cm_to_df(paste0(measure_name, '_pos'),
                            pos_cm, 
                            pos_mcnemar$p.value)
  
  neg_cm = table(
    factor(pull(hist_neg, 'cut') == hist_neg$HIST, levels = c(T,F)),
    factor(pull(hist_neg, visual_measure_name) == hist_neg$HIST, levels = c(T,F))
  )
  neg_mcnemar = mcnemar.test(neg_cm)
  neg_df = convert_cm_to_df(paste0(measure_name, '_neg'),
                            neg_cm, 
                            neg_mcnemar$p.value)
  
  res_df = rbind(pos_df, neg_df)

  return(
    list(pos_cm = pos_cm, 
         pos_mcnemar = pos_mcnemar, 
         neg_cm = neg_cm, 
         neg_mcnemar = neg_mcnemar,
         df = res_df)
  )
}

## endo ####
endo = df %>%
  filter(CANC == 'endo')

#fdg
fdg = mcnemar_test(measure_name = 'FDG_SUV', 
             visual_measure_name = 'VIS_FDG', 
             cutpoint = cuts$optimal_cutpoint[cuts$name == 'FDG_SUV'],
             df = endo)
fdg_star = mcnemar_test(measure_name = 'FDG_STAR', 
                        visual_measure_name = 'VIS_FDG', 
                        cutpoint = cuts$optimal_cutpoint[cuts$name == 'FDG_STAR'],
                        df = endo)
fdg_ntr = mcnemar_test(measure_name = 'FDG_NTR', 
                       visual_measure_name = 'VIS_FDG', 
                       cutpoint = cuts$optimal_cutpoint[cuts$name == 'FDG_NTR'],
                       df = endo)

#fec
fec = mcnemar_test(measure_name = 'FEC_SUV', 
                   visual_measure_name = 'VIS_FEC', 
                   cutpoint = cuts$optimal_cutpoint[cuts$name == 'FEC_SUV'],
                   df = endo)
fec_star = mcnemar_test(measure_name = 'FEC_STAR', 
                        visual_measure_name = 'VIS_FEC', 
                        cutpoint = cuts$optimal_cutpoint[cuts$name == 'FEC_STAR'],
                        df = endo)
fec_ntr = mcnemar_test(measure_name = 'FEC_NTR', 
                       visual_measure_name = 'VIS_FEC', 
                       cutpoint = cuts$optimal_cutpoint[cuts$name == 'FEC_NTR'],
                       df = endo)

#adc
adc = mcnemar_test(measure_name = 'ADC', 
             visual_measure_name = 'VIS_MRI', 
             cutpoint = cuts$optimal_cutpoint[cuts$name == 'ADC'],
             df = endo)
adc_ntr = mcnemar_test(measure_name = 'ADC_NTR', 
             visual_measure_name = 'VIS_MRI', 
             cutpoint = cuts$optimal_cutpoint[cuts$name == 'ADC_NTR'],
             df = endo)

endo_results_df = rbind(
  fdg$df,
  fdg_ntr$df,
  fdg_star$df,
  fec$df,
  fec_ntr$df,
  fec_star$df,
  adc$df,
  adc_ntr$df
)

write.csv(endo_results_df, './figures/endo_mcnemar_test_results.csv')

# cer #### 
cer = df %>%
  filter(CANC == 'cer')

#fdg
fdg = mcnemar_test(measure_name = 'FDG_SUV', 
                   visual_measure_name = 'VIS_FDG', 
                   cutpoint = cuts$optimal_cutpoint[cuts$name == 'FDG_SUV'],
                   df = cer)
fdg_star = mcnemar_test(measure_name = 'FDG_STAR', 
                        visual_measure_name = 'VIS_FDG', 
                        cutpoint = cuts$optimal_cutpoint[cuts$name == 'FDG_STAR'],
                        df = cer)
fdg_ntr = mcnemar_test(measure_name = 'FDG_NTR', 
                       visual_measure_name = 'VIS_FDG', 
                       cutpoint = cuts$optimal_cutpoint[cuts$name == 'FDG_NTR'],
                       df = cer)

#fec
fec = mcnemar_test(measure_name = 'FEC_SUV', 
                   visual_measure_name = 'VIS_FEC', 
                   cutpoint = cuts$optimal_cutpoint[cuts$name == 'FEC_SUV'],
                   df = cer)
fec_star = mcnemar_test(measure_name = 'FEC_STAR', 
                        visual_measure_name = 'VIS_FEC', 
                        cutpoint = cuts$optimal_cutpoint[cuts$name == 'FEC_STAR'],
                        df = cer)
fec_ntr = mcnemar_test(measure_name = 'FEC_NTR', 
                       visual_measure_name = 'VIS_FEC', 
                       cutpoint = cuts$optimal_cutpoint[cuts$name == 'FEC_NTR'],
                       df = cer)

#adc
adc = mcnemar_test(measure_name = 'ADC', 
                   visual_measure_name = 'VIS_MRI', 
                   cutpoint = cuts$optimal_cutpoint[cuts$name == 'ADC'],
                   df = cer)
adc_ntr = mcnemar_test(measure_name = 'ADC_NTR', 
                       visual_measure_name = 'VIS_MRI', 
                       cutpoint = cuts$optimal_cutpoint[cuts$name == 'ADC_NTR'],
                       df = cer)

cer_results_df = rbind(
  fdg$df,
  fdg_ntr$df,
  fdg_star$df,
  fec$df,
  fec_ntr$df,
  fec_star$df,
  adc$df,
  adc_ntr$df
)

write.csv(cer_results_df, './figures/cer_mcnemar_test_results.csv')
