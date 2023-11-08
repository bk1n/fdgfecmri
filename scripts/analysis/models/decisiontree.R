library(pacman)
p_load(tidyverse, caret, rpart, partykit)

fdg = read_csv('./data/quant_allPooled.csv')

sum(rowSums(!is.na(select(fdg, FDG_NTR, ADC_NTR))) > 1)
sum(rowSums(!is.na(select(fdg, FDG_SUV, ADC))) > 1)
# 71 FDG NTR + ADC NTR
# sum(rowSums(!is.na(select(fdg, FDG_SUV, ADC_NTR))) > 1)

fdg_dt = fdg %>%
  filter(CANC == 'endo') %>%
  select(-c(CANC, region)) %>%
  select(HIST, FDG_SUV, FDG_NTR, ADC, ADC_NTR) %>%
  mutate(r = rowSums(is.na(.))) %>%
  filter(r == 0) %>%
  select(-r) %>%
  mutate(HIST = factor(HIST))

colSums(!is.na(fdg_dt))

fit_tree = function(maxdepth){
  train_ind = createDataPartition(fdg_dt$HIST, p = 0.7, times = 1, list = F)
  
  # build training and test sets
  dt_train = fdg_dt[train_ind,]
  dt_test = fdg_dt[-train_ind,]
  
  fit = ctree(formula = HIST ~ .,
              data = dt_train,
              control = ctree_control(maxdepth = maxdepth,
                                      logmincriterion = log(0),
                                      intersplit = T))
  
  preds_class = predict(fit, dt_test, type = 'response')
  obs = dt_test$HIST

  tp = sum(preds_class == 1 & obs == 1)
  fp = sum(preds_class == 1 & obs == 0)
  fn = sum(preds_class == 0 & obs == 1)
  tn = sum(preds_class == 0 & obs == 0)
  
  acc = mean(preds_class == obs)
  #sensitivty, specifity
  sens = tp / (tp + fn)
  spec = tn / (tn + fp)
  
  #fpr, fnr
  fnr = 1 - sens
  fpr = 1 - spec
  
  #plr, nlr
  plr = sens / fpr
  nlr = fnr / spec
  
  #f1 score
  f1 = 2*tp/(2*tp + fp + fn)
  
  return(
    data.frame(maxdepth = maxdepth,
               acc = acc,
             sens = sens,
             spec = spec,
             f1 = f1)
  )
}

set.seed(123)
res = lapply(sample(1:8, 2000, replace =T), fit_tree)

r = do.call(rbind, res)

r %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = c('acc', 'sens', 'spec', 'f1')) %>%
  mutate(maxdepth = factor(maxdepth)) %>%
  ggplot(aes(x = name,
             y = value,
             fill = maxdepth)) +
  geom_boxplot()

plot(fit)
