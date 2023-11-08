library(pacman)
p_load(tidyverse, caret, randomForest)

fdg = read_csv('./data/quant_allPooled.csv')

# class imblanace issue w rf - https://stats.stackexchange.com/questions/242833/is-random-forest-a-good-option-for-unbalanced-data-classification

fdg_rf  = fdg %>%
  filter(CANC == 'endo') %>%
  select(-c(CANC, region)) %>%
  select(HIST, ADC_NTR, FDG_NTR) %>%
  filter(!is.na(ADC_NTR) & !is.na(FDG_NTR)) %>%
  mutate(HIST = factor(HIST))

fit_rf = function(...) {  
  train_ind = createDataPartition(fdg_rf$HIST, p = 0.7, times = 1, list = F)
  
  # build training and test sets
  rf_train = fdg_rf[train_ind,]
  rf_test = fdg_rf[-train_ind,]
  
  rf = randomForest(HIST ~ ADC_NTR + FDG_NTR,
               data = rf_train,
               proximity = T,
               na.action = na.omit)
  
  p1 = predict(rf, rf_test)
  obs = rf_test$HIST
  
  tp = sum(p1 == 1 & obs == 1)
  fp = sum(p1 == 1 & obs == 0)
  fn = sum(p1 == 0 & obs == 1)
  tn = sum(p1 == 0 & obs == 0)
  
  acc = mean(p1 == obs)
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
    data.frame(acc = acc,
               sens = sens,
               spec = spec,
               f1 = f1)
    )
}

set.seed(123)
res = lapply(1:2000, fit_rf)

r = do.call(rbind, res)

r %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = c('acc', 'sens', 'spec', 'f1')) %>%
  ggplot(aes(x = name,
             y = value)) +
  geom_boxplot()
