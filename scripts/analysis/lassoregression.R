library(pacman)
p_load(tidyverse, caret, glmnet, parallel)

fdg = readRDS('./data/fdgpet_allPooled.rds')

# filter data so there are no missing SUV values
fdg_lr = fdg %>%
  filter(CANC == 'endo') %>%
  dplyr::select(-c(CANC, region, PT_SUV)) %>%
  dplyr::select(-c(LA, ADC,SNSA)) %>%
  filter(!is.na(HIST),
         !is.na(SUV))

runLogitLASSO = function(x){
  # define training and test set indices with 70% 30% split
  train_ind = createDataPartition(fdg_lr$HIST, p = 0.7, times = 1, list = F)
  
  # build training and test sets
  lr_train = fdg_lr[train_ind,]
  lr_test = fdg_lr[-train_ind,]
  
  train_hist = lr_train$HIST
  test_hist = lr_test$HIST
  lr.train.mat = lr_train %>%
    dplyr::select(-HIST) %>%
    as.matrix()
  lr.test.mat = lr_test %>%
    dplyr::select(-HIST) %>%
    as.matrix()
  
  # cross validate model to find lambda that minimises MSE
  cv.glm_opt = cv.glmnet(x = lr.train.mat, y = as.numeric(train_hist), alpha = 1, family = 'binomial')
  best_lambda = cv.glm_opt$lambda.min
  
  # build final model with this data!
  glm.model = glmnet(x = lr.train.mat, y = train_hist, alpha = 1, lambda = best_lambda, family = 'binomial')

  # calculate predictions for test data and define actual obs
  preds = predict(glm.model, newx = lr.test.mat, type="response", s="lambda.min")
  obs = as.integer(as.character(test_hist))
  
  #func to calc f1
  calc_f1 = function(preds, obs, cutoff) {
    preds_class = ifelse(preds > cutoff, 1, 0)
    
    # calculate tp, fp, fn, tn matrix
    tp = sum(preds_class == 1 & obs == 1)
    fp = sum(preds_class == 1 & obs == 0)
    fn = sum(preds_class == 0 & obs == 1)
    tn = sum(preds_class == 0 & obs == 0)
    
    #f1 score
    f1 = 2*tp/(2*tp + fp + fn)
    
    return(data.frame(cutoff, f1))
  }
  
  # calculate optimal cutoff based on f1
  f1_byCutoff = do.call(rbind,lapply(seq(1/1000, 1, length.out = 1000), calc_f1, preds=preds, obs=obs))
  optimal_cutoff = median(f1_byCutoff$cutoff[which(f1_byCutoff$f1 == max(f1_byCutoff$f1))])
  
  # categorise predictions based on > or < 0.5 
  preds_class = ifelse(preds > optimal_cutoff, 1, 0)
  # calculate accuracy
  acc = mean(preds_class == obs)
  # calculate rmse
  rmse = sqrt(mean(((preds - obs)**2)))
  # calculate tp, fp, fn, tn matrix
  tp = sum(preds_class == 1 & obs == 1)
  fp = sum(preds_class == 1 & obs == 0)
  fn = sum(preds_class == 0 & obs == 1)
  tn = sum(preds_class == 0 & obs == 0)
  
  #confusion matrix
  confusion_mat = data.frame(c(tp, fn), c(fp, tn))
  colnames(confusion_mat) = c('malignant', 'benign')
  rownames(confusion_mat) = c('positive', 'negative')
  
  #sensitivty, specifity
  sens = tp / (tp + fn)
  spec = tn / (tn + fn)
  
  #fpr, fnr
  fnr = 1 - sens
  fpr = 1 - spec
  
  #plr, nlr
  plr = sens / fpr
  nlr = fnr / spec
  
  #f1 score
  f1 = 2*tp/(2*tp + fp + fn)
  
  #model lambda
  best_lambda = glm.model$lambda
  
  co = coef(glm.model)
  n_suv = 0
  n_sa = 0
  n_ntr = 0
  n_star = 0
  if(co['SUV',] > 0) n_suv = 1
  else if(co['SA',] > 0) n_sa = 1
  else if(co['NTR',] > 0) n_ntr = 1
  else if(co['STAR'] > 0) n_star = 1
  
  res.model = data.frame(optimal_cutoff, acc, rmse, f1, sens, spec, fnr, fpr, plr, nlr, best_lambda, n_suv, n_sa, n_ntr, n_star)
  
  return(res.model)
}

set.seed = 123
res.model = as.data.frame(do.call(rbind, mclapply(c(rep(NA, 5000)), runLogitLASSO))) # use do.call to rbind our resulting list of data.frames
res.model.avg = t(as.matrix(signif(colMeans(select(res.model, optimal_cutoff:best_lambda), na.rm = T),3)))
res.model.avg.param_counts = t(as.matrix(colSums(select(res.model, n_suv:n_star))))
res.model.avg = as.data.frame(cbind(res.model.avg, res.model.avg.param_counts))
print('----------- LASSO Regression Model Performance after 5000 resamples')
print(res.model.avg)
print('-----------')


write.csv(res.model.avg, './figures/tables/lasso.model_performance.csv')

# see - imbalanced dataset - consider the F1 score - and calculate specificty and sensitivity!
# https://stats.stackexchange.com/questions/449220/probability-cut-off-value-for-logistic-regression