library(pacman)
p_load(tidyverse, caret, parallel)

fdg = readRDS('./data/fdgpet_allPooled.rds')

# filter data so there are no missing SUV values
fdg_lr = fdg %>%
  filter(CANC == 'endo') %>%
  dplyr::select(-c(CANC, region, PT_SUV)) %>%
  dplyr::select(-c(LA, ADC,SNSA)) %>%
  filter(!is.na(HIST),
         !is.na(SUV)) %>%
  mutate(HIST = as.numeric(as.character(HIST)))

runLogitRegr = function(x){
  # define training and test set indices with 70% 30% split
  train_ind = createDataPartition(fdg_lr$HIST, p = 0.7, times = 1, list = F)
  
  # build training and test sets
  lr_train = fdg_lr[train_ind,]
  lr_test = fdg_lr[-train_ind,]
  
  # cross validate model to find lambda that minimises MSE
  lr.model = glm(HIST ~ SUV, 
                    lr_train,
                    family = binomial())
  
  # calculate predictions for test data and define actual obs
  preds = as.numeric(predict(lr.model, newdata = lr_test, type="response"))
  obs = as.integer(as.character(lr_test$HIST))
  
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
    
    spec = tn / (tn + fp)
    
    return(data.frame(cutoff, spec))
  }
  
  # calculate optimal cutoff based on f1
  
  f1_byCutoff = na.omit(do.call(rbind,lapply(seq(1/1000, 1, length.out = 1000), calc_f1, preds=preds, obs=obs)))
  optimal_cutoff = median(f1_byCutoff$cutoff[which(f1_byCutoff$spec == max(f1_byCutoff$spec))])
  
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
  spec = tn / (tn + fp)
  
  #fpr, fnr
  fnr = 1 - sens
  fpr = 1 - spec
  
  #plr, nlr
  plr = sens / fpr
  nlr = fnr / spec
  
  #f1 score
  f1 = 2*tp/(2*tp + fp + fn)
  
  #SUV cutoff
  inverse.logit = function(model, p) {(log(p/(1-p)) - as.numeric(model$coefficients[1]))/as.numeric(model$coefficients[2])}
  optimal_SUVcutoff = inverse.logit(lr.model, optimal_cutoff)
  
  res.model = data.frame(optimal_cutoff, optimal_SUVcutoff, acc, rmse, f1, sens, spec, fnr, fpr, plr, nlr)
  
  return(res.model)
}

set.seed = 123
res.model = as.data.frame(do.call(rbind, lapply(c(rep(NA, 50)), runLogitRegr))) # use do.call to rbind our resulting list of data.frames
res.model.avg = signif(colMeans(res.model, na.rm = T),3)
res.model.avg = data.frame(t(res.model.avg))
print('----------- Linear Regression Model Performance after 5000 resamples (F1 score)')
print(res.model.avg)
print('-----------')

write.csv(res.model.avg, './figures/tables/lr.model_performance.csv')

lr.model.full = glm(HIST ~ SUV, 
                    fdg_lr,
                    family = binomial())

## plot logistic regression
SUV = seq(0, to = 30, by = 0.1)
x = data.frame(SUV)
y = as.numeric(predict(lr.model.full, newdata = x, type = 'response'))

annot = data.frame(x = rep(15,5),
                   y = seq(0.6, 0.4, length.out = 5),
                   label = c(
                     paste0('Model Performance:'),
                     paste0('Optimal SUV cutoff: ', res.model.avg$optimal_SUVcutoff),
                     paste0('F1 score: ', res.model.avg$f1),
                     paste0('Sensitivity: ', res.model.avg$sens*100, '%'),
                     paste0('Specificity: ', res.model.avg$spec*100, '%'))
)

g_lr = ggplot(fdg_lr, 
       aes(x = SUV,
           y = HIST)) +
  stat_smooth(method="glm", 
              se=F,
              formula = y~x,
              method.args = list(family=binomial)) +
  geom_point(position = position_jitter(height = 0.01), alpha = 0.5) + 
  geom_vline(xintercept = res.model.avg$optimal_SUVcutoff,
             color = 'red',
             linetype = 'dashed') +
  xlab('SUVmax') +
  ylab(expression(Pr~(Malignant))) +
  geom_text(data = annot,
            aes(x = x,
            y = y,
            label=label), hjust = 0) +
  theme_classic()
g_lr

ggsave('./figures/lr.model_plt.png',
       g_lr,
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)
