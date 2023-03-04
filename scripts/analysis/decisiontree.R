library(pacman)
p_load(tidyverse, caret, rpart, partykit)

fdg = readRDS('./data/fdgpet_allPooled.rds')

fdg_dt = fdg %>%
  filter(CANC == 'endo') %>%
  select(-c(CANC, region)) %>%
  select(SUV, NTR, ADC, HIST) %>%
  filter(!is.na(SUV) & !is.na(NTR) & !is.na(ADC))

colSums(!is.na(fdg_dt))

set.seed(123)
train_ind = createDataPartition(fdg_dt$HIST, p = 0.8, times = 1, list = F)

# build training and test sets
dt_train = fdg_dt[train_ind,]
dt_test = fdg_dt[-train_ind,]

fit = ctree(formula = HIST ~ .,
            data = dt_train,
            control = ctree_control(minsplit = 0,
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

print(
  data.frame(acc = acc,
           sens = sens,
           spec = spec,
           f1 = f1)
)

plot(fit)
