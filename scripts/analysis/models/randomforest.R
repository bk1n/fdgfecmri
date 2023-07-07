library(pacman)
p_load(tidyverse, caret, randomForest)

fdg = readRDS('./data/fdgpet_allPooled.rds')

# class imblanace issue w rf - https://stats.stackexchange.com/questions/242833/is-random-forest-a-good-option-for-unbalanced-data-classification

fdg_rf  = fdg %>%
  filter(CANC == 'endo') %>%
  select(-c(CANC, region)) 

set.seed(123)
train_ind = createDataPartition(fdg_rf$HIST, p = 0.7, times = 1, list = F)

# build training and test sets
rf_train = fdg_rf[train_ind,]
rf_test = fdg_rf[-train_ind,]

rf = randomForest(HIST ~ SUV + NTR + SA + LA,
             data = rf_train,
             proximity = T,
             na.action = na.omit)

p1 = predict(rf, rf_test)
confusionMatrix(p1, rf_test$HIST)
