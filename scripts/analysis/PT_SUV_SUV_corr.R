
library(pacman)
p_load(tidyverse, infer)

fdg = readRDS('./data/fdgpet_allPooled.rds')

ggplot(fdg,
       aes(x = PT_SUV,
           y = SUV)) +
  geom_point()