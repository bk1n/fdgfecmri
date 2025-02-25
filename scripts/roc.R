library(pROC)
library(tidyverse)

source("scripts/tex.R")

data <- read.csv("./data/processed_data.csv")

get_proc <- function(measure, data) {
  dir <- if (measure %in% c("ADC", "ADC_NTR")) "<=" else ">="

  d <- data[c(measure, "HIST")]
  d <- d[rowSums(is.na(d)) == 0, ]

  proc <- cutpointr::roc(d, x = !!measure, class = HIST, pos_class = 1, neg_class = 0, direction = dir)
  proc$name <- measure
  return(proc)
}

measures <- c(
  "FDG_SUV", "FDG_STAR", "FDG_NTR",
  "FEC_SUV", "FEC_STAR", "FEC_NTR",
  "ADC", "ADC_NTR"
)

# endo ----
endo_df <- data %>%
  filter(CANC == "endo")

# run roc
roc_plt_df <- do.call(rbind, lapply(measures, get_proc, data = endo_df))

# mutate roc df
roc_plt <- roc_plt_df %>%
  mutate(measure_type = case_when(
    grepl("FDG", name) ~ "FDG-PET/CT",
    grepl("FEC", name) ~ "FEC-PET/CT",
    grepl("ADC", name) ~ "DW-MRI",
    name == "LR (SUV)" ~ "LR"
  )) %>%
  mutate(
    name = gsub("FDG_", "", name),
    name = gsub("FEC_", "", name),
    name = gsub("ADC_", "", name)
  ) %>%
  mutate(name = case_when(name == "SUV" ~ "SUVmax",
    name == "SA" ~ "SA (mm)",
    name == "LA" ~ "LA (mm)",
    name == "ADC" ~ "ADCmean",
    .default = name
  )) %>%
  mutate(measure_type = factor(measure_type, levels = c("FDG-PET/CT", "FEC-PET/CT", "DW-MRI"))) %>%
  mutate(name = factor(name, levels = c("SUVmax", "STAR", "NTR", "ADCmean")))

# plot roc
g_roc <- ggplot(
  roc_plt,
  aes(
    x = fpr,
    y = 1 - fnr,
    color = name
  )
) +
  geom_line(
    alpha = 0.7,
    linewidth = 1
  ) +
  geom_line(
    inherit.aes = F,
    data = data.frame(x = seq(0, 1), y = seq(0, 1)),
    aes(
      x = x,
      y = y
    ),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  theme_classic() +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  scale_color_brewer(
    name = "Quantitative \nMeasure",
    palette = "Set1",
    direction = -1
  ) +
  scale_linetype_manual(
    name = "Measure Type",
    values = c("solid", "twodash", "dotdash", "dashed")
  ) +
  facet_wrap(~measure_type, labeller = label_bquote(.(c(tex$fdg, tex$fec, tex$mri))))

ggsave("outputs/ROC_endo.png",
  g_roc,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

# cervical ----
cer_df <- data %>%
  filter(CANC == "cer")

# run roc
roc_plt_df <- do.call(rbind, lapply(measures, get_proc, data = cer_df))

# mutate roc df
roc_plt <- roc_plt_df %>%
  mutate(measure_type = case_when(
    grepl("FDG", name) ~ "FDG-PET/CT",
    grepl("FEC", name) ~ "FEC-PET/CT",
    grepl("ADC", name) ~ "DW-MRI",
    name == "LR (SUV)" ~ "LR"
  )) %>%
  mutate(
    name = gsub("FDG_", "", name),
    name = gsub("FEC_", "", name),
    name = gsub("ADC_", "", name)
  ) %>%
  mutate(name = case_when(name == "SUV" ~ "SUVmax",
    name == "SA" ~ "SA (mm)",
    name == "LA" ~ "LA (mm)",
    name == "ADC" ~ "ADCmean",
    .default = name
  )) %>%
  mutate(measure_type = factor(measure_type, levels = c("FDG-PET/CT", "FEC-PET/CT", "DW-MRI"))) %>%
  mutate(name = factor(name, levels = c("SUVmax", "STAR", "NTR", "ADCmean")))

# plot roc
g_roc <- ggplot(
  roc_plt,
  aes(
    x = fpr,
    y = 1 - fnr,
    color = name
  )
) +
  geom_line(
    alpha = 0.7,
    linewidth = 1
  ) +
  geom_line(
    inherit.aes = F,
    data = data.frame(x = seq(0, 1), y = seq(0, 1)),
    aes(
      x = x,
      y = y
    ),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  theme_classic() +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  scale_color_brewer(
    name = "Quantitative \nMeasure",
    palette = "Set1",
    direction = -1
  ) +
  scale_linetype_manual(
    name = "Measure Type",
    values = c("solid", "twodash", "dotdash", "dashed")
  ) +
  facet_wrap(~measure_type, labeller = label_bquote(.(c(tex$fdg, tex$fec, tex$mri))))

ggsave("outputs/ROC_cer.png",
  g_roc,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)
