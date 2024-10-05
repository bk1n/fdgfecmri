library(tidyverse)
library(ggsignif)
library(ggpubr)

data <- read.csv("./data/processed_data.csv")

# cervical vs endometrial ----
g_fdg <- ggplot(
  data,
  aes(
    x = CANC,
    y = FDG_SUV
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1),
    alpha = 0.5
  ) +
  geom_signif(
    comparisons = list(c("cer", "endo")),
    test = wilcox.test,
    map_signif_level = function(p) paste0("p = ", signif(p, 3))
  ) +
  scale_x_discrete(label = c("Cervical", "Endometrial")) +
  ylab("FDG SUVmax") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

g_fec <- ggplot(
  data,
  aes(
    x = CANC,
    y = FEC_SUV
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1),
    alpha = 0.5
  ) +
  geom_signif(
    comparisons = list(c("cer", "endo")),
    annotations = c(paste0("p = ", signif(wilcox.test(FEC_SUV ~ CANC, data = data, na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)$p.value, 3), ""))
  ) +
  scale_x_discrete(label = c("Cervical", "Endometrial")) +
  ylab("FEC SUVmax") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

g_adc <- ggplot(
  data,
  aes(
    x = CANC,
    y = ADC
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1),
    alpha = 0.5
  ) +
  geom_signif(
    comparisons = list(c("cer", "endo")),
    annotations = c(paste("p =", signif(wilcox.test(ADC ~ CANC, data = data, na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)$p.value, 3)))
  ) +
  scale_x_discrete(label = c("Cervical", "Endometrial")) +
  ylab("ADCmean") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

g <- ggarrange(g_fdg, g_fec, g_adc, nrow = 1)

ggsave("outputs/boxplot_raw_quants_byCancer.png",
  g,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

# region ----
# fdg
res.aov <- aov(FDG_SUV ~ region, data)
summary(res.aov)
TukeyHSD(res.aov)

g_fdg <- ggplot(
  data,
  aes(
    x = region,
    y = FDG_SUV
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1),
    alpha = 0.5
  ) +
  scale_x_discrete(label = c("Left Pelvis", "Para-aortic", "Right Pelvis")) +
  geom_signif(
    comparisons = list(
      c("paln", "lp"),
      c("rp", "lp"),
      c("rp", "paln")
    ),
    annotations = c(paste("p =", signif(TukeyHSD(res.aov)$region[1:3, "p adj"], digits = 2))),
    y_position = c(30, 38, 34)
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  ylab("FDG SUVmax") +
  xlab("")

# fec
res.aov <- aov(FEC_SUV ~ region, data)
summary(res.aov)
TukeyHSD(res.aov)

g_fec <- ggplot(
  data,
  aes(
    x = region,
    y = FEC_SUV
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1),
    alpha = 0.5
  ) +
  scale_x_discrete(label = c("Left Pelvis", "Para-aortic", "Right Pelvis")) +
  geom_signif(
    comparisons = list(
      c("paln", "lp"),
      c("rp", "lp"),
      c("rp", "paln")
    ),
    annotations = c(paste("p =", signif(TukeyHSD(res.aov)$region[1:3, "p adj"], digits = 2))),
    y_position = c(9, 11, 10)
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  ylab("FEC SUVmax") +
  xlab("")

# adc
res.aov <- aov(ADC ~ region, data)
summary(res.aov)
TukeyHSD(res.aov)

g_adc <- ggplot(
  data,
  aes(
    x = region,
    y = ADC
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1),
    alpha = 0.5
  ) +
  scale_x_discrete(label = c("Left Pelvis", "Para-aortic", "Right Pelvis")) +
  geom_signif(
    comparisons = list(
      c("paln", "lp"),
      c("rp", "lp"),
      c("rp", "paln")
    ),
    annotations = c(paste("p =", signif(TukeyHSD(res.aov)$region[1:3, "p adj"], digits = 2))),
    y_position = c(2000, 2400, 2200)
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  ylab("ADCmean") +
  xlab("")

g <- ggarrange(g_fdg, g_fec, g_adc,
  nrow = 1,
  ncol = 3
)

ggsave("outputs/boxplot_raw_quants_byRegion.png",
  g,
  height = 5,
  width = 5,
  units = "in",
  dpi = 300
)

# all quants ----
quant_row_1 <- c("FDG_SUV", "FDG_NTR", "FDG_STAR", "ADC")
quant_row_2 <- c("FEC_SUV", "FEC_NTR", "FEC_STAR", "ADC_NTR")
pt_features <- c("FDG_SUV_PT", "FEC_SUV_PT", "ADC_PT")

canc <- c("cer", "endo")

plot_myCustomGG <- function(feature, show_opt_cut = F) {
  plt_data <- data %>%
    filter(CANC == canc)

  # convert hist to factor
  plt_data <- mutate(plt_data, HIST = factor(HIST, levels = c(0, 1)))

  ylabs <- c(
    "FDG_SUV_PT" = expression(FDG ~ SUV[max]),
    "FEC_SUV_PT" = expression(FEC ~ SUV[max]),
    "ADC_PT" = expression(ADC[mean]),
    "FDG_SUV" = expression(FDG ~ SUV[max]),
    "FDG_SA" = expression(FDG ~ SA ~ (mm)),
    "FDG_LA" = expression(FDG ~ LA ~ (mm)),
    "FDG_NTR" = expression(FDG ~ NTR),
    "FDG_STAR" = expression(FDG ~ STAR),
    "FDG_SNSA" = expression(FDG ~ SNSA),
    "FEC_SUV" = expression(FEC ~ SUV[max]),
    "FEC_SA" = expression(FEC ~ SA ~ (mm)),
    "FEC_LA" = expression(FEC ~ LA ~ (mm)),
    "FEC_NTR" = expression(FEC ~ NTR),
    "FEC_STAR" = expression(FEC ~ STAR),
    "FEC_SNSA" = expression(FEC ~ SNSA),
    "ADC" = expression(ADC[mean]),
    "ADC_NTR" = expression(ADC[mean] ~ NTR)
  )

  ylabs <- c(
    "FDG_SUV_PT" = "FDG SUVmax",
    "FEC_SUV_PT" = "FEC SUVmax",
    "ADC_PT" = "ADCmean",
    "FDG_SUV" = "FDG SUVmax",
    "FDG_SA" = "FDG SA (mm)",
    "FDG_LA" = "FDG LA (mm)",
    "FDG_NTR" = "FDG SUVmax NTR",
    "FDG_STAR" = "FDG STAR",
    "FDG_SNSA" = "FDG SNSA",
    "FEC_SUV" = "FEC SUVmax",
    "FEC_SA" = "FEC SA (mm)",
    "FEC_LA" = "FEC LA (mm)",
    "FEC_NTR" = "FEC SUVmax NTR",
    "FEC_STAR" = "FEC STAR",
    "FEC_SNSA" = "FEC SNSA",
    "ADC" = "ADCmean",
    "ADC_NTR" = "ADCmean NTR"
  )

  print(feature)

  # get yposition for signif bar
  signif_ypos <- max(plt_data[feature], na.rm = T) + (max(plt_data[feature], na.rm = T) * .02)

  # get n observations for each
  n_obs <- plt_data %>%
    select(any_of(feature), HIST) %>%
    filter(!is.na(unname(unlist(as.vector(.[feature]))))) %>%
    group_by(HIST) %>%
    count() %>%
    mutate(n = paste0("n=", n))

  # get yposition for n
  n_ypos <- data.frame(ypos = c(signif_ypos * 1.2, signif_ypos * 1.2))

  n_obs <- cbind(n_obs, n_ypos)

  # plot g
  g <- ggplot(
    plt_data,
    aes(
      x = HIST,
      y = !!rlang::sym(feature)
    )
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.1), alpha = .7) +
    geom_signif(
      comparisons = list(c("0", "1")),
      test = wilcox.test,
      map_signif_level = function(p) paste0("p = ", signif(p, 3)),
      y_position = c(signif_ypos)
    ) +
    ylab(ylabs[feature]) +
    xlab("") +
    scale_x_discrete(label = c("B", "M")) +
    geom_text(
      data = n_obs,
      aes(
        x = HIST,
        y = ypos,
        label = n
      )
    ) +
    theme_classic()

  if (show_opt_cut) {
    oc <- opt_cut %>%
      filter(name == ylabs[feature]) %>%
      pull(optimal_cutpoint)
    g <- g +
      geom_hline(
        yintercept = oc,
        color = "red",
        linetype = "dashed",
        linewidth = 1
      )
  }

  return(g)
}

opt_cut <- read.csv("outputs/tables/optimal_cutpoints_avg.csv")

canc <- "cer"
gg_cer_pt <- lapply(pt_features, plot_myCustomGG)
gg_cer_r1 <- lapply(quant_row_1, plot_myCustomGG)
gg_cer_r2 <- lapply(quant_row_2, plot_myCustomGG)

canc <- "endo"
gg_endo_pt <- lapply(pt_features, plot_myCustomGG)
gg_endo_r1 <- lapply(quant_row_1, plot_myCustomGG, show_opt_cut = T)
gg_endo_r2 <- lapply(quant_row_2, plot_myCustomGG, show_opt_cut = T)

## pt quants ----
g <- ggarrange(
  plotlist = c(gg_endo_pt, gg_cer_pt),
  ncol = length(pt_features),
  nrow = 2,
  labels = c("A", rep("", length(pt_features) - 1), "B")
)
ggsave("outputs/boxplot_quants_pt.png",
  g,
  width = 6,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

## endo quants ----
g <- ggarrange(
  plotlist = c(gg_endo_r1, gg_endo_r2),
  ncol = length(gg_endo_r1),
  nrow = 2,
  labels = c("A", rep("", length(gg_endo_r1) - 2), "C", "B")
)
ggsave("outputs/boxplot_quants_endo.png",
  g,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

## cer quants ----
g <- ggarrange(
  plotlist = c(gg_cer_r1, gg_cer_r2),
  ncol = length(gg_cer_r1),
  nrow = 2,
  labels = c("A", rep("", length(gg_cer_r1) - 2), "C", "B")
)
ggsave("outputs/boxplot_quants_cer.png",
  g,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)
