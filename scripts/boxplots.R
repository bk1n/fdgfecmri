library(tidyverse)
library(ggsignif)
library(ggpubr)
library(latex2exp)

source("scripts/tex.R")

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
  ylab(tex$fdg_suvmax) +
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
    test = wilcox.test,
    map_signif_level = function(p) paste0("p = ", signif(p, 3))
  ) +
  scale_x_discrete(label = c("Cervical", "Endometrial")) +
  ylab(tex$fec_suvmax) +
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
    test = wilcox.test,
    map_signif_level = function(p) paste0("p = ", signif(p, 3))
  ) +
  scale_x_discrete(label = c("Cervical", "Endometrial")) +
  ylab(tex$mri_adc) +
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
  ylab(tex$fdg_suvmax) +
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
  ylab(tex$fec_suvmax) +
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
  ylab(tex$mri_adc) +
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

  print(feature)

  # get yposition for signif bar
  signif_ypos <- max(plt_data[feature], na.rm = T) * 1.08

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
      xmin = 1,
      xmax = 2,
      annotations = paste0("p=", signif(padj[[feature]], 3)),
      y_position = c(signif_ypos)
    ) +
    ylab(ylabs_short[feature]) +
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
      select(measure, beta, median_oc, cancer) %>%
      filter(measure == feature, beta %in% c(0.5, 1, 2), cancer == canc) %>%
      mutate(beta = as.factor(beta))
    g <- g +
      geom_hline(
        data = oc,
        aes(yintercept = median_oc, color = beta, linetype = beta),
        # linetype = "dashed",
        alpha = .7,
        linewidth = 1
      ) +
      labs(color = TeX("$\\beta"), linetype = TeX("$\\beta"))
  }

  return(g)
}

run_wilcox <- function(feature) {
  plt_data <- data %>%
    filter(CANC == canc)
  w <- wilcox.test(plt_data[[feature]] ~ HIST, plt_data)
  return(w$p.value)
}
opt_cut <- read.csv("outputs/tables/diagnostic_performance.csv")

canc <- "cer"
r1 <- lapply(quant_row_1, run_wilcox)
names(r1) <- quant_row_1
r2 <- lapply(quant_row_2, run_wilcox)
names(r2) <- quant_row_2
pt <- lapply(pt_features, run_wilcox)
names(pt) <- pt_features

pval <- c(r1, r2, pt)
padj <- p.adjust(pval, method = "BH")

gg_cer_pt <- lapply(pt_features, plot_myCustomGG)
gg_cer_r1 <- lapply(quant_row_1, plot_myCustomGG)
gg_cer_r2 <- lapply(quant_row_2, plot_myCustomGG)

canc <- "endo"
r1 <- lapply(quant_row_1, run_wilcox)
names(r1) <- quant_row_1
r2 <- lapply(quant_row_2, run_wilcox)
names(r2) <- quant_row_2
pt <- lapply(pt_features, run_wilcox)
names(pt) <- pt_features

pval <- c(r1, r2, pt)
padj <- p.adjust(pval, method = "BH")

gg_endo_pt <- lapply(pt_features, plot_myCustomGG)
gg_endo_r1 <- lapply(quant_row_1, plot_myCustomGG, show_opt_cut = T)
gg_endo_r2 <- lapply(quant_row_2, plot_myCustomGG, show_opt_cut = T)

## pt quants ----
g <- ggarrange(
  plotlist = c(gg_endo_pt, gg_cer_pt),
  ncol = length(pt_features),
  nrow = 2,
  labels = c("a", rep("", length(pt_features) - 1), "b")
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
  align = "v",
  labels = c("a", rep("", length(gg_endo_r1) - 2), "c", "b"),
  common.legend = T
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
  labels = c("a", rep("", length(gg_cer_r1) - 2), "c", "b")
)
ggsave("outputs/boxplot_quants_cer.png",
  g,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)
