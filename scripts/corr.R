library(tidyverse)
library(ComplexHeatmap)
library(Hmisc)

source("scripts/tex.R")

data <- read.csv("./data/processed_data.csv")

d <- data %>%
  filter(CANC == "endo") %>%
  select(FDG_SUV, FEC_SUV, ADC, ADC_PT, FDG_SUV_PT, FEC_SUV_PT)

pairwise_corr <- rcorr(as.matrix(d), type = "spearman")
sig <- pairwise_corr$P
sig[is.na(sig)] <- 0
sig <- p.adjust(sig, method = "fdr") %>% matrix(ncol = length(colnames(pairwise_corr$P)))

ylabs <- c(
  "FDG_SUV_PT" = tex$fdg_suvmax_pt,
  "FEC_SUV_PT" = tex$fec_suvmax_pt,
  "ADC_PT" = tex$mri_adc_pt,
  "FDG_SUV" = tex$fdg_suvmax,
  "FDG_SA" = "FDG SA (mm)",
  "FDG_LA" = "FDG LA (mm)",
  "FDG_NTR" = tex$fdg_suvmax_ntr,
  "FDG_STAR" = tex$fdg_suvmax_star,
  "FDG_SNSA" = "FDG SNSA",
  "FEC_SUV" = tex$fec_suvmax,
  "FEC_SA" = "FEC SA (mm)",
  "FEC_LA" = "FEC LA (mm)",
  "FEC_NTR" = tex$fec_suvmax_ntr,
  "FEC_STAR" = tex$fec_suvmax_star,
  "FEC_SNSA" = "FEC SNSA",
  "ADC" = tex$mri_adc,
  "ADC_NTR" = tex$mri_adc_ntr
)

col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

hm <- Heatmap(pairwise_corr$r,
  col = col_fun,
  row_labels = ylabs[as.character(rownames(pairwise_corr$r))],
  column_labels = ylabs[as.character(colnames(pairwise_corr$r))],
  heatmap_legend_param = list(title = "R"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (i != j) {
      if (sig[i, j] < 0.001) {
        text <- paste0(signif(pairwise_corr$r[i, j], 2), "***")
        grid.text(text, x, y, gp = gpar(fontsize = 10))
      } else if (sig[i, j] < 0.01) {
        text <- paste0(signif(pairwise_corr$r[i, j], 2), "**")
        grid.text(text, x, y, gp = gpar(fontsize = 10))
      } else if (sig[i, j] < 0.05) {
        text <- paste0(signif(pairwise_corr$r[i, j], 2), "*")
        grid.text(text, x, y, gp = gpar(fontsize = 10))
      }
    } else {
      grid.rect(x, y, w, h, gp = gpar(fill = "grey", col = NA))
    }
  }
)

png("outputs/pairwise_corr_raw.png",
  width = 7,
  height = 6.5,
  units = "in",
  res = 300
)
draw(hm)
dev.off()
