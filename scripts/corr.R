library(tidyverse)
library(ComplexHeatmap)
library(Hmisc)

data <- read.csv("./data/processed_data.csv")

d <- data %>%
  filter(CANC == "endo") %>%
  select(FDG_SUV, FEC_SUV, ADC, ADC_PT, FDG_SUV_PT, FEC_SUV_PT)

pairwise_corr <- rcorr(as.matrix(d), type = "spearman")
sig <- pairwise_corr$P
sig[is.na(sig)] <- 0
sig <- p.adjust(sig, method = "fdr") %>% matrix(ncol = length(colnames(pairwise_corr$P)))

ylabs <- c(
  "FDG_SUV_PT" = "FDG SUVmax (PT)",
  "FEC_SUV_PT" = "FEC SUVmax (PT)",
  "ADC_PT" = "ADCmean (PT)",
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

colnames(pairwise_corr$r) <- ylabs[as.character(colnames(pairwise_corr$r))]
rownames(pairwise_corr$r) <- ylabs[as.character(rownames(pairwise_corr$r))]

col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

hm <- Heatmap(pairwise_corr$r,
  col = col_fun,
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
  width = 6,
  height = 5.5,
  units = "in",
  res = 300
)
draw(hm)
dev.off()
