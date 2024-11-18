# calculate ICC for each method
library(tidyverse)
library(irr)

data <- as_tibble(read.csv("./data/fdgpet_data.csv", header = T))

data_endo <- data %>% dplyr::filter(CANC_ENDO == 1)
data_cer <- data %>% dplyr::filter(CANC_CER == 1)

# endo ----
# fdg
lp <- data_endo[c("FDG_SUV_LP_1_CR", "FDG_SUV_LP_2_CR", "FDG_SUV_LP")]
rp <- data_endo[c("FDG_SUV_RP_1_CR", "FDG_SUV_RP_2_CR", "FDG_SUV_RP")]
pa <- data_endo[c("FDG_SUV_PALN_1_CR", "FDG_SUV_PALN_2_CR", "FDG_SUV_PALN")]

colnames(lp) <- colnames(rp) <- colnames(pa) <- NULL

x <- rbind(lp, rp, pa)
x <- log(x)
hist(as.matrix(x))

fdg_icc <- icc(x)

# fec
lp <- data_endo[c("FEC_SUV_LP_1_CR", "FEC_SUV_LP_2_CR", "FEC_SUV_LP")]
rp <- data_endo[c("FEC_SUV_RP_1_CR", "FEC_SUV_RP_2_CR", "FEC_SUV_RP")]
pa <- data_endo[c("FEC_SUV_PALN_1_CR", "FEC_SUV_PALN_2_CR", "FEC_SUV_PALN")]

colnames(lp) <- colnames(rp) <- colnames(pa) <- NULL

x <- rbind(lp, rp, pa)
x <- log(x)
hist(as.matrix(x))

fec_icc <- icc(x)

# mri
lp <- data_endo[c("ADC_LP_1_CR", "ADC_LP_2_CR")]
rp <- data_endo[c("ADC_RP_1_CR", "ADC_RP_2_CR")]
pa <- data_endo[c("ADC_PALN_1_CR", "ADC_PALN_2_CR")]

colnames(lp) <- colnames(rp) <- colnames(pa) <- NULL

x <- rbind(lp, rp, pa)
adc_icc <- icc(x)

# pt
pt <- data_endo[c("ADC_PT_1_CR", "ADC_PT_2_CR")]
pt_endo_icc <- icc(pt)


# df
df_endo <- data.frame(
    Method = c("FDG PET/CT", "FEC PET/CT", "DW-MRI", "DW-MRI (PT)"),
    Cancer = "Endometrial",
    ICC = c(signif(fdg_icc$value, 3), signif(fec_icc$value, 3), signif(adc_icc$value, 3), signif(pt_endo_icc$value, 3)),
    CI95 = c(
        paste0(signif(fdg_icc$lbound, 3), "-", signif(fdg_icc$ubound, 3)),
        paste0(signif(fec_icc$lbound, 3), "-", signif(fec_icc$ubound, 3)),
        paste0(signif(adc_icc$lbound, 3), "-", signif(adc_icc$ubound, 3)),
        paste0(signif(pt_endo_icc$lbound, 3), "-", signif(pt_endo_icc$ubound, 3))
    )
)

# cervical ----
# fdg
lp <- data_cer[c("FDG_SUV_LP_1_CR", "FDG_SUV_LP_2_CR", "FDG_SUV_LP")]
rp <- data_cer[c("FDG_SUV_RP_1_CR", "FDG_SUV_RP_2_CR", "FDG_SUV_RP")]
pa <- data_cer[c("FDG_SUV_PALN_1_CR", "FDG_SUV_PALN_2_CR", "FDG_SUV_PALN")]

colnames(lp) <- colnames(rp) <- colnames(pa) <- NULL

x <- rbind(lp, rp, pa)
x <- log(x)
hist(as.matrix(x))

fdg_icc <- icc(x)

# fec
lp <- data_cer[c("FEC_SUV_LP_1_CR", "FEC_SUV_LP_2_CR", "FEC_SUV_LP")]
rp <- data_cer[c("FEC_SUV_RP_1_CR", "FEC_SUV_RP_2_CR", "FEC_SUV_RP")]
pa <- data_cer[c("FEC_SUV_PALN_1_CR", "FEC_SUV_PALN_2_CR", "FEC_SUV_PALN")]

colnames(lp) <- colnames(rp) <- colnames(pa) <- NULL

x <- rbind(lp, rp, pa)
x <- log(x)
hist(as.matrix(x))

fec_icc <- icc(x)

# mri
lp <- data_cer[c("ADC_LP_1_CR", "ADC_LP_2_CR")]
rp <- data_cer[c("ADC_RP_1_CR", "ADC_RP_2_CR")]
pa <- data_cer[c("ADC_PALN_1_CR", "ADC_PALN_2_CR")]

colnames(lp) <- colnames(rp) <- colnames(pa) <- NULL

x <- rbind(lp, rp, pa)
adc_icc <- icc(x)

# pt
pt <- data_cer[c("ADC_PT_1_CR", "ADC_PT_2_CR")]
pt_cer_icc <- icc(pt)

# df
df_cer <- data.frame(
    Method = c("FDG PET/CT", "FEC PET/CT", "DW-MRI", "DW-MRI (PT)"),
    Cancer = "Cervical",
    ICC = c(signif(fdg_icc$value, 3), signif(fec_icc$value, 3), signif(adc_icc$value, 3), signif(pt_cer_icc$value, 3)),
    CI95 = c(
        paste0(signif(fdg_icc$lbound, 3), "-", signif(fdg_icc$ubound, 3)),
        paste0(signif(fec_icc$lbound, 3), "-", signif(fec_icc$ubound, 3)),
        paste0(signif(adc_icc$lbound, 3), "-", signif(adc_icc$ubound, 3)),
        paste0(signif(pt_cer_icc$lbound, 3), "-", signif(pt_cer_icc$ubound, 3))
    )
)

df <- rbind(df_endo, df_cer)
write.csv(df, "outputs/tables/icc.csv")
