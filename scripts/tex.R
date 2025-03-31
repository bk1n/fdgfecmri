library(latex2exp)

tex <- list(
    "fdg" = TeX("$\\lbrack^{18}F\\rbrack$FDG PET/CT"),
    "fec" = TeX("$\\lbrack^{18}F\\rbrack$FEC PET/CT"),
    "mri" = TeX("DW-MRI"),
    "fdg_suvmax" = TeX("$\\lbrack^{18}F\\rbrack$FDG PET/CT $SUV_{max}$"),
    "fec_suvmax" = TeX("$\\lbrack^{18}F\\rbrack$FEC PET/CT $SUV_{max}$"),
    "mri_adc" = TeX("DW-MRI $ADC_{mean}$"),
    "fdg_suvmax_pt" = TeX("$\\lbrack^{18}F\\rbrack$FDG PET/CT $SUV_{max}$ (PT)"),
    "fec_suvmax_pt" = TeX("$\\lbrack^{18}F\\rbrack$FEC PET/CT $SUV_{max}$ (PT)"),
    "mri_adc_pt" = TeX("DW-MRI $ADC_{mean}$ (PT)"),
    "fdg_suvmax_ntr" = TeX("$\\lbrack^{18}F\\rbrack$FDG PET/CT $SUV_{max}$ NTR"),
    "fec_suvmax_ntr" = TeX("$\\lbrack^{18}F\\rbrack$FEC PET/CT $SUV_{max}$ NTR"),
    "mri_adc_ntr" = TeX("DW-MRI $ADC_{mean}$ NTR"),
    "fdg_suvmax_star" = TeX("$\\lbrack^{18}F\\rbrack$FDG PET/CT $SUV_{max}$ STAR"),
    "fec_suvmax_star" = TeX("$\\lbrack^{18}F\\rbrack$FEC PET/CT $SUV_{max}$ STAR")
)

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

ylabs_norm <- c(
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

# tex <- list(
#     "fdg" = TeX("$\\lbrack^{18}F\\rbrack$FDG PET/CT"),
#     "fec" = TeX("$\\lbrack^{18}F\\rbrack$FEC PET/CT"),
#     "mri" = TeX("DW-MRI"),
#     "fdg_suvmax" = TeX("$SUV^{FDG}_{max}$"),
#     "fec_suvmax" = TeX("$SUV^{FEC}_{max}$"),
#     "mri_adc" = TeX("$ADC_{mean}$"),
#     "fdg_suvmax_pt" = TeX("$SUV^{FDG}_{max}$ (PT)"),
#     "fec_suvmax_pt" = TeX("$SUV^{FEC}_{max}$ (PT)"),
#     "mri_adc_pt" = TeX("$ADC_{mean}$ (PT)"),
#     "fdg_suvmax_ntr" = TeX("$SUV^{FDG}_{max}$ NTR"),
#     "fec_suvmax_ntr" = TeX("$SUV^{FEC}_{max}$ NTR"),
#     "mri_adc_ntr" = TeX("$ADC_{mean}$ NTR"),
#     "fdg_suvmax_star" = TeX("$SUV^{FDG}_{max}$ STAR"),
#     "fec_suvmax_star" = TeX("$SUV^{FEC}_{max}$ STAR")
# )
