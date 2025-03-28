# estimates diagnostic performance of measures by all quantitative measures
# compares to visual assessment by expert readers

library(tidyverse)
library(cutpointr)
library(ggpubr)
library(caret)
library(boot)

source("scripts/tex.R")

# filter data ----
data <- read.csv("./data/processed_data.csv")

data <- data %>%
    select(
        FDG_SUV, FDG_STAR, FDG_NTR,
        FEC_SUV, FEC_STAR, FEC_NTR,
        ADC, ADC_NTR,
        VIS_FDG, VIS_FEC, VIS_MRI,
        HIST, region, CANC
    )

data <- data %>% dplyr::filter(CANC == "endo") # endo only
# convert visual measures to binary
data <- data %>%
    mutate(across(VIS_FDG:VIS_MRI, .fn = function(x) if_else(x >= 5, 1, 0)))

measures <- c("FDG_SUV", "FDG_STAR", "FDG_NTR", "FEC_SUV", "FEC_STAR", "FEC_NTR", "ADC", "ADC_NTR")

cv.fun <- function(data, train_idx, measure, visual_measure) {
    m <- measure
    vm <- visual_measure

    train <- data[train_idx, ]
    test <- data[-train_idx, ]

    train_x <- train[[m]]
    train_y <- train[["HIST"]]

    dir <- if (grepl("ADC", m)) "<=" else ">="

    oc <- cutpointr(
        x = train_x,
        class = train_y,
        na.rm = T,
        method = maximize_metric,
        metric = F1_score,
        pos_class = 1,
        neg_class = 0,
        direction = dir
    )
    oc <- oc$optimal_cutpoint

    test_y <- test[["HIST"]]

    test_oc <- if (dir == ">=") as.numeric(test[[m]] >= oc) else as.numeric(test[[m]] <= oc)
    test_vm <- test[[vm]]

    c(
        optimal_cutpoint = oc,
        visual = mean(test_y == test_vm),
        quant = mean(test_y == test_oc)
    )
}

r <- lapply(measures, function(m) {
    # filter data to subset
    vm <- if (grepl("FDG", m)) "VIS_FDG" else if (grepl("FEC", m)) "VIS_FEC" else "VIS_MRI"
    x <- data %>% dplyr::filter(!is.na(!!rlang::sym(m)) & !is.na(!!rlang::sym(vm)))

    # 100x 10-fold CV
    cv <- caret::createMultiFolds(y = x[["HIST"]], k = 4, times = 100)
    z <- lapply(cv, function(cv_idx) {
        cv.fun(data = x, train_idx = cv_idx, measure = m, visual_measure = vm)
    })
    z <- as.data.frame(do.call(rbind, z))
    z$measure <- m
    z
})
res <- do.call(rbind, r)
d <- res %>%
    group_by(measure) %>%
    mutate(diff = visual - quant) %>%
    filter(measure == "FDG_SUV")


x <- d$visual
x <- rnorm(100000)
r <- (length(x) / 2) - (1.96 * sqrt(length(x)) / 2)
s <- 1 + (length(x) / 2) + (1.96 * sqrt(length(x)) / 2)

lci <- x[order(x)][round(r)]
uci <- x[order(x)][round(s)]
median(x)

hist(x)
abline(v = c(lci, median(x), uci))

100 * ((1 - 0.95) / 2)

hist(res$quant[res$measure == "FDG_SUV"])
hist(res$visual[res$measure == "FDG_SUV"])



# get bootstrap dataset
# for each measure
# determine optimal cutpoints
# calculate performance metrics of cutpoints on OOB samples
# calculate performance metrics

# opt cut func ----
get_opt_cut <- function(data, measure) {
    # filter df to patients w both visual and quants


    print(paste("Getting cutpoint for:", measure))

    x <- df[[measure]]
    class <- df$HIST
    if (grepl("ADC", measure)) {
        dir <- "<="
    } else {
        dir <- ">="
    }
    opt_cut <- cutpointr(
        x = x,
        class = class,
        na.rm = T,
        method = maximize_metric,
        metric = F1_score,
        pos_class = 1,
        neg_class = 0,
        direction = dir,
        boot_runs = 2000,
        boot_stratify = T
    )
    summary.opt_cut <- summary(opt_cut)

    boot_data <- summary.opt_cut$cutpointr[[1]]$boot[[1]]
    boot_data$name <- measure
    boot_data <- relocate(boot_data, "name")
    boot_data <- select(boot_data, -where(is.list))

    return(list(
        opt_cut = opt_cut,
        summary = summary.opt_cut,
        boot_data = boot_data,
        df = df
    ))
}

# run opt cut ----
run_opt_cut <- F
if (run_opt_cut) {
    set.seed(42)

    run_cols <- c(
        "FDG_SUV", "FDG_STAR", "FDG_NTR",
        "FEC_SUV", "FEC_STAR", "FEC_NTR",
        "ADC", "ADC_NTR"
    )

    cutpoint_res <- lapply(run_cols, get_opt_cut, data = data)
    names(cutpoint_res) <- run_cols

    plt_df <- do.call(rbind, lapply(cutpoint_res, function(r) r$boot_data))

    ## save opt cut res ----
    write.csv(plt_df, "outputs/tables/optimal_cutpoints.csv", row.names = F)
} else {
    plt_df <- read.csv("outputs/tables/optimal_cutpoints.csv")
}

## plt opt cut res boxplots ----
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

# format plt_df
plt_df_tidy <- plt_df %>%
    # select measures
    select(name, optimal_cutpoint, AUC_oob, sensitivity_oob, specificity_oob) %>%
    # convert names
    mutate(name = ylabs_norm[name]) %>%
    # get measure type
    mutate(measure_type = case_when(
        grepl("FDG", name) ~ "FDG-PET/CT",
        grepl("FEC", name) ~ "FEC-PET/CT",
        grepl("ADC", name) ~ "DW-MRI",
        grepl("LR", name) ~ "LR"
    )) %>%
    # mutate measure type
    mutate(
        measure_type = factor(measure_type, levels = c("FDG-PET/CT", "FEC-PET/CT", "DW-MRI")),
        name = factor(name, levels = c(
            "FDG SUVmax", "FDG STAR", "FDG SUVmax NTR",
            "FEC SUVmax", "FEC STAR", "FEC SUVmax NTR",
            "ADCmean", "ADCmean NTR"
        ))
    )

# get mean of each measure
mean_df <- plt_df %>%
    group_by(name) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
    ungroup()

mean_df_tidy <- plt_df_tidy %>%
    group_by(name, measure_type) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
    ungroup()

# write mean results
write.csv(mean_df_tidy, "outputs/tables/optimal_cutpoints_avg.csv", row.names = F)

# plot boxplots
g_auc <- ggplot(
    plt_df_tidy,
    aes(
        x = name,
        y = AUC_oob
    )
) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(
        data = mean_df,
        size = 3,
        shape = 3,
        color = "red"
    ) +
    theme_classic() +
    xlab("") +
    ylab("Area Under Curve (AUC)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.x = element_blank()) +
    facet_grid(~measure_type, space = "free", scales = "free_x")

g_sens <- ggplot(
    plt_df_tidy,
    aes(
        x = name,
        y = sensitivity_oob
    )
) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(
        data = mean_df,
        size = 3,
        shape = 3,
        color = "red"
    ) +
    theme_classic() +
    xlab("") +
    ylab("Sensitivity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.x = element_blank()) +
    facet_grid(~measure_type, space = "free", scales = "free_x")

g_spec <- ggplot(
    plt_df_tidy,
    aes(
        x = name,
        y = specificity_oob
    )
) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(
        data = mean_df,
        size = 3,
        shape = 3,
        color = "red"
    ) +
    theme_classic() +
    xlab("") +
    ylab("Specificity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(~measure_type, space = "free", scales = "free_x")

g <- ggarrange(g_auc, g_sens, g_spec, nrow = 3, heights = c(0.8, 0.8, 1))

# refit optimal cutpoints & calculate diagnostic performance ----
get_diagnostic_performance <- function(measure,
                                       data) {
    # given measure, get optimal cutpoint
    # mean_df = optimal cutpoints
    cut <- mean_df$optimal_cutpoint[mean_df$name == measure]

    # get vis measure
    if (!(measure %in% c("VIS_FDG", "VIS_FEC", "VIS_MRI"))) {
        vis_measure <- if (grepl("FDG", measure)) "VIS_FDG" else if (grepl("FEC", measure)) "VIS_FEC" else "VIS_MRI"
        df <- data %>%
            dplyr::filter(!is.na(!!rlang::sym(measure)) & !is.na(!!rlang::sym(vis_measure)))
    } else {
        vis_measure <- measure
        alt_measure <- if (grepl("FDG", measure)) "FDG_SUV" else if (grepl("FEC", measure)) "FEC_SUV" else "ADC"
        df <- data %>%
            dplyr::filter(!is.na(!!rlang::sym(vis_measure)) & !is.na(!!rlang::sym(alt_measure)))
    }

    # get x (measure) and y (HIST)
    x <- df[[measure]]
    y <- df[["HIST"]]

    # convert x
    if (!(measure %in% c("VIS_FDG", "VIS_FEC", "VIS_MRI"))) {
        x <- if (vis_measure != "VIS_MRI") as.numeric(x >= cut) else as.numeric(x <= cut)
    } else {
        x <- as.numeric(x >= 5)
    }

    cm <- table(measure = x, hist = y)

    # hist = cols, measure = rows
    stopifnot(length(x) == length(y))

    num_regions <- length(x)
    tp <- cm["1", "1"]
    fn <- cm["0", "1"]
    fp <- cm["1", "0"]
    tn <- cm["0", "0"]

    sens <- tp / (tp + fn)
    spec <- tn / (tn + fp)
    ppv <- tp / (tp + fp)
    npv <- tn / (tn + fn)
    f1 <- 2 * tp / (2 * tp + fp + fn)

    return(list(
        num_regions = num_regions,
        tp = tp,
        fn = fn,
        tn = tn,
        fp = fp,
        sens = sens,
        spec = spec,
        ppv = ppv,
        npv = npv,
        f1 = f1
    ))
}

measures <- c(
    "FDG_SUV", "FDG_STAR", "FDG_NTR",
    "FEC_SUV", "FEC_STAR", "FEC_NTR",
    "ADC", "ADC_NTR",
    "VIS_FDG", "VIS_FEC", "VIS_MRI"
)

performance <- lapply(measures, get_diagnostic_performance, data = data)
names(performance) <- measures

performance_df <- do.call(rbind, performance)

write.csv(performance_df, "outputs/tables/diagnostic_performance_refit.csv", row.names = T)

# mcnemars test ----
mcnemar_test <- function(measure, data) {
    # given measure, get optimal cutpoint
    # mean_df = optimal cutpoints
    cut <- mean_df$optimal_cutpoint[mean_df$name == measure]

    # get vis measure
    vis_measure <- if (grepl("FDG", measure)) "VIS_FDG" else if (grepl("FEC", measure)) "VIS_FEC" else "VIS_MRI"
    df <- data %>%
        dplyr::filter(!is.na(!!rlang::sym(measure)) & !is.na(!!rlang::sym(vis_measure)))

    # get x (measure) and y (HIST)
    q <- df[[measure]]
    v <- df[[vis_measure]]
    y <- df[["HIST"]]

    # convert q & v
    q <- if (vis_measure != "VIS_MRI") as.numeric(q >= cut) else as.numeric(q <= cut)
    v <- as.numeric(v >= 5)

    pos <- list(
        q = q[y == 1],
        v = v[y == 1],
        y = y[y == 1]
    )

    neg <- list(
        q = q[y == 0],
        v = v[y == 0],
        y = y[y == 0]
    )

    pos_cm <- table(
        quant = (pos$q == pos$y),
        vis = (pos$v == pos$y)
    )

    neg_cm <- table(
        quant = (neg$q == neg$y),
        vis = (neg$v == neg$y)
    )


    return(
        list(
            measure = measure,
            vis_measure = vis_measure,
            num_regions = length(q),
            sens_mcnemar = mcnemar.test(pos_cm)$p.value,
            spec_mcnemar = mcnemar.test(neg_cm)$p.value
        )
    )
}

measures <- c(
    "FDG_SUV", "FDG_STAR", "FDG_NTR",
    "FEC_SUV", "FEC_STAR", "FEC_NTR",
    "ADC", "ADC_NTR"
)

mcnemar_res <- lapply(measures, mcnemar_test, data = data)
names(mcnemar_res) <- measures

mcnemar_res_df <- do.call(rbind, lapply(mcnemar_res, data.frame))

mcnemar_res_df$sens_padj <- p.adjust(mcnemar_res_df$sens_mcnemar)
mcnemar_res_df$spec_padj <- p.adjust(mcnemar_res_df$spec_mcnemar)

## save results
write.csv(mcnemar_res_df, "outputs/tables/diagnostic_performance_mcnemar_results.csv", row.names = F)
