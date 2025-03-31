# estimates diagnostic performance of measures by all quantitative measures
# compares to visual assessment by expert readers

library(tidyverse)
library(cutpointr)
library(ggpubr)
library(caret)
library(boot)

source("scripts/tex.R")

out_path <- "outputs/"

run_pipeline <- F

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

# helper funcs ----
metric_translator <- c(
    "f1" = "F1-score",
    "fbeta" = TeX("$F_{\\beta}$"),
    "npv" = "NPV",
    "ppv" = "PPV",
    "sens" = "Sensitivity",
    "spec" = "Specificity"
)

# custom metric function
fbeta <- function(tp, fp, tn, fn, beta, ...) {
    se <- tp / (tp + fn)
    ppv <- tp / (tp + fp)
    fbeta <- (1 + beta^2) * (2 * se * ppv) / (beta^2 * ppv + se)
    fbeta
}

calculate_metrics <- function(hist, preds, beta = 1) {
    cm <- table(
        hist = factor(hist, levels = c(0, 1)),
        pred = factor(preds, levels = c(0, 1))
    )
    tp <- cm[2, 2]
    fp <- cm[1, 2]
    tn <- cm[1, 1]
    fn <- cm[2, 1]
    sens <- tp / (tp + fn)
    spec <- tn / (tn + fp)
    ppv <- tp / (tp + fp)
    npv <- tn / (tn + fn)
    acc <- (tp + tn) / (tp + tn + fp + fn)
    f1 <- (2 * sens * ppv) / (ppv + sens)
    fbeta <- (1 + beta^2) * ((sens * ppv) / (beta^2 * ppv + sens))
    data.frame(
        tp, fn, tn, fp,
        sens, spec, ppv, npv, f1, fbeta
    )
}

# for a given dataset with train_idx, split data and run cutpoint optimisation on measure
# returns predictions
cv.fun <- function(data, train_idx, measure, visual_measure, beta) {
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
        metric = fbeta, # TODO - select appropriate metric
        pos_class = 1,
        neg_class = 0,
        direction = dir,
        beta = beta
    )
    oc <- oc$optimal_cutpoint

    test_y <- test[["HIST"]]

    test_oc <- if (dir == ">=") as.numeric(test[[m]] >= oc) else as.numeric(test[[m]] <= oc)
    test_vm <- test[[vm]]

    c(oc = oc, hist = test_y, visual = test_vm, quant = test_oc)
}

# apply cv.fun across a set of measures
# returns binded data.frame of predictions
cv <- function(measures, beta = 1) {
    r <- lapply(measures, function(m) {
        # filter data to subset
        vm <- if (grepl("FDG", m)) "VIS_FDG" else if (grepl("FEC", m)) "VIS_FEC" else "VIS_MRI"
        x <- data %>% dplyr::filter(!is.na(!!rlang::sym(m)) & !is.na(!!rlang::sym(vm)))

        # 100x 10-fold CV
        cv <- caret::createFolds(y = x[["HIST"]], k = nrow(x), returnTrain = T) # LOO
        z <- lapply(cv, function(cv_idx) {
            cv.fun(data = x, train_idx = cv_idx, measure = m, visual_measure = vm, beta = beta)
        })
        z <- as.data.frame(do.call(rbind, z))
        z$measure <- m
        z$beta <- beta
        z
    })
    res <- do.call(rbind, r)
    res
}

# bootstrap predictions to get confidence intervals for different metrics
# x should be a data.frame with hist, quant and beta columns (for fbeta)
boot.fun <- function(x, i) {
    x <- x[i, ]
    qm_r <- unlist(calculate_metrics(x$hist, x$quant, unique(x$beta))[1, , drop = T])
    vm_r <- unlist(calculate_metrics(x$hist, x$visual, unique(x$beta))[1, , drop = T])
    names(qm_r) <- paste0("qm.", names(qm_r))
    names(vm_r) <- paste0("vm.", names(vm_r))
    c(qm_r, vm_r)
}

get.ci <- function(x, w, method = "perc") {
    if (all(is.na(x$t[, w])) | all(x$t[, w] == 1, na.rm = T) | all(x$t[, w] == 0, na.rm = T)) {
        return(data.frame(lwr = NA, upr = NA))
    }
    b1 <- boot.ci(x, index = w, type = method)
    ## extract info for all CI types
    tab <- t(sapply(b1[-(1:3)], function(x) tail(c(x), 2)))
    ## combine with metadata: CI method, index
    tab <- as.data.frame(tab)
    colnames(tab) <- c("lwr", "upr")
    tab
}

bootstrap <- function(x) {
    bs <- boot::boot(x, boot.fun,
        R = 1000,
        sim = "ordinary"
    )

    ci <- do.call(rbind, lapply(1:length(bs$t0), function(i) get.ci(bs, w = i, method = "perc")))
    ci$estimate <- bs$t0
    ci$metric <- names(bs$t0)
    ci$num_regions <- nrow(x)
    ci$measure <- unique(x$measure)
    ci$beta <- unique(x$beta)

    ci
}

# run pipeline

# For each sample j in {1:,...,n}:
#   1. Identify cut-off c_j on all samples except j: {1,...,j-1,j+1,...,n}
#   2. If x_j >= O_j, Ŷ = 1 otherwise Ŷ = 0 (or <= for ADC measures)
#   3. Repeat for all samples to obtain predictions Ŷ = {Ŷ_1,...,Ŷ_n}
# Bootstrap Ŷ 1000 times to estimate confidence intervals:
#   1. From original predictions {Ŷ₁,...,Ŷₙ}, draw n samples with replacement
#   2. Calculate performance metrics on this bootstrap sample
#   3. Repeat 1000 times to generate empirical sampling distribution
#   4. Extract 2.5th and 97.5th percentiles for 95% CI of each metric
# Compare Ŷ to Y (predicted vs. actual values)
# Calculate performance metrics:
#   Sensitivity = P(Ŷ = 1 | Y = 1) = TP/(TP+FN)
#   Specificity = P(Ŷ = 0 | Y = 0) = TN/(TN+FP)
#   PPV = P(Y = 1 | Ŷ = 1) = TP/(TP+FP)
#   NPV = P(Y = 0 | Ŷ = 0) = TN/(TN+FN)
#   Fbeta = 2*(Sens*PPV)/(Sens+PPV) = 2TP/(2TP+FP+FN)
run <- function(measures, beta = 1) {
    preds <- cv(measures, beta = beta)

    oc <- preds %>%
        group_by(measure) %>%
        summarise(mean_oc = mean(oc), median_oc = median(oc))

    bs.metrics <- lapply(split(preds, preds$measure), function(r) {
        bs <- bootstrap(r)
        bs
    })
    bs.metrics <- do.call(rbind, bs.metrics)
    rownames(bs.metrics) <- NULL

    metrics <- bs.metrics %>%
        dplyr::select(estimate, lwr, upr, metric, measure, num_regions, beta) %>%
        pivot_wider(names_from = metric, values_from = c(estimate, lwr, upr))

    # McNemar's test on predictions vs actual
    mcnemar <- lapply(split(preds, preds$measure), function(r) {
        pos_cm <- table(
            qm = factor(as.numeric(r$quant[r$hist == 1] == 1), levels = c(0, 1)),
            vm = factor(as.numeric(r$visual[r$hist == 1] == 1), levels = c(0, 1))
        )
        neg_cm <- table(
            qm = factor(as.numeric(r$quant[r$hist == 0] == 0), levels = c(0, 1)),
            vm = factor(as.numeric(r$visual[r$hist == 0] == 0), levels = c(0, 1))
        )
        data.frame(
            sens_p = mcnemar.test(pos_cm)$p.value,
            spec_p = mcnemar.test(neg_cm)$p.value
        )
    })
    mcnemar <- do.call(rbind, mcnemar)
    mcnemar$sens_padj <- p.adjust(mcnemar$sens_p, method = "fdr")
    mcnemar$spec_padj <- p.adjust(mcnemar$spec_p, method = "fdr")
    mcnemar$measure <- rownames(mcnemar)

    metrics <- inner_join(metrics, mcnemar, by = "measure")
    metrics <- inner_join(oc, metrics, by = "measure")
    metrics
}


# run pipeline, get full results ----
measures <- c("FDG_SUV", "FDG_STAR", "FDG_NTR", "FEC_SUV", "FEC_STAR", "FEC_NTR", "ADC", "ADC_NTR")
if (run_pipeline) {
    beta_res <- lapply(2^seq(-4, 5), run, measures = measures)
    beta_res <- do.call(rbind, beta_res)
    beta_res <- beta_res %>%
        mutate(measure = factor(measure, levels = measures)) %>%
        arrange(measure)

    write.csv(beta_res, file.path(out_path, "tables", "diagnostic_performance.csv"))
} else {
    beta_res <- read.csv(file.path(out_path, "tables", "diagnostic_performance.csv"))
}

hist(log2(beta_res$beta))

# prettify results
x <- beta_res %>%
    filter(beta == 1)

## qm
qm <- x %>%
    select(estimate_qm.sens:estimate_qm.fbeta) %>%
    mutate(across(everything(), function(x) gsub(" ", "", format(round(x * 100, 1), nsmall = 1)))) %>%
    as.data.frame()
qm.lwr <- x %>%
    select(lwr_qm.sens:lwr_qm.fbeta) %>%
    mutate(across(everything(), function(x) gsub(" ", "", format(round(x * 100, 1), nsmall = 1)))) %>%
    as.data.frame()
qm.upr <- x %>%
    select(upr_qm.sens:upr_qm.fbeta) %>%
    mutate(across(everything(), function(x) gsub(" ", "", format(round(x * 100, 1), nsmall = 1)))) %>%
    as.data.frame()

for (i in 1:nrow(qm)) {
    for (j in 1:ncol(qm)) {
        qm[i, j] <- paste0(qm[i, j], " (", qm.lwr[i, j], "-", qm.upr[i, j], ")")
    }
}

## vm
vm <- x %>%
    select(estimate_vm.sens:estimate_vm.fbeta) %>%
    mutate(across(everything(), function(x) gsub(" ", "", format(round(x * 100, 1), nsmall = 1)))) %>%
    as.data.frame()
vm.lwr <- x %>%
    select(lwr_vm.sens:lwr_vm.fbeta) %>%
    mutate(across(everything(), function(x) gsub(" ", "", format(round(x * 100, 1), nsmall = 1)))) %>%
    as.data.frame()
vm.upr <- x %>%
    filter(beta == 1) %>%
    select(upr_vm.sens:upr_vm.fbeta) %>%
    mutate(across(everything(), function(x) gsub(" ", "", format(round(x * 100, 1), nsmall = 1)))) %>%
    as.data.frame()

for (i in 1:nrow(vm)) {
    for (j in 1:ncol(vm)) {
        vm[i, j] <- paste0(vm[i, j], " (", vm.lwr[i, j], "-", vm.upr[i, j], ")")
    }
}

pretty <- cbind(
    x[, c("measure", "median_oc")],
    x[, c("estimate_qm.tp", "estimate_qm.fn", "estimate_qm.tn", "estimate_qm.fp")],
    qm,
    x[, c("estimate_vm.tp", "estimate_vm.fn", "estimate_vm.tn", "estimate_vm.fp")],
    vm
)
write.csv(pretty, file.path(out_path, "tables", "diagnostic_performance_pretty.csv"))

# plot diagnostic performance for differing levels of beta values
plt <- beta_res %>%
    dplyr::select(measure, beta, mean_oc, median_oc, contains("estimate")) %>%
    pivot_longer(cols = contains("estimate")) %>%
    filter(!grepl("_vm", name)) %>%
    mutate(name = gsub("estimate_qm.", "", name)) %>%
    filter(!name %in% c("tp", "fp", "tn", "fn"))

sig <- rbind(
    beta_res %>% select(measure, beta, estimate_qm.sens, sens_padj) %>% rename(value = estimate_qm.sens, padj = sens_padj) %>% mutate(name = "sens"),
    beta_res %>% select(measure, beta, estimate_qm.spec, spec_padj) %>% rename(value = estimate_qm.spec, padj = spec_padj) %>% mutate(name = "spec")
) %>%
    mutate(sig = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*"))

vm <- beta_res %>%
    dplyr::select(measure, beta, mean_oc, median_oc, contains("estimate")) %>%
    pivot_longer(cols = contains("estimate")) %>%
    filter(grepl("_vm", name)) %>%
    mutate(name = gsub("estimate_vm.", "", name)) %>%
    filter(!name %in% c("tp", "fp", "tn", "fn")) %>%
    distinct(measure, beta, name, value) %>%
    filter(name != "fbeta")

plt.ci <- beta_res %>%
    dplyr::select(measure, beta, colnames(.)[grepl("lwr|upr", colnames(.))]) %>%
    pivot_longer(cols = colnames(.)[grepl("lwr|upr", colnames(.))]) %>%
    filter(!grepl("_vm", name)) %>%
    mutate(ci_type = str_split_i(name, "_", 1)) %>%
    mutate(name = gsub("lwr_qm.|upr_qm.", "", name)) %>%
    pivot_wider(names_from = ci_type, values_from = value) %>%
    filter(!name %in% c("tp", "fp", "tn", "fn"))

g <- ggpubr::ggarrange(
    ggplot(plt, aes(x = log2(signif(beta^2, 2)), y = value)) +
        # geom_bar(stat = "identity") +
        geom_ribbon(data = plt.ci, aes(y = NULL, ymin = lwr, ymax = upr), fill = "black", alpha = .2) +
        geom_line(aes(group = measure)) +
        geom_hline(data = vm, aes(yintercept = value), color = "red", linetype = "dashed") +
        geom_text(data = sig, aes(y = value + .2, label = sig), angle = 90, nudge_x = .1) +
        geom_point() +
        facet_grid(
            rows = vars(name),
            cols = vars(measure),
            scales = "free_x",
            labeller = label_bquote(cols = .(ylabs_short[match(levels(factor(plt$measure)), names(ylabs_short))]), rows = .(metric_translator[match(levels(factor(plt$name)), names(metric_translator))]))
        ) +
        labs(x = TeX("$log_2(\\beta^2)"), y = "") +
        theme_bw() +
        scale_y_continuous(breaks = c(0, 0.5, 1)),
    ggplot(plt, aes(x = log2(signif(beta^2, 2)), y = median_oc)) +
        geom_line(aes(group = measure)) +
        geom_point() +
        facet_wrap(~measure, nrow = 1, scales = "free_y", labeller = label_bquote(cols = .(ylabs_short[match(levels(factor(plt$measure)), names(ylabs_short))]))) +
        labs(x = TeX("$log_2(\\beta^2)"), y = "Median cut-off") +
        theme_bw(),
    nrow = 2,
    heights = c(.8, .2),
    labels = "auto"
)

ggsave(file.path(out_path, "diagnostic_performance_bybeta.png"), g, width = 14.4, height = 10, units = "in", dpi = 600)
