library("ggplot2")
library("ggdist")
library("ggpointdensity")
library("cowplot")
library("baseqtl")
library("viridis")
library("yardstick")
library("dplyr")
library("ggrepel")
library("argparse")
library("ROCR")
library("pROC")

theme_set(theme_bw())

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "GT",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 0.01,
    type = "double"
)

maxRhat <- 1.1 ## from baseqtl-paper repo
minEff <- 500 ## from stan docs (-ish)

args <- parser$parse_args()
model <- args[["model"]]
tol <- args[["tolerance"]]

lm_df <- readRDS(sprintf("rds/%s/lm-filtering.rds", model))
fpath <- sprintf("fig_%1.0e", tol)
# fpath <- "fig"
source("src/functions.R")
mkfigdir(fpath, model)

# "optimizing",
methods <- c("vb", "sampling")
dfs <- lapply(
    methods,
    function(method) {
        cat(method, "\n")
        mtol <- if (method == "vb") sprintf("vb_%1.0e", tol) else method
        combfile <- sprintf("rds/%s/%s_combined.rds", model, mtol)
        if (file.exists(combfile)) {
            return(readRDS(combfile))
        }
    }
)

dfs[] <- lapply(
    dfs,
    function(df) {
        ## TRUE means that the HPD interval is one side of zero (sig)
        ## FALSE means it overlaps zero (null)
        df$null.95 <- sign(df$"2.5%") == sign(df$"97.5%")
        df$null.90 <- sign(df$"5.0%") == sign(df$"95.0%")
        df$null.80 <- sign(df$"10.0%") == sign(df$"90.0%")
        df$null.70 <- sign(df$"15.0%") == sign(df$"85.0%")
        df$null.60 <- sign(df$"20.0%") == sign(df$"80.0%")
        df$null.50 <- sign(df$"25.0%") == sign(df$"75.0%")
        df$null.40 <- sign(df$"30.0%") == sign(df$"70.0%")
        df$null.30 <- sign(df$"35.0%") == sign(df$"65.0%")
        df$null.20 <- sign(df$"40.0%") == sign(df$"60.0%")
        df$null.10 <- sign(df$"45.0%") == sign(df$"55.0%")
        df
    }
)

names(dfs) <- methods
by <- if (model == "GT") {
    c(
        "test", "gene", "snp", "n_tot", "n_ase",
        "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
    )
} else {
    c(
        "test", "gene", "snp", "condition",
        "n_tot", "n_ase", "mean_ase", "sd_ase", "mean_count", "sd_count",
        "n_wt", "n_het", "n_hom"
    )
}
mdf <- merge(dfs[["vb"]], dfs[["sampling"]], by = by, suffix = c(".vb", ".hmc"))

################################################################################
##
## lm/glm point estimates
##
################################################################################
mdf <- mdf[mdf$n_eff.hmc > minEff & mdf$Rhat < maxRhat, ]
## from PSIS paper, arxiv 1507.02646
# mdf <- mdf[mdf$khat < 0.7, ]

rownames(lm_df) <- NULL
mdf_lm <- merge(mdf, lm_df[lm_df$method == "lm", ],
    by = c("snp", "gene"),
    suffixes = c(".hmc", ".lm")
)
mdf_glm <- merge(mdf, lm_df[lm_df$method == "glm", ],
    by = c("snp", "gene"),
    suffixes = c(".hmc", ".glm")
)


##### plots of estimates and coefs/pvals... not needed but useful maybe

lim <- range(c(mdf_lm$mean.hmc, mdf_lm$coef))
g <- ggplot(mdf_lm) +
    aes(mean.hmc, coef, colour = -log10(pval)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_viridis() +
    lims(x = lim, y = lim) +
    labs(x = "BaseQTL estimate (HMC)", y = "lm estimate")
#  + scale
ggsave(
    sprintf("%s/%s/estimates/lm-hmc-95.pdf", fpath, model),
    width = 5, height = 5
)

lim <- range(c(mdf_glm$mean.hmc, mdf_glm$coef))
g <- ggplot(mdf_glm) +
    aes(mean.hmc, coef, colour = -log10(pval)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_viridis() +
    lims(x = lim, y = lim) +
    labs(x = "BaseQTL estimate (HMC)", y = "glm estimate")
#  + scale
ggsave(
    sprintf("%s/%s/estimates/glm-hmc-95.pdf", fpath, model),
    width = 5, height = 5
)

lim <- range(lm_df$coef)
g <- ggplot() +
    aes(
        lm_df[lm_df$method == "lm", "coef"],
        lm_df[lm_df$method == "glm", "coef"]
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_pointdensity(size = 0.5, alpha = 0.7) +
    scale_colour_viridis() +
    lims(x = lim, y = lim) +
    labs(x = "lm estimate", y = "glm estimate")
#  + scale
ggsave(
    sprintf("%s/%s/estimates/lm-glm-coef.pdf", fpath, model),
    width = 5, height = 5
)

lim <- range(lm_df$padj)
g <- ggplot() +
    aes(
        lm_df[lm_df$method == "lm", "padj"],
        lm_df[lm_df$method == "glm", "padj"]
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_pointdensity(size = 0.5, alpha = 0.7) +
    scale_colour_viridis() +
    scale_x_log10() +
    scale_y_log10() +
    # lims(x = lim, y = lim) +
    labs(x = "lm pval (FDR)", y = "glm pval (FDR)")
#  + scale
ggsave(
    sprintf("%s/%s/estimates/lm-glm-pval.pdf", fpath, model),
    width = 5, height = 5
)

################################################################################
##
## AUPR/AUROC
##
################################################################################

lm_dfs <- lapply(
    c("lm", "glm"),
    function(lmod) {
        lm_df[lm_df$method == lmod, ]
    }
)
mlm_df <- merge(
    lm_dfs[[1]], lm_dfs[[2]],
    by = c("gene", "snp"),
    suffixes = c(".lm", ".glm")
)

df <- merge(
    mdf, mlm_df,
    by = c("snp", "gene"),
    suffixes = c(".hmc", ".lm")
)

get_sens_spec <- function(truth,
                                                    prob,
                                                    thresholds,
                                                    which,
                                                    comparison = `<`) {
    fun <- match.fun(paste(which, "vec", sep = "_"))
    comparison <- match.fun(comparison)
    sapply(
        thresholds,
        function(threshold) {
            fun(
                truth = factor(truth, levels = c(TRUE, FALSE)),
                estimate = factor(
                    comparison(prob, threshold),
                    levels = c(TRUE, FALSE)
                )
            )
        }
    )
}
get_time <- function(times_hmc,
                                         times_approx,
                                         prob,
                                         thresholds,
                                         comparison = `<`) {
    comparison <- match.fun(comparison)
    sapply(
        thresholds,
        function(threshold) {
            sum(times_hmc[comparison(prob, threshold)]) + sum(times_approx)
        }
    )
}

## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t, ".hmc")
    thresholds_lm <- 10^seq(0, -5, length.out = 200)
    thresholds_glm <- 10^seq(0, -300, length.out = 200)
    thresholds_vb <- seq(min(df$prob.vb), max(df$prob.vb), length.out = 200)

    sens_lm <- get_sens_spec(
        df[[bc]],
        df[["padj.lm"]],
        thresholds_lm,
        which = "sens",
        comparison = `<`
    )
    sens_glm <- get_sens_spec(
        df[[bc]],
        df[["padj.glm"]],
        thresholds_glm,
        which = "sens",
        comparison = `<`
    )
    sens_vb <- get_sens_spec(
        df[[bc]],
        df[["prob.vb"]],
        thresholds_vb,
        which = "sens",
        comparison = `>`
    )

    spec_lm <- get_sens_spec(
        df[[bc]],
        df[["padj.lm"]],
        thresholds_lm,
        which = "spec",
        comparison = `<`
    )
    spec_glm <- get_sens_spec(
        df[[bc]],
        df[["padj.glm"]],
        thresholds_glm,
        which = "spec",
        comparison = `<`
    )
    spec_vb <- get_sens_spec(
        df[[bc]],
        df[["prob.vb"]],
        thresholds_vb,
        which = "spec",
        comparison = `>`
    )

    time_lm <- get_time(
        df$time.hmc,
        df$time.lm,
        df[["padj.lm"]],
        thresholds_lm
    )
    time_glm <- get_time(
        df$time.hmc,
        df$time.glm,
        df[["padj.glm"]],
        thresholds_glm
    )
    time_vb <- get_time(
        df$time.hmc,
        df$time.vb,
        df[["prob.vb"]],
        thresholds_vb,
        comparison = `>`
    )
    # if (t == 99) stop()
    min(time_vb[sens_vb == 1] / sum(df$time.hmc))

    gtime <- ggplot() +
        geom_path(aes(sens_lm, time_lm, colour = "lm")) +
        geom_path(aes(sens_glm, time_glm, colour = "glm")) +
        geom_path(aes(sens_vb, time_vb, colour = "ADVI")) +
        geom_hline(
            yintercept = sum(df$time.hmc),
            linetype = "dashed"
        ) +
        annotate(
            geom = "text",
            x = 0.5,
            y = sum(df$time.hmc),
            label = "Total time without screening",
            vjust = -0.3,
        ) +
        # geom_text_repel(
        #     aes(label = format(thresholds, digits=2)),
        #     vjust = 1, hjust = 1,
        #     segment.color = "grey70"
        # ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        # scale_y_log10(
        #     limits = c(1, max(sum(df$time.hmc), time_vb, time_lm, time_glm))
        # ) +
        ylim(0, max(sum(df$time.hmc), time_vb, time_lm, time_glm)) +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom")
    ggsave(
        sprintf("%s/%s/time/time_vs_sens_all_%s.pdf", fpath, model, t),
        width = 4, height = 4.5
    )


    g <- ggplot() +
        geom_path(aes(spec_lm, time_lm, colour = "lm")) +
        geom_path(aes(spec_glm, time_glm, colour = "glm")) +
        geom_path(aes(spec_vb, time_vb, colour = "ADVI")) +
        geom_hline(
            yintercept = sum(df$time.hmc),
            linetype = "dashed"
        ) +
        annotate(
            geom = "text",
            x = 0.5,
            y = sum(df$time.hmc),
            label = "Total time without screening",
            vjust = -0.3,
        ) +
        # geom_text_repel(
        #     aes(label = thresholds), vjust = 1, hjust = 1,
        #     segment.color = "grey70"
        # ) +
        labs(x = "Specificity", y = "Total time (s)") +
        scale_y_log10(
            limits = c(1, max(sum(df$time.hmc), time_vb, time_lm, time_glm))
        ) +
        # ylim(0, max(sum(df$time.hmc), time_vb, time_lm, time_glm)) +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom")
    ggsave(
        sprintf("%s/%s/time/time_vs_spec_%s.pdf", fpath, model, t),
        width = 4, height = 4
    )

    pred_lm <- prediction(
        predictions = 1 - df[["padj.lm"]],
        labels = factor(df[[bc]], levels = c(TRUE, FALSE))
    )
    pred_glm <- prediction(
        predictions = 1 - df[["padj.glm"]],
        labels = factor(df[[bc]], levels = c(TRUE, FALSE))
    )
    pred_vb <- prediction(
        predictions = df[["prob.vb"]],
        labels = factor(df[[bc]], levels = c(TRUE, FALSE))
    )

    perf_auroc_lm <- performance(pred_lm, "auc")@y.values[[1]]
    perf_auroc_glm <- performance(pred_glm, "auc")@y.values[[1]]
    perf_auroc_vb <- performance(pred_vb, "auc")@y.values[[1]]

    perf_aupr_lm <- performance(pred_lm, "aucpr")@y.values[[1]]
    perf_aupr_glm <- performance(pred_glm, "aucpr")@y.values[[1]]
    perf_aupr_vb <- performance(pred_vb, "aucpr")@y.values[[1]]

    perf_pr_lm <- performance(pred_lm, "prec", "rec")
    perf_pr_glm <- performance(pred_glm, "prec", "rec")
    perf_pr_vb <- performance(pred_vb, "prec", "rec")

    perf_roc_lm <- performance(pred_lm, "tpr", "fpr")
    perf_roc_glm <- performance(pred_glm, "tpr", "fpr")
    perf_roc_vb <- performance(pred_vb, "tpr", "fpr")

    g <- ggplot() +
        geom_path(
            aes(
                perf_pr_lm@x.values[[1]], perf_pr_lm@y.values[[1]],
                colour = sprintf("lm;\nAUPRC: %0.3f", perf_aupr_lm)
            )
        ) +
        geom_path(
            aes(
                perf_pr_glm@x.values[[1]], perf_pr_glm@y.values[[1]],
                colour = sprintf("glm;\nAUPRC: %0.3f", perf_aupr_glm)
            )
        ) +
        geom_path(
            aes(
                perf_pr_vb@x.values[[1]], perf_pr_vb@y.values[[1]],
                colour = sprintf("ADVI;\nAUPRC: %0.3f", perf_aupr_vb)
            )
        ) +
        lims(x = 0:1, y = 0:1) +
        labs(x = "Precision", y = "Recall") +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom")
    ggsave(
        sprintf("%s/%s/roc/pr_all_%s.pdf", fpath, model, t),
        width = 4, height = 4
    )

    groc <- ggplot() +
        geom_path(
            aes(
                perf_roc_lm@x.values[[1]], perf_roc_lm@y.values[[1]],
                colour = sprintf("lm;\nAUROC: %0.3f", perf_auroc_lm)
            )
        ) +
        geom_path(
            aes(
                perf_roc_glm@x.values[[1]], perf_roc_glm@y.values[[1]],
                colour = sprintf("glm;\nAUROC: %0.3f", perf_auroc_glm)
            )
        ) +
        geom_path(
            aes(
                perf_roc_vb@x.values[[1]], perf_roc_vb@y.values[[1]],
                colour = sprintf("ADVI;\nAUROC: %0.3f", perf_auroc_vb)
            )
        ) +
        lims(x = 0:1, y = 0:1) +
        labs(x = "Sensitivity", y = "Specificity") +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom")
    ggsave(
        sprintf("%s/%s/roc/roc_all_%s.pdf", fpath, model, t),
        width = 4, height = 4.5
    )

    plot_with_legend_below(groc, gtime, labels = "AUTO")
    ggsave(
        sprintf("%s/%s/roc/time_roc_all_%s.pdf", fpath, model, t),
        width = 8, height = 4.5
    )
}
