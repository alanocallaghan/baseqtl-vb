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
source("src/functions.R")
mkfigdir(fpath, model)

# "optimizing",
methods <- c("vb", "sampling")
dfs <- lapply(
    methods,
    function(method) {
        cat(method, "\n")
        mtol <- mtol(method, tol)
        combfile <- sprintf("rds/%s/%s_combined.rds", model, mtol)
        print(combfile)
        if (file.exists(combfile)) {
            return(readRDS(combfile))
        }
    }
)
dfs[] <- lapply(dfs, add_nulls)
dfs[] <- lapply(dfs, function(df) {
    df[df$seed == 42, ]
})

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
df$time.hmc <- df$time.hmc / 3600
df$time.vb <- df$time.vb / 3600
df$time.glm <- df$time.glm / 3600
df$time.lm <- df$time.lm / 3600


get_sens_spec <- function(
        truth,
        prob,
        thresholds,
        which,
        comparison = `<`
    ) {
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
get_time <- function(
        times_hmc,
        times_approx,
        prob,
        thresholds,
        comparison = `<`
    ) {
    comparison <- match.fun(comparison)
    sapply(
        thresholds,
        function(threshold) {
            (sum(times_hmc[comparison(prob, threshold)]) + sum(times_approx))
        }
    )
}

p_threshold <- 99
p_column <- paste0("null.", p_threshold, ".hmc")
thresholds_lm <- 10^seq(0, -5, length.out = 200)
thresholds_glm <- 10^seq(0, -300, length.out = 200)
thresholds_vb <- seq(min(df$PEP.vb), max(df$PEP.vb), length.out = 200)

sens_lm <- get_sens_spec(
    df[[p_column]],
    df[["padj.lm"]],
    thresholds_lm,
    which = "sens",
    comparison = `<`
)
sens_glm <- get_sens_spec(
    df[[p_column]],
    df[["padj.glm"]],
    thresholds_glm,
    which = "sens",
    comparison = `<`
)
sens_vb <- get_sens_spec(
    df[[p_column]],
    df[["PEP.vb"]],
    thresholds_vb,
    which = "sens",
    comparison = `>`
)

spec_lm <- get_sens_spec(
    df[[p_column]],
    df[["padj.lm"]],
    thresholds_lm,
    which = "spec",
    comparison = `<`
)
spec_glm <- get_sens_spec(
    df[[p_column]],
    df[["padj.glm"]],
    thresholds_glm,
    which = "spec",
    comparison = `<`
)
spec_vb <- get_sens_spec(
    df[[p_column]],
    df[["PEP.vb"]],
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
    df[["PEP.vb"]],
    thresholds_vb,
    comparison = `>`
)

gtime <- ggplot() +
    geom_path(aes(sens_lm, time_lm, colour = "lm"), na.rm = TRUE) +
    geom_path(aes(sens_glm, time_glm, colour = "glm"), na.rm = TRUE) +
    geom_path(aes(sens_vb, time_vb, colour = "ADVI"), na.rm = TRUE) +
    geom_hline(
        yintercept = sum(df$time.hmc),
        linetype = "dashed"
    ) +
    annotate(
        geom = "text",
        x = 0.5,
        y = sum(df$time.hmc),
        label = "Total time\nwithout screening",
        vjust = -0.3,
    ) +
    labs(x = "Sensitivity", y = "Total time (hr)") +
    # ylim(0, 2700) +
    ylim(0, max(sum(df$time.hmc), time_vb, time_lm, time_glm) * 1.1) +
    scale_colour_brewer(palette = "Set2", name = NULL) +
    theme(legend.position = "bottom")
# ggsave(
#     sprintf("%s/%s/time/time_vs_sens_all_%s.pdf", fpath, model, p_threshold),
#     width = 4, height = 4.5
# )


pred_lm <- prediction(
    predictions = 1 - df[["padj.lm"]],
    labels = factor(df[[p_column]], levels = c(TRUE, FALSE))
)
pred_glm <- prediction(
    predictions = 1 - df[["padj.glm"]],
    labels = factor(df[[p_column]], levels = c(TRUE, FALSE))
)
pred_vb <- prediction(
    predictions = df[["PEP.vb"]],
    labels = factor(df[[p_column]], levels = c(TRUE, FALSE))
)

perf_auroc_lm <- performance(pred_lm, "auc")@y.values[[1]]
perf_auroc_glm <- performance(pred_glm, "auc")@y.values[[1]]
perf_auroc_vb <- performance(pred_vb, "auc")@y.values[[1]]

perf_roc_lm <- performance(pred_lm, "tpr", "fpr")
perf_roc_glm <- performance(pred_glm, "tpr", "fpr")
perf_roc_vb <- performance(pred_vb, "tpr", "fpr")



groc <- ggplot() +
    geom_path(
        aes(
            perf_roc_lm@x.values[[1]], perf_roc_lm@y.values[[1]],
            colour = sprintf("lm;\nAUROC: %0.3f", perf_auroc_lm)
        ), na.rm = TRUE
    ) +
    geom_path(
        aes(
            perf_roc_glm@x.values[[1]], perf_roc_glm@y.values[[1]],
            colour = sprintf("glm;\nAUROC: %0.3f", perf_auroc_glm)
        ), na.rm = TRUE
    ) +
    geom_path(
        aes(
            perf_roc_vb@x.values[[1]], perf_roc_vb@y.values[[1]],
            colour = sprintf("ADVI;\nAUROC: %0.3f", perf_auroc_vb)
        ), na.rm = TRUE
    ) +
    lims(x = 0:1, y = 0:1) +
    labs(x = "Sensitivity", y = "Specificity") +
    scale_colour_brewer(palette = "Set2", name = NULL) +
    theme(legend.position = "bottom")

plot_with_legend_below(groc, gtime, labels = "AUTO")
ggsave(
    sprintf("%s/%s/roc/time_roc_all_%s.pdf", fpath, model, p_threshold),
    width = 5.5, height = 3.5
)




## specificity

# g <- ggplot() +
#     geom_path(aes(spec_lm, time_lm, colour = "lm"), na.rm = TRUE) +
#     geom_path(aes(spec_glm, time_glm, colour = "glm"), na.rm = TRUE) +
#     geom_path(aes(spec_vb, time_vb, colour = "ADVI"), na.rm = TRUE) +
#     geom_hline(
#         yintercept = sum(df$time.hmc),
#         linetype = "dashed"
#     ) +
#     annotate(
#         geom = "text",
#         x = 0.5,
#         y = sum(df$time.hmc),
#         label = "Total time\nwithout screening",
#         vjust = -0.3,
#     ) +
#     labs(x = "Specificity", y = "Total time (hr)") +
#     scale_y_log10(
#         limits = c(1, max(sum(df$time.hmc), time_vb, time_lm, time_glm))
#     ) +
#     scale_colour_brewer(palette = "Set2", name = NULL) +
#     theme(legend.position = "bottom")
# ggsave(
#     sprintf("%s/%s/time/time_vs_spec_%s.pdf", fpath, model, p_threshold),
#     width = 4, height = 4
# )


##### plots of estimates and coefs/pvals... not needed but useful maybe



# mdf_lm <- merge(mdf, lm_df[lm_df$method == "lm", ],
#     by = c("snp", "gene"),
#     suffixes = c(".hmc", ".lm")
# )
# mdf_glm <- merge(mdf, lm_df[lm_df$method == "glm", ],
#     by = c("snp", "gene"),
#     suffixes = c(".hmc", ".glm")
# )




# lim <- range(c(mdf_lm$mean.hmc, mdf_lm$coef))
# g <- ggplot(mdf_lm) +
#     aes(mean.hmc, coef, colour = -log10(pval)) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_point(size = 0.5, alpha = 0.7) +
#     scale_colour_viridis() +
#     lims(x = lim, y = lim) +
#     labs(x = "BaseQTL estimate (HMC)", y = "lm estimate")
# ggsave(
#     sprintf("%s/%s/estimates/lm-hmc-95.pdf", fpath, model),
#     width = 5, height = 5
# )

# lim <- range(c(mdf_glm$mean.hmc, mdf_glm$coef))
# g <- ggplot(mdf_glm) +
#     aes(mean.hmc, coef, colour = -log10(pval)) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_point(size = 0.5, alpha = 0.7) +
#     scale_colour_viridis() +
#     lims(x = lim, y = lim) +
#     labs(x = "BaseQTL estimate (HMC)", y = "glm estimate")
# ggsave(
#     sprintf("%s/%s/estimates/glm-hmc-95.pdf", fpath, model),
#     width = 5, height = 5
# )

# lim <- range(lm_df$coef)
# g <- ggplot() +
#     aes(
#         lm_df[lm_df$method == "lm", "coef"],
#         lm_df[lm_df$method == "glm", "coef"]
#     ) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_pointdensity(size = 0.5, alpha = 0.7) +
#     scale_colour_viridis() +
#     lims(x = lim, y = lim) +
#     labs(x = "lm estimate", y = "glm estimate")
# ggsave(
#     sprintf("%s/%s/estimates/lm-glm-coef.pdf", fpath, model),
#     width = 5, height = 5
# )

# lim <- range(lm_df$padj)
# g <- ggplot() +
#     aes(
#         lm_df[lm_df$method == "lm", "padj"],
#         lm_df[lm_df$method == "glm", "padj"]
#     ) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_pointdensity(size = 0.5, alpha = 0.7) +
#     scale_colour_viridis() +
#     scale_x_log10() +
#     scale_y_log10() +
#     # lims(x = lim, y = lim) +
#     labs(x = "lm pval (FDR)", y = "glm pval (FDR)")
# ggsave(
#     sprintf("%s/%s/estimates/lm-glm-pval.pdf", fpath, model),
#     width = 5, height = 5
# )



##### AUPR plots... also not necessary but maybe nice

# perf_aupr_lm <- performance(pred_lm, "aucpr")@y.values[[1]]
# perf_aupr_glm <- performance(pred_glm, "aucpr")@y.values[[1]]
# perf_aupr_vb <- performance(pred_vb, "aucpr")@y.values[[1]]

# perf_pr_lm <- performance(pred_lm, "prec", "rec")
# perf_pr_glm <- performance(pred_glm, "prec", "rec")
# perf_pr_vb <- performance(pred_vb, "prec", "rec")

# g <- ggplot() +
#     geom_path(
#         aes(
#             perf_pr_lm@x.values[[1]], perf_pr_lm@y.values[[1]],
#             colour = sprintf("lm;\nAUPRC: %0.3f", perf_aupr_lm)
#         )
#     ) +
#     geom_path(
#         aes(
#             perf_pr_glm@x.values[[1]], perf_pr_glm@y.values[[1]],
#             colour = sprintf("glm;\nAUPRC: %0.3f", perf_aupr_glm)
#         )
#     ) +
#     geom_path(
#         aes(
#             perf_pr_vb@x.values[[1]], perf_pr_vb@y.values[[1]],
#             colour = sprintf("ADVI;\nAUPRC: %0.3f", perf_aupr_vb)
#         )
#     ) +
#     lims(x = 0:1, y = 0:1) +
#     labs(x = "Precision", y = "Recall") +
#     scale_colour_brewer(palette = "Set2", name = NULL) +
#     theme(legend.position = "bottom")
# ggsave(
#     sprintf("%s/%s/roc/pr_all_%s.pdf", fpath, model, p_threshold),
#     width = 4, height = 4
# )
