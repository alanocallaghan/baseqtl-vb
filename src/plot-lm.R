library("ggplot2")
library("ggdist")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
library("geomtextpath")
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
fpath <- "fig"

################################################################################
##
## Read/merge
##
################################################################################

# "optimizing",
methods <- c("sampling", "vb")
dfs <- lapply(methods,
    function(method) {
        cat(method, "\n")
        mtol <- if (method == "vb") sprintf("vb_%1.0e", tol) else method
        combfile <- sprintf("rds/%s/%s_combined.rds", model, mtol)
        if (file.exists(combfile)) {
            return(readRDS(combfile))
        }
    }
)


dfs[] <- lapply(dfs,
    function(df) {
        ## TRUE means that the HPD interval is one side of zero (sig)
        ## FALSE means it overlaps zero (null)
        df$null.99 <- sign(df$"0.5%") == sign(df$"99.5%")
        df$null.95 <- sign(df$"2.5%") == sign(df$"97.5%")
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
## Discrepancy
##
################################################################################

cmdf <- mdf
cmdf <- cmdf %>% mutate(discrepancy = mean.hmc - mean.vb)

x <- cmdf %>%
    arrange(-abs(discrepancy)) %>%
    top_n(50, abs(discrepancy)) %>%
    select(
        vb = mean.vb,
        vb_low = `5.0%.vb`,
        vb_high = `95.0%.vb`,
        hmc = mean.hmc,
        hmc_low = `5.0%.hmc`,
        hmc_high = `95.0%.hmc`,
        discrepancy = discrepancy,
        snp,
        gene
    )

saveRDS(x,
    sprintf("rds/%s_discrepancies_vb_%1.0e.rds", model, tol)
)

cmdf <- cmdf %>%
    mutate(
        disc_sc_sdh = discrepancy / sd.hmc,
        disc_sc_sdv = discrepancy / sd.vb,
        disc_sc_meanh = discrepancy / mean.hmc,
        disc_sc_meanv = discrepancy / mean.vb
    )

mname <- "ADVI"
method <- "vb"
r <- range(c(cmdf$mean.hmc, cmdf$mean.vb))
## not se mean but other hmc se
diag_vars <- c(
    "gene",
    "Rhat",
    "khat",
    "converged",
    "niter",
    "n_eff.hmc", "time.hmc", "n_ase",
    "se_mean.hmc", "sd.hmc",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)
cmdf$converged <- factor(cmdf$converged)
for (type in c("discrepancy")) {
# for (type in c("discrepancy", "disc_sc_sdh", "disc_sc_sdh", "disc_sc_meanh", "disc_sc_meanv")) {
    for (x in diag_vars) {    
        geom <- if (is.numeric(cmdf[[x]])) {
            geom_point(size = 0.8) 
        } else if (is.logical(cmdf[[x]]) || is.character(cmdf[[x]]) || is.factor(cmdf[[x]])) {
            geom_boxplot(fill = "grey80")
        }
        g <- ggplot(cmdf) +
            aes_string(x, sprintf("abs(%s)", type)) +
            geom +
            geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
            # scale_y_log10() +
            labs(x = x, y = "Discrepancy")

        ggsave(
            sprintf("%s/%s/diag/%s_%s_%s.png", fpath, model, type, gsub("\\.", "_", x), method),
            width = 7, height = 7
        )
    }
}

g <- ggplot(cmdf) +
    aes(n_eff.hmc, discrepancy) +
    geom_point(size = 0.7, alpha = 0.7) +
    geom_vline(xintercept = 500, linetype = "dashed") +
    labs(x = "Effective sample size", y = "Discrepancy (ADVI - HMC)")
ggsave(
    sprintf("%s/%s/diag/disc_ESS_vb.png", fpath, model),
    width = 5, height = 5
)
g <- ggplot(cmdf) +
    aes(niter, discrepancy) +
    geom_point(size = 0.7, alpha = 0.7) +
    labs(x = "Number of iterations before convergence", y = "Discrepancy (ADVI - HMC)")
ggsave(
    sprintf("%s/%s/diag/disc_niter_vb.png", fpath, model),
    width = 5, height = 5
)

g <- ggplot(cmdf) +
    aes(khat, discrepancy) +
    geom_point(size = 0.7, alpha = 0.7) +
    geom_vline(xintercept = 0.7, linetype = "dashed") +
    labs(x = expression(hat(K)), y = "Discrepancy (ADVI - HMC)")
ggsave(
    sprintf("%s/%s/diag/disc_khat_vb.png", fpath, model),
    width = 5, height = 5
)

g <- ggplot(cmdf) +
    aes(time.hmc, discrepancy) +
    geom_point() +
    labs(x = "Time taken for HMC (s)", y = "Discrepancy (ADVI - HMC)")
ggsave(
    sprintf("%s/%s/diag/time_vb.png", fpath, model),
    width = 7, height = 7
)



################################################################################
##
## Time
##
################################################################################

cmdf %>%
    group_by(gene) %>%
    summarise(time = sum(time.hmc), nsnps = n()) %>%
    ggplot() +
    aes(nsnps, time) +
    labs(x = "Number of SNPs in cis window", y = "Time (s)") +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) +
    scale_x_log10() +
    scale_y_log10() -> g
ggsave(
    sprintf("%s/%s/time/nsnp.png", fpath, model),
    width = 5, height = 5
)



g <- ggplot(cmdf) +
    aes(time.hmc, time.vb) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_pointdensity(size = 0.8) +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_viridis() +
    theme(legend.position = "none") +
    labs(x = "Time (s) for HMC", y = "Time (s) for ADVI")
ggsave(
    sprintf("%s/%s/time/time_hmc_vs_vb.png", fpath, model),
    width = 5, height = 5
)

cn <- lapply(dfs, colnames)
cn <- Reduce(intersect, cn)
dfs_int <- lapply(dfs, function(x) x[, cn])
dfs_int[] <- lapply(
    names(dfs_int),
    function(n) {
        dfs_int[[n]]$method <- n
        dfs_int[[n]]
    }
)
df_int <- do.call(rbind, dfs_int)

g <- ggplot(df_int) +
    aes(x = time, colour = method) +
    geom_density() +
    geom_vline(
        aes(
            xintercept = mean(time[method == "vb"]),
            colour = "vb"
        ),
        linetype = "dashed"
    ) +
    geom_vline(
        aes(
            xintercept = mean(time[method == "sampling"]),
            colour = "sampling"
        ),
        linetype = "dashed"
    ) +
    scale_x_log10(name = "Time (s)") +
    scale_colour_brewer(palette = "Set1", name = "Method") +
    ylab("Density") +
    theme(legend.position = "bottom")
ggsave(
    sprintf("%s/%s/time/time_comparison.png", fpath, model),
    width = 5, height = 5
)





################################################################################
##
## Point estimates
##
################################################################################


cmdf <- cmdf[cmdf$n_eff.hmc > minEff & cmdf$Rhat < maxRhat, ]
## from PSIS paper, arxiv 1507.02646
# cmdf <- cmdf[cmdf$khat < 0.7, ]

gdf <- group_by(cmdf, gene)

lab_str <- "HMC: %s\nADVI: %s\n"
cmdf$nullstr95 <- sprintf(lab_str,
    ifelse(cmdf$null.95.hmc, "yes", "no"),
    ifelse(cmdf$null.95.vb, "yes", "no")
)
cmdf$nullstr99 <- sprintf(lab_str,
    ifelse(cmdf$null.99.hmc, "yes", "no"),
    ifelse(cmdf$null.99.vb, "yes", "no")
)

null_levs <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "no",  "yes", "yes")
)
null_ord <- null_levs[c(1, 3, 2, 4)]
scale <- scale_colour_manual(
    values = setNames(c("#fb9a99", "#e31a1c", "#a6cee3", "#1f78b4"), null_levs),
    name = "Significance"
)
cmdf$nullstr95 <- factor(cmdf$nullstr95, levels = null_ord)
cmdf$nullstr99 <- factor(cmdf$nullstr99, levels = null_ord)



mdfs <- cmdf
mdfs <- cmdf[cmdf$mean.vb < 2, ]
lim <- range(c(mdfs$mean.hmc, mdfs$mean.vb))

g <- ggplot(mdfs[order(mdfs$nullstr95), ]) +
    aes(mean.hmc, mean.vb, colour = nullstr95) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "HMC estimate", y = "ADVI estimate") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = "bottom")
ggsave(
    sprintf("%s/%s/estimates/point-estimates-95.png", fpath, model),
    width = 5, height = 5
)

g <- ggplot(mdfs[order(mdfs$nullstr99), ]) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "HMC estimate", y = "ADVI estimate") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = "bottom")
ggsave(
    sprintf("%s/%s/estimates/point-estimates-99.png", fpath, model),
    width = 5, height = 5
)

mdfss <- mdfs[!(mdfs$null.95.hmc == mdfs$null.95.vb), ]

g <- ggplot(mdfss) +
    aes(mean.hmc, mean.vb, colour = nullstr95) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "HMC estimate", y = "ADVI estimate") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = "bottom")
ggsave(
    sprintf("%s/%s/estimates/point-estimates-diff-95.png", fpath, model),
    width = 5, height = 5
)
mdfss <- mdfs[!(mdfs$null.99.hmc == mdfs$null.99.vb), ]

g <- ggplot(mdfss) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "HMC estimate", y = "ADVI estimate") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = "bottom")
ggsave(
    sprintf("%s/%s/estimates/point-estimates-diff-99.png", fpath, model),
    width = 5, height = 5
)




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
    sprintf("%s/%s/estimates/lm-hmc-95.png", fpath, model),
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
    sprintf("%s/%s/estimates/glm-hmc-95.png", fpath, model),
    width = 5, height = 5
)


lim <- range(lm_df$coef)
g <- ggplot() +
    aes(
        lm_df[lm_df$method=="lm", "coef"],
        lm_df[lm_df$method=="glm", "coef"]
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_pointdensity(size = 0.5, alpha = 0.7) +
    scale_colour_viridis() +
    lims(x = lim, y = lim) +
    labs(x = "lm estimate", y = "glm estimate")
    #  + scale
ggsave(
    sprintf("%s/%s/estimates/lm-glm-coef.png", fpath, model),
    width = 5, height = 5
)

lim <- range(lm_df$padj)
g <- ggplot() +
    aes(
        lm_df[lm_df$method=="lm", "padj"],
        lm_df[lm_df$method=="glm", "padj"]
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
    sprintf("%s/%s/estimates/lm-glm-pval.png", fpath, model),
    width = 5, height = 5
)


################################################################################
##
## AUPR/AUROC
##
################################################################################


lm_dfs <- lapply(c("lm", "glm"),
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

    g <- ggplot() +
        geom_path(aes(sens_lm, time_lm, colour = "lm")) +
        geom_path(aes(sens_glm, time_glm, colour = "glm")) +
        geom_path(aes(sens_vb, time_vb, colour = "ADVI")) +
        geom_texthline(
            yintercept = sum(df$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        # geom_text_repel(
        #     aes(label = format(thresholds, digits=2)),
        #     vjust = 1, hjust = 1,
        #     segment.color = "grey70"
        # ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(df$time.hmc), time_vb, time_lm, time_glm)) +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom")
    ggsave(
        sprintf("%s/%s/time/time_vs_sens_all_%s.png", fpath, model, t),
        width = 5, height = 5
    )

    g <- ggplot() +
        geom_path(aes(spec_lm, time_lm, colour = "lm")) +
        geom_path(aes(spec_glm, time_glm, colour = "glm")) +
        geom_path(aes(spec_vb, time_vb, colour = "ADVI")) +
        geom_texthline(
            yintercept = sum(df$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        # geom_text_repel(
        #     aes(label = thresholds), vjust = 1, hjust = 1,
        #     segment.color = "grey70"
        # ) +
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(df$time.hmc), time_vb, time_lm, time_glm)) +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom")
    ggsave(
        sprintf("%s/%s/time/time_vs_spec_%s.png", fpath, model, t),
        width = 5, height = 5
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
                colour = sprintf("lm; AUPRC: %0.3f", perf_aupr_lm)
            )
        ) +
        geom_path(
            aes(
                perf_pr_glm@x.values[[1]], perf_pr_glm@y.values[[1]],
                colour = sprintf("glm; AUPRC: %0.3f", perf_aupr_glm)
            )
        ) +
        geom_path(
            aes(
                perf_pr_vb@x.values[[1]], perf_pr_vb@y.values[[1]],
                colour = sprintf("ADVI; AUPRC: %0.3f", perf_aupr_vb)
            )
        ) +
        lims(x = 0:1, y = 0:1) +
        labs(x = "Precision", y = "Recall") +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom")
    ggsave(
        sprintf("%s/%s/roc/pr_all_%s.png", fpath, model, t),
        width = 5, height = 5
    )

    g <- ggplot() +
        geom_path(
            aes(
                perf_roc_lm@x.values[[1]], perf_roc_lm@y.values[[1]],
                colour = sprintf("lm; AUROC: %0.3f", perf_auroc_lm)
            )
        ) +
        geom_path(
            aes(
                perf_roc_glm@x.values[[1]], perf_roc_glm@y.values[[1]],
                colour = sprintf("glm; AUROC: %0.3f", perf_auroc_glm)
            )
        ) +
        geom_path(
            aes(
                perf_roc_vb@x.values[[1]], perf_roc_vb@y.values[[1]],
                colour = sprintf("ADVI; AUROC: %0.3f", perf_auroc_vb)
            )
        ) +
        lims(x = 0:1, y = 0:1) +
        labs(x = "Sensitivity", y = "Specificity") +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom")
    ggsave(
        sprintf("%s/%s/roc/roc_all_%s.png", fpath, model, t),
        width = 5, height = 5
    )
}
