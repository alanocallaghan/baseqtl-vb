library("argparse")
library("ggplot2")
library("ggdist")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
library("geomtextpath")
library("dplyr")
library("ggrepel")
library("ROCR")

theme_set(theme_bw())

parser <- ArgumentParser()
parser$add_argument( 
    "-m", "--model",
    default = "noGT",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-2,
    type = "double"
)

maxRhat <- 1.1 ## from baseqtl-paper repo
minEff <- 500 ## from stan docs (-ish)

args <- parser$parse_args()
tol <- args[["tolerance"]]
model <- args[["model"]]

fpath <- sprintf("fig_%1.0e", tol)
sapply(file.path(fpath, model, c("diag", "time", "estimates", "roc")), function(p) {
    dir.create(p, showWarnings = FALSE, recursive = TRUE)
})
# fpath <- "fig"

# "optimizing",
methods <- c("vb", "sampling")
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

cmdf <- mdf
cmdf <- cmdf %>% mutate(discrepancy = mean.hmc - mean.vb)
cmdf <- cmdf[cmdf$discrepancy < 10, ]

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
    # "niter",
    "n_eff.hmc", "time.hmc", "n_ase",
    "se_mean.hmc", "sd.hmc",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)
cmdf$converged <- factor(cmdf$converged)
for (type in c("discrepancy")) {
# for (type in c("discrepancy", "disc_sc_sdh", "disc_sc_sdh", "disc_sc_meanh", "disc_sc_meanv")) {
    for (x in diag_vars) {    
        geom <- if (is.numeric(cmdf[[x]])) {
            geom_point(size = 0.8, alpha = 0.6) 
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

# cmdf <- cmdf[cmdf$n_eff.hmc > 500, ]
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
    aes(khat, discrepancy) +
    geom_point(size = 0.7, alpha = 0.7) +
    geom_vline(xintercept = 0.7, linetype = "dashed") +
    labs(x = expression(hat(K)), y = "Discrepancy (ADVI - HMC)")
ggsave(
    sprintf("%s/%s/diag/disc_khat_vb.png", fpath, model),
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
    aes(time.hmc, discrepancy) +
    geom_point() +
    labs(x = "Time taken for HMC (s)", y = "Discrepancy (ADVI - HMC)")
ggsave(
    sprintf("%s/%s/diag/time_vb.png", fpath, model),
    width = 7, height = 7
)

cmdf %>%
    group_by(gene) %>%
    summarise(time = sum(time.hmc), nsnps = n()) %>%
    ggplot() +
    aes(nsnps, time) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) +
    scale_x_log10() +
    scale_y_log10() -> g
ggsave(
    sprintf("%s/%s/time/nsnp.png", fpath, model),
    width = 5, height = 5
)


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
mdfs <- cmdf[abs(cmdf$mean.vb) < 2, ]
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
stop()

## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t, ".hmc")
    levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
    sens_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            sens_vec(
                truth = factor(cmdf[[bc]], levels = c(TRUE, FALSE)),
                estimate = factor(cmdf[[column]], levels = c(TRUE, FALSE))
            )
        }
    )
    spec_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            spec_vec(
                truth = factor(cmdf[[bc]], levels = c(TRUE, FALSE)),
                estimate = factor(cmdf[[column]], levels = c(TRUE, FALSE))
            )
        }
    )
    time_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            sum(cmdf$time.hmc[cmdf[[column]]]) + sum(cmdf$time.vb)
        }
    )
    g <- ggplot() +
        aes(sens_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(cmdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        # geom_text_repel(
        #     aes(label = levs), vjust = 1, hjust = 1,
        #     segment.color = "grey70"
        # ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb))
        # +
        # ggtitle("Sensitivity vs total time for ADVI")
    ggsave(
        sprintf("%s/%s/time/time_vs_sens_vb_%s.png", fpath, model, t),
        width = 5, height = 5
    )

    g <- ggplot() +
        aes(spec_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(cmdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        # geom_text_repel(
        #     aes(label = levs), vjust = 1, hjust = 1,
        #     segment.color = "grey70"
        # ) +
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb))
        # +
        # ggtitle("Specificity vs total time for ADVI")
    ggsave(
        sprintf("%s/%s/time/time_vs_spec_vb_%s.png", fpath, model, t),
        width = 5, height = 5
    )
    # if (t == 99) stop()
    # min(time_vb[sens_vb == 1] / sum(cmdf$time.hmc))

    pred_vb <- prediction(
        predictions = cmdf[["prob"]],
        labels = factor(cmdf[[bc]], levels = c(TRUE, FALSE))
    )

    perf_auroc_vb <- performance(pred_vb, "auc")@y.values[[1]]
    perf_aupr_vb <- performance(pred_vb, "aucpr")@y.values[[1]]
    perf_pr_vb <- performance(pred_vb, "prec", "rec")
    perf_roc_vb <- performance(pred_vb, "tpr", "fpr")

    g <- ggplot() +
        geom_path(
            aes(
                perf_pr_vb@x.values[[1]], perf_pr_vb@y.values[[1]]
            )
        ) +
        lims(x = 0:1, y = 0:1) +
        labs(x = "Precision", y = "Recall") +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom") +
        ggtitle(sprintf("ADVI; AUPRC: %0.3f", perf_aupr_vb))
    ggsave(
        sprintf("%s/%s/roc/pr_vb_%s.png", fpath, model, t),
        width = 5, height = 5
    )

    g <- ggplot() +
        geom_path(
            aes(
                perf_roc_vb@x.values[[1]], perf_roc_vb@y.values[[1]]
            )
        ) +
        lims(x = 0:1, y = 0:1) +
        labs(x = "Precision", y = "Recall") +
        scale_colour_brewer(palette = "Set2", name = NULL) +
        theme(legend.position = "bottom") +
        ggtitle(sprintf("ADVI; AUROC: %0.3f", perf_auroc_vb))
    ggsave(
        sprintf("%s/%s/roc/roc_vb_%s.png", fpath, model, t),
        width = 5, height = 5
    )


    # g <- ggplot() +
    #     aes(1 - spec_vb, sens_vb) +
    #     geom_line() +
    #     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    #     labs(x = "1 - specificity", y = "Sensitivity") +
    #     ggtitle("ADVI") +
    #     lims(x = 0:1, y = 0:1)
    # ggsave(
    #     sprintf("%s/%s/roc/roc_vb_%s.png", fpath, model, t),
    #     width = 5, height = 5
    # )
}


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






# gdf <- group_by(cmdf, gene)

# ## by snp
# for (t in c(99, 95)) {
#     bc <- paste0("null.", t, ".hmc")
#     levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
#     sens_vb <- sapply(
#         levs,
#         function(x) {
#             column <- sprintf("null.%s.vb", x)
#             gdf %>%
#                 summarise(t = any(.data[[bc]]), e = any(.data[[column]])) %>%
#                 sens(
#                     ## TRUE means that the HPD interval is one side of zero (sig)
#                     ## FALSE means it overlaps zero (null)
#                     truth = factor(t, levels = c(TRUE, FALSE)),
#                     estimate = factor(e, levels = c(TRUE, FALSE))
#                 ) %>%
#                 pull(.estimate)
#         }
#     )
#     spec_vb <- sapply(
#         levs,
#         function(x) {
#             column <- sprintf("null.%s.vb", x)
#             gdf %>%
#                 summarise(t = any(.data[[bc]]), e = any(.data[[column]])) %>%
#                 spec(
#                     ## TRUE means that the HPD interval is one side of zero (sig)
#                     ## FALSE means it overlaps zero (null)
#                     truth = factor(t, levels = c(TRUE, FALSE)),
#                     estimate = factor(e, levels = c(TRUE, FALSE))
#                 ) %>%
#                 pull(.estimate)
#         }
#     )
#     time_vb <- sapply(
#         levs,
#         function(x) {
#             column <- sprintf("null.%s.vb", x)
#             gd <- gdf %>%
#                 summarise(e = any(.data[[column]]))
#             sum(cmdf$time.hmc[cmdf$gene %in% gd$gene[gd$e]]) + sum(cmdf$time.vb)
#         }
#     )
#     g <- ggplot() +
#         aes(sens_vb, time_vb) +
#         geom_line() +
#         geom_texthline(
#             yintercept = sum(cmdf$time.hmc),
#             label = "Total time without screening", 
#             vjust = -0.2,
#             linetype = "dashed"
#         ) +
#         geom_text_repel(
#             aes(label = levs), vjust = 1, hjust = 1,
#             segment.color = "grey70"
#         ) +
#         labs(x = "Sensitivity", y = "Total time (s)") +
#         ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
#         ggtitle("Sensitivity vs total time for ADVI")
#     ggsave(
#         sprintf("%s/%s/time/time_vs_sens_gene_vb_%s.png", fpath, model, t),
#         width = 5, height = 5
#     )

#     g <- ggplot() +
#         aes(spec_vb, time_vb) +
#         geom_line() +
#         geom_texthline(
#             yintercept = sum(cmdf$time.hmc),
#             label = "Total time without screening", 
#             vjust = -0.2,
#             linetype = "dashed"
#         ) +
#         geom_text_repel(
#             aes(label = levs), vjust = 1, hjust = 1,
#             segment.color = "grey70"
#         ) +
#         labs(x = "Specificity", y = "Total time (s)") +
#         ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
#         ggtitle("Specificity vs total time for ADVI")
#     ggsave(
#         sprintf("%s/%s/time/time_vs_spec_gene_vb_%s.png", fpath, model, t),
#         width = 5, height = 5
#     )

#     g <- ggplot() +
#         aes(1 - spec_vb, sens_vb) +
#         geom_line() +
#         geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#         labs(x = "1 - specificity", y = "Sensitivity") +
#         ggtitle("ADVI") +
#         lims(x = 0:1, y = 0:1)
#     ggsave(
#         sprintf("%s/%s/roc/roc_vb_gene_%s.png", fpath, model, t),
#         width = 5, height = 5
#     )
# }
