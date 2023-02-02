library("argparse")
library("ggplot2")
library("cowplot")
library("ggdist")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
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

max_rhat <- 1.1 ## from baseqtl-paper repo
min_eff <- 500 ## from stan docs (-ish)

args <- parser$parse_args()
tol <- args[["tolerance"]]
model <- args[["model"]]

fpath <- sprintf("fig_%1.0e", tol)
source("src/functions.R")
mkfigdir(fpath, model)


if (model == "GT") {
    methods <- c("vb", "sampling")
} else {
    methods <- c("vb", "sampling")
}
dfs <- lapply(
    methods,
    function(method) {
        cat(method, "\n")
        mtol <- mtol(method, tol)
        combfile <- sprintf("rds/%s/%s_combined.rds", model, mtol)
        if (file.exists(combfile)) {
            return(readRDS(combfile))
        } else {
            stop("File not found:", combfile)
        }
    }
)

dfs[] <- lapply(dfs, add_nulls)


names(dfs) <- methods
by <- if (model == "GT") {
    c(
        "test", "gene", "snp", "n_tot", "n_ase", "p_het", "n_tot",
        "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
    )
} else {
    c(
        "test", "gene", "snp", "condition",
        "n_tot", "n_ase", "mean_ase", "sd_ase", "mean_count", "sd_count",
        "n_wt", "n_het", "n_hom", "p_het", "n_tot"
    )
}
dfs[["vb"]] <- dfs[["vb"]][dfs[["vb"]]$seed == "7", ]
df_vb_hmc <- merge(
    dfs[["vb"]], dfs[["sampling"]],
    by = by, suffix = c(".vb", ".hmc"))
print(dim(df_vb_hmc))



df_vb_hmc <- df_vb_hmc %>% mutate(
    discrepancy = mean.hmc - mean.vb,
    relative_discrepancy = discrepancy / ((mean.hmc + mean.vb) / 2),
    hpd.width.99.vb = abs(`0.5%.vb` - `99.5%.vb`),
    hpd.width.99.hmc = abs(`0.5%.hmc` - `99.5%.hmc`),
    discrepancy_99hpdi_width = hpd.width.99.hmc - hpd.width.99.vb,
    relative_discrepancy_99hpdi_width = discrepancy_99hpdi_width /
        ((hpd.width.99.hmc + hpd.width.99.vb) / 2)
)
df_vb_hmc <- df_vb_hmc[abs(df_vb_hmc$discrepancy) < 10, ]
eff_filter <- df_vb_hmc$n_eff.hmc > min_eff
rhat_filter <- df_vb_hmc$Rhat < max_rhat
df_vb_hmc <- df_vb_hmc[eff_filter & rhat_filter, ]
## from PSIS paper, arxiv 1507.02646
## khat filter of 0.7


if (model == "GT") {
    disc_df <- df_vb_hmc %>%
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
} else {
    disc_df <- df_vb_hmc %>%
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
            condition = condition,
            snp = snp,
            gene = gene
        )
}

saveRDS(
    disc_df,
    sprintf("rds/%s_discrepancies_vb_%1.0e.rds", model, tol)
)


lab_str <- "HMC: %s\nADVI: %s\n"
df_vb_hmc$nullstr99 <- sprintf(
    lab_str,
    ifelse(df_vb_hmc$null.99.hmc, "yes", "no"),
    ifelse(df_vb_hmc$null.99.vb, "yes", "no")
)

null_levs <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "no", "yes", "yes")
)
null_ord <- null_levs[c(1, 3, 2, 4)]


scale <- scale_colour_manual(
    values = setNames(
        c("#fb9a99", "#e31a1c", "#a6cee3", "#1f78b4"),
        null_levs
    ),
    drop = TRUE,
    name = "Significant"
)
df_vb_hmc$nullstr99 <- factor(df_vb_hmc$nullstr99, levels = null_ord)


mdf_filtered_outliers <- df_vb_hmc[abs(df_vb_hmc$mean.vb) < 2, ]
lim <- range(c(mdf_filtered_outliers$mean.hmc, mdf_filtered_outliers$mean.vb))
modname <- if (model == "GT") {
    "with known genotypes"
} else {
    "with unknown genotypes"
}


mdf_nullord <- mdf_filtered_outliers[order(mdf_filtered_outliers$nullstr99), ]
gp99 <- ggplot(mdf_nullord) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(shape = 16, size = 0.75, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "HMC estimate", y = "ADVI estimate") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = "bottom")
# ggsave(
#     sprintf("%s/%s/estimates/point-estimates-99.pdf", fpath, model),
#     width = 5, height = 5
# )
make_crosstab(
    x = ifelse(mdf_filtered_outliers$null.99.hmc, "Significant", "Null"),
    y = ifelse(mdf_filtered_outliers$null.99.vb, "Significant", "Null"),
    xn = "HMC",
    yn = "ADVI",
    caption = paste(
        "Confusion matrix of HMC and ADVI significance calls at 99\\%",
        "threshold for BaseQTL", modname
    ),
    label = sprintf("tab:%s-xtab-99", model),
    file = sprintf("table/%s-xtab-99.tex", model)
)
make_crosstab(
    x = ifelse(mdf_filtered_outliers$null.99.hmc, "Significant", "Null"),
    y = ifelse(mdf_filtered_outliers$null.99.vb, "Significant", "Null"),
    xn = "HMC",
    yn = "ADVI",
    prop = TRUE,
    caption = paste(
        "Frequencies of HMC and ADVI significance calls at 99\\% threshold",
        "for BaseQTL", modname
    ),
    label = sprintf("tab:%s-xtab-prop-99", model),
    file = sprintf("table/%s-xtab-prop-99.tex", model)
)


ind_discord <- !(mdf_filtered_outliers$null.99.hmc == mdf_filtered_outliers$null.99.vb)
mdf_discordant <- mdf_filtered_outliers[ind_discord, ]


gpd99 <- ggplot(mdf_discordant) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(shape = 16, size = 0.75, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "HMC estimate", y = "ADVI estimate") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = "bottom")
ggsave(
    sprintf("%s/%s/estimates/point-estimates-diff-99.pdf", fpath, model),
    width = 5, height = 5
)


g <- plot_with_legend_below(
    gp99 + annotate(
        geom = "text",
        x = -1.5, y = 1.5,
        hjust = -0.2, vjust = -0.2,
        fontface = "bold",
        label = "A"
    ),
    gpd99 + annotate(
        geom = "text",
        x = -1.5, y = 1.5,
        hjust = -0.2, vjust = -0.2,
        fontface = "bold",
        label = "B"
    ) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        ),
    labels = NULL,
    common_x = TRUE,
    xlab = "HMC estimate"
)

ggsave(
    sprintf("%s/%s/estimates/point-estimates-both-99.pdf", fpath, model),
    width = 5.5, height = 3.5
)


level <- 0.99
## by snp
bc <- paste0("null.99.hmc")

peps <- seq(1, 0, by = -0.01)
peps <- peps[peps != 1]
sens_vb <- sapply(
    peps,
    function(x) {
        sens_vec(
            truth = factor(df_vb_hmc[[bc]], levels = c(TRUE, FALSE)),
            estimate = factor(df_vb_hmc$PEP.vb > x, levels = c(TRUE, FALSE))
        )
    }
)

spec_vb <- sapply(
    peps,
    function(x) {
        spec_vec(
            truth = factor(df_vb_hmc[[bc]], levels = c(TRUE, FALSE)),
            estimate = factor(df_vb_hmc$PEP.vb > x, levels = c(TRUE, FALSE))
        )
    }
)
time_vb <- sapply(
    peps,
    function(x) {
        sum(df_vb_hmc$time.hmc[df_vb_hmc$PEP.vb > x]) + sum(df_vb_hmc$time.vb)
    }
)
sens_spec_time_df <- data.frame(
    pep = peps,
    sens = sens_vb,
    spec = spec_vb,
    time = time_vb
)

total_time_without_screening <- sum(df_vb_hmc$time.hmc)
ind_sens_1 <- which(sens_spec_time_df$sens == 1)[[1]]
time_saved <- (total_time_without_screening - sens_spec_time_df$time[ind_sens_1]) / total_time_without_screening
prob_used <- sens_spec_time_df$pep[ind_sens_1]
print(paste("Time reduction:", time_saved))
print(paste("PEP used:", prob_used))


gtime <- ggplot(sens_spec_time_df) +
    aes(sens, time / 3600) +
    geom_line(na.rm = TRUE) +
    geom_hline(
        yintercept = total_time_without_screening / 3600,
        linetype = "dashed"
    ) +
    annotate(
        geom = "text",
        x = 0.96,
        y = total_time_without_screening / 3600,
        label = "Total time\nwithout screening",
        vjust = -0.2,
    ) +
    ylim(0, max(total_time_without_screening, sens_spec_time_df$time_vb) / 3600) +
    labs(x = "Sensitivity", y = "Total time (hr)")
# ggsave(
#     sprintf("%s/%s/time/time_vs_sens_vb_%s.pdf", fpath, model, level),
#     width = 5, height = 5
# )

prob <- df_vb_hmc[["PEP.vb"]] %||% df_vb_hmc[["PEP"]]
pred_vb <- prediction(
    predictions = prob,
    labels = factor(df_vb_hmc[[bc]], levels = c(TRUE, FALSE))
)


perf_auroc_vb <- performance(pred_vb, "auc")@y.values[[1]]
perf_roc_vb <- performance(pred_vb, "tpr", "fpr")

g <- ggplot() +
    geom_path(
        aes(
            perf_roc_vb@x.values[[1]], perf_roc_vb@y.values[[1]]
        )
    ) +
    lims(x = 0:1, y = 0:1) +
    labs(x = "True positive rate", y = "False positive rate") +
    scale_colour_brewer(palette = "Set2", name = NULL) +
    theme(legend.position = "bottom") +
    ggtitle(sprintf("ADVI; AUROC: %0.3f", perf_auroc_vb))
# ggsave(
#     sprintf("%s/%s/roc/roc_vb_%s.pdf", fpath, model, level),
#     width = 5, height = 5
# )
groc <- g + ggtitle(NULL)
# ggsave(
#     sprintf("%s/%s/roc/roc_vb_%s.pdf", fpath, model, level),
#     width = 3.5, height = 4
# )

p <- plot_grid(groc, gtime, labels = "AUTO")
ggsave(
    sprintf("%s/%s/time/time_roc_%s.pdf", fpath, model, level),
    width = 5.5, height = 3.5
)
print("Success!")










################################################################################
## Unused graphs
################################################################################


# ## time comparison
# g <- ggplot(df_vb_hmc) +
#     aes(time.hmc, time.vb) +
#     geom_pointdensity(shape = 16, size = 0.8) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     scale_x_log10() +
#     scale_y_log10() +
#     scale_colour_viridis() +
#     theme(legend.position = "none") +
#     labs(x = "Time (s) for HMC", y = "Time (s) for ADVI")
# ggsave(
#     sprintf("%s/%s/time/time_hmc_vs_vb.pdf", fpath, model),
#     width = 5, height = 5
# )

# cn <- lapply(dfs, colnames)
# cn <- Reduce(intersect, cn)
# dfs_int <- lapply(dfs, function(x) x[, cn])
# dfs_int[] <- lapply(
#     names(dfs_int),
#     function(n) {
#         dfs_int[[n]]$method <- n
#         dfs_int[[n]]
#     }
# )
# df_int <- do.call(rbind, dfs_int)


# g <- ggplot(df_int) +
#     aes(x = time, colour = method) +
#     geom_density() +
#     geom_vline(
#         aes(
#             xintercept = mean(time[method == "vb"]),
#             colour = "vb"
#         ),
#         linetype = "dashed"
#     ) +
#     geom_vline(
#         aes(
#             xintercept = mean(time[method == "sampling"]),
#             colour = "sampling"
#         ),
#         linetype = "dashed"
#     ) +
#     scale_x_log10(name = "Time (s)") +
#     scale_colour_brewer(palette = "Set1", name = "Method") +
#     ylab("Density") +
#     theme(legend.position = "bottom")
# ggsave(
#     sprintf("%s/%s/time/time_comparison.pdf", fpath, model),
#     width = 5, height = 5
# )


# ## comparing absolute and relative differences
# mdf_filtered_outliers <- df_vb_hmc[abs(df_vb_hmc$mean.vb) < 2, ]
# lim <- range(c(mdf_filtered_outliers$mean.hmc, mdf_filtered_outliers$mean.vb))



# ## point estimates with HPD shown
# mdf_discordant <- mdf_discordant %>% mutate(
#     hpd_width_hmc = `99.5%.hmc` - `0.5%.hmc`,
#     hpd_width_vb = `99.5%.vb` - `0.5%.vb`,
#     hpd_width_ratio = hpd_width_hmc / hpd_width_vb
# )

# limv <- range(c(
#     mdf_discordant[["0.5%.hmc"]], mdf_discordant[["99.5%.hmc"]],
#     mdf_discordant[["0.5%.vb"]], mdf_discordant[["99.5%.vb"]]
# ))
# g <- ggplot(mdf_discordant) +
#     geom_pointrange(
#         aes(
#             x = mean.hmc,
#             xmin = `0.5%.hmc`,
#             xmax = `99.5%.hmc`,
#             colour = nullstr99,
#             y = mean.vb
#         ),
#         alpha = 0.2,
#         fatten = 1,
#         pch = 16
#     ) +
#     geom_pointrange(
#         aes(
#             x = mean.hmc,
#             y = mean.vb,
#             ymin = `0.5%.vb`,
#             colour = nullstr99,
#             ymax = `99.5%.vb`
#         ),
#         alpha = 0.2,
#         fatten = 1,
#         pch = 16
#     ) +
#     lims(x = limv, y = limv) +
#     labs(x = "HMC estimate", y = "ADVI estimate") +
#     theme(legend.position = "below") +
#     # guides(colour = guide_legend(override.aes = list(size = 2))) +
#     scale
# ggsave(
#     sprintf("%s/%s/estimates/point-estimates-hpd-99.pdf", fpath, model),
#     width = 4, height = 4
# )


# seqr <- seq(lim[[1]], lim[[2]], length.out = 200)
# grid <- expand.grid(x = seqr, y = seqr)
# grid$absdiff <- sapply(seq_len(nrow(grid)), function(i) {
#     abs(grid[i, 1] - grid[i, 2])
# })
# grid$reldiff1 <- sapply(seq_len(nrow(grid)), function(i) {
#     abs((grid[i, 1] - grid[i, 2]) / ((grid[i, 1] + grid[i, 2]) / 2))
# })
# grid$reldiff2 <- sapply(seq_len(nrow(grid)), function(i) {
#     abs((grid[i, 1] - grid[i, 2]) / grid[i, 1])
# })
# ga <- ggplot(grid) +
#     aes(x = x, y = y, fill = absdiff) +
#     geom_tile() +
#     scale_fill_viridis(name = "Absolute difference") +
#     theme(legend.position = "bottom")
# gr1 <- ggplot(grid) +
#     aes(x = x, y = y, fill = reldiff1) +
#     geom_tile() +
#     scale_fill_viridis(name = "Relative (to mean) absolute difference") +
#     theme(legend.position = "bottom")
# gr2 <- ggplot(grid) +
#     aes(x = x, y = y, fill = reldiff2) +
#     geom_tile() +
#     scale_fill_viridis(name = "Relative (to x) absolute difference") +
#     theme(legend.position = "bottom")

# pp <- cowplot::plot_grid(ga, gr1, gr2, nrow = 1, labels = "AUTO")
# ggsave(
#     sprintf("%s/%s/estimates/abs-rel-diff-comp.pdf", fpath, model),
#     width = 18, height = 6
# )


# dord <- order(abs(mdf_filtered_outliers$discrepancy))
# mdf_dord <- mdf_filtered_outliers[dord, ]
# pd1 <- ggplot(mdf_dord) +
#     aes(mean.hmc, mean.vb, colour = abs(discrepancy)) +
#     geom_abline(
#         slope = 1, intercept = 0, linetype = "dashed", colour = "grey80"
#     ) +
#     geom_point(shape = 16, size = 0.5, alpha = 0.7) +
#     lims(x = lim, y = lim) +
#     labs(x = "HMC estimate", y = "ADVI estimate") +
#     scale_colour_viridis(name = "Absolute difference") +
#     guides(colour = guide_legend(override.aes = list(size = 2))) +
#     theme(legend.position = "bottom")

# rdord <- order(abs(mdf_filtered_outliers$relative_discrepancy))
# mdf_rdord <- mdf_filtered_outliers[rdord, ]
# pd2 <- ggplot(mdf_rdord) +
#     aes(mean.hmc, mean.vb, colour = abs(relative_discrepancy)) +
#     geom_abline(
#         slope = 1, intercept = 0, linetype = "dashed", colour = "grey80"
#     ) +
#     geom_point(shape = 16, size = 0.75, alpha = 0.7) +
#     lims(x = lim, y = lim) +
#     labs(x = "HMC estimate", y = "ADVI estimate") +
#     scale_colour_viridis(name = "Relative absolute difference", trans = "log10") +
#     guides(colour = guide_legend(override.aes = list(size = 2))) +
#     theme(legend.position = "bottom")
# pp <- cowplot::plot_grid(pd1, pd2, labels = "AUTO")
# ggsave(
#     sprintf("%s/%s/estimates/point-estimates-rel-abs-disc.pdf", fpath, model),
#     width = 12, height = 6
# )

# g <- ggplot(sens_spec_time_df) +
#     aes(spec, time) +
#     geom_line() +
#     geom_hline(
#         yintercept = total_time_without_screening,
#         linetype = "dashed"
#     ) +
#     annotate(
#         geom = "text",
#         x = median(sens_spec_time_df$spec_vb),
#         y = total_time_without_screening / 3600,
#         label = "Total time\nwithout screening",
#         vjust = -0.2,
#     ) +
#     labs(x = "Specificity", y = "Total time (hr)") +
#     ylim(0, max(total_time_without_screening, time_vb))
# ggsave(
#     sprintf("%s/%s/time/time_vs_spec_vb_%s.pdf", fpath, model, level),
#     width = 5, height = 5
# )


# perf_aupr_vb <- performance(pred_vb, "aucpr")@y.values[[1]]
# perf_pr_vb <- performance(pred_vb, "prec", "rec")


# g <- ggplot() +
#     geom_path(
#         aes(
#             perf_pr_vb@x.values[[1]], perf_pr_vb@y.values[[1]]
#         )
#     ) +
#     lims(x = 0:1, y = 0:1) +
#     labs(x = "Precision", y = "Recall") +
#     scale_colour_brewer(palette = "Set2", name = NULL) +
#     theme(legend.position = "bottom") +
#     ggtitle(sprintf("ADVI; AUPRC: %0.3f", perf_aupr_vb))
# ggsave(
#     sprintf("%s/%s/roc/pr_vb_%s.pdf", fpath, model, level),
#     width = 4, height = 4
# )






################################################################################
#
# Pathfinder figs
#
################################################################################

# stop()
# if (model == "GT") {
#     df_pf_hmc <- merge(dfs[["pathfinder"]], dfs[["sampling"]], by = by, suffix = c(".vb", ".hmc"))
#     dfs[["pathfinder_parallel"]] <- dfs[["pathfinder_parallel"]][, !duplicated(colnames(dfs[["pathfinder_parallel"]]))]
#     df_pfp_hmc <- merge(dfs[["pathfinder_parallel"]], dfs[["sampling"]], by = by, suffix = c(".vb", ".hmc"))
#     lab_str <- "HMC: %s\nPathfinder: %s\n"
#     df_pf_hmc$nullstr99 <- sprintf(lab_str,
#         ifelse(df_pf_hmc$null.99.hmc, "yes", "no"),
#         ifelse(df_pf_hmc$null.99.vb, "yes", "no")
#     )
#     null_levs <- sprintf(
#         lab_str,
#         c("no", "yes", "yes", "no"),
#         c("no", "no",  "yes", "yes")
#     )
#     null_ord <- null_levs[c(1, 3, 2, 4)]

#     scale <- scale_colour_manual(
#         values = setNames(c("#fb9a99", "#e31a1c", "#a6cee3", "#1f78b4"), null_levs),
#         drop = TRUE,
#         name = "Significance"
#     )

#     df_pf_hmc$nullstr99 <- factor(df_pf_hmc$nullstr99, levels = null_ord)


#     mdf_filtered_outliers <- df_pf_hmc[abs(df_pf_hmc$mean.vb) < 2, ]
#     lim <- range(c(mdf_filtered_outliers$mean.hmc, mdf_filtered_outliers$mean.vb))

#     gp99 <- ggplot(mdf_filtered_outliers[order(mdf_filtered_outliers$nullstr99), ]) +
#         aes(mean.hmc, mean.vb, colour = nullstr99) +
#         geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#         geom_point(shape = 16, size = 0.75, alpha = 0.7) +
#         lims(x = lim, y = lim) +
#         labs(x = "HMC estimate", y = "Pathfinder estimate") +
#         scale +
#         guides(colour = guide_legend(override.aes = list(size = 2))) +
#         theme(legend.position = "bottom")
#     ggsave(
#         sprintf("%s/%s/estimates/point-estimates-99-pf.pdf", fpath, model),
#         width = 6, height = 6
#     )
#     make_crosstab(
#         x = ifelse(mdf_filtered_outliers$null.99.hmc, "Significant", "Null"),
#         y = ifelse(mdf_filtered_outliers$null.99.vb, "Significant", "Null"),
#         xn = "HMC",
#         yn = "Pathfinder",
#         caption = sprintf("Confusion matrix of HMC and Pathfinder significance calls at 99\\%% threshold for BaseQTL %s", modname),
#         label = sprintf("tab:%s-xtab-99-pf", model),
#         file = sprintf("table/%s-xtab-99-pf.tex", model)
#     )
#     make_crosstab(
#         x = ifelse(mdf_filtered_outliers$null.99.hmc, "Significant", "Null"),
#         y = ifelse(mdf_filtered_outliers$null.99.vb, "Significant", "Null"),
#         xn = "HMC",
#         yn = "ADVI",
#         prop = TRUE,
#         caption = sprintf("Frequencies of HMC and Pathfinder significance calls at 99\\%% threshold for BaseQTL %s", modname),
#         label = sprintf("tab:%s-xtab-prop-99-pf", model),
#         file = sprintf("table/%s-xtab-prop-99-pf.tex", model)
#     )





#     lab_str <- "HMC: %s\nPathfinder: %s\n"
#     df_pfp_hmc$nullstr99 <- sprintf(lab_str,
#         ifelse(df_pfp_hmc$null.99.hmc, "yes", "no"),
#         ifelse(df_pfp_hmc$null.99.vb, "yes", "no")
#     )
#     df_pfp_hmc$nullstr99 <- factor(df_pfp_hmc$nullstr99, levels = null_ord)




#     mdf_filtered_outliers <- df_pfp_hmc[abs(df_pfp_hmc$mean.vb) < 2, ]
#     lim <- range(c(mdf_filtered_outliers$mean.hmc, mdf_filtered_outliers$mean.vb))

#     gp99 <- ggplot(mdf_filtered_outliers[order(mdf_filtered_outliers$nullstr99), ]) +
#         aes(mean.hmc, mean.vb, colour = nullstr99) +
#         geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#         geom_point(shape = 16, size = 0.75, alpha = 0.7) +
#         lims(x = lim, y = lim) +
#         labs(x = "HMC estimate", y = "Pathfinder estimate") +
#         scale +
#         guides(colour = guide_legend(override.aes = list(size = 2))) +
#         theme(legend.position = "bottom")
#     ggsave(
#         sprintf("%s/%s/estimates/point-estimates-99-pfp.pdf", fpath, model),
#         width = 6, height = 6
#     )
#     make_crosstab(
#         x = ifelse(mdf_filtered_outliers$null.99.hmc, "Significant", "Null"),
#         y = ifelse(mdf_filtered_outliers$null.99.vb, "Significant", "Null"),
#         xn = "HMC",
#         yn = "Pathfinder",
#         caption = sprintf("Confusion matrix of HMC and Pathfinder significance calls at 99\\%% threshold for BaseQTL %s", modname),
#         label = sprintf("tab:%s-xtab-99-pfp", model),
#         file = sprintf("table/%s-xtab-99-pfp.tex", model)
#     )
#     make_crosstab(
#         x = ifelse(mdf_filtered_outliers$null.99.hmc, "Significant", "Null"),
#         y = ifelse(mdf_filtered_outliers$null.99.vb, "Significant", "Null"),
#         xn = "HMC",
#         yn = "ADVI",
#         prop = TRUE,
#         caption = sprintf("Frequencies of HMC and Pathfinder significance calls at 99\\%% threshold for BaseQTL %s", modname),
#         label = sprintf("tab:%s-xtab-prop-99-pfp", model),
#         file = sprintf("table/%s-xtab-prop-99-pfp.tex", model)
#     )

#     ## time comparison
#     g <- ggplot(df_pfp_hmc) +
#         aes(time.hmc, time.vb) +
#         geom_pointdensity(shape = 16, size = 0.8) +
#         geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#         scale_x_log10() +
#         scale_y_log10() +
#         scale_colour_viridis() +
#         theme(legend.position = "none") +
#         labs(x = "Time (s) for HMC", y = "Time (s) for Pathfinder")
#     ggsave(
#         sprintf("%s/%s/time/time_hmc_vs_pfp.pdf", fpath, model),
#         width = 5, height = 5
#     )
# }

