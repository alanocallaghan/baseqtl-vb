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
    "-t", "--tolerance",
    default = 1e-3,
    type = "double"
)

max_rhat <- 1.1 ## from baseqtl-paper repo
min_eff <- 500 ## from stan docs (-ish)

args <- parser$parse_args()
tol <- args[["tolerance"]]

fpath <- sprintf("fig_%1.0e", tol)
source("src/functions.R")


methods <- c("vb", "sampling")
models <- c("GT", "noGT")

dfs_both <- lapply(models,
    function(model) {
        cat(model, "\n")
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
        df_vb_hmc <- merge(
            dfs[["vb"]], dfs[["sampling"]],
            by = by, suffix = c(".vb", ".hmc"))
        print(dim(df_vb_hmc))
        df_vb_hmc$model <- model
        df_vb_hmc
    }
)
common_cols <- do.call(intersect, lapply(dfs_both, colnames))
dfs_both_sub <- lapply(dfs_both, function(x) x[, common_cols])
df_vb_hmc <- do.call(rbind, dfs_both_sub)


df_vb_hmc <- df_vb_hmc %>% mutate(
    discrepancy = mean.hmc - mean.vb,
    relative_discrepancy = discrepancy / ((mean.hmc + mean.vb) / 2),
    hpd.width.99.vb = abs(`0.5%.vb` - `99.5%.vb`),
    hpd.width.99.hmc = abs(`0.5%.hmc` - `99.5%.hmc`),
    discrepancy_99hpdi_width = hpd.width.99.hmc - hpd.width.99.vb,
    relative_discrepancy_99hpdi_width = discrepancy_99hpdi_width /
        ((hpd.width.99.hmc + hpd.width.99.vb) / 2)
)
# df_vb_hmc <- df_vb_hmc[abs(df_vb_hmc$discrepancy) < 10, ]
eff_filter <- df_vb_hmc$n_eff.hmc > min_eff
rhat_filter <- df_vb_hmc$Rhat < max_rhat
df_vb_hmc <- df_vb_hmc[eff_filter & rhat_filter, ]
## from PSIS paper, arxiv 1507.02646
## khat filter of 0.7

df_vb_hmc <- df_vb_hmc[df_vb_hmc$seed.vb == "7", ]

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


mdf_nullord <- mdf_filtered_outliers[order(mdf_filtered_outliers$nullstr99), ]
gp99_gt <- ggplot(mdf_nullord[mdf_nullord$model == "GT", ]) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(shape = 16, size = 0.75, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "HMC estimate", y = "ADVI estimate") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = "bottom")
gp99_nogt <- ggplot(mdf_nullord[mdf_nullord$model == "noGT", ]) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(shape = 16, size = 0.75, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "HMC estimate", y = "ADVI estimate") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = "bottom")



g <- plot_with_legend_below(
    gp99_gt, gp99_nogt
    labels = "AUTO"
)

# g <- cowplot::plot_grid(gp99_gt, gp99_nogt, labels="AUTO")
ggsave(g,
    file = sprintf("%s/point-estimates-99-GT-noGT.pdf", fpath),
    width = 5.5, height = 3.5
)
