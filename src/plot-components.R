library("ggplot2")
library("ggdist")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
library("dplyr")
library("argparse")

theme_set(theme_bw())

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "noGT",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-3,
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

# "optimizing",
methods <- c("vb", "sampling")
dfs <- lapply(
    methods,
    function(method) {
        cat(method, "\n")
        mtol <- mtol(method, tol)
        combfile <- sprintf("rds/%s/components/%s_combined.rds", model, mtol)
        if (file.exists(combfile)) {
            return(readRDS(combfile))
        }
    }
)

# dfs[] <- lapply(dfs,
#     function(df) {
#         ## TRUE means that the HPD interval is one side of zero (sig)
#         ## FALSE means it overlaps zero (null)
#         df$null.95 <- sign(df$"2.5%") == sign(df$"97.5%")
#         df$null.90 <- sign(df$"5.0%") == sign(df$"95.0%")
#         df$null.80 <- sign(df$"10.0%") == sign(df$"90.0%")
#         df$null.70 <- sign(df$"15.0%") == sign(df$"85.0%")
#         df$null.60 <- sign(df$"20.0%") == sign(df$"80.0%")
#         df$null.50 <- sign(df$"25.0%") == sign(df$"75.0%")
#         df$null.40 <- sign(df$"30.0%") == sign(df$"70.0%")
#         df$null.30 <- sign(df$"35.0%") == sign(df$"65.0%")
#         df$null.20 <- sign(df$"40.0%") == sign(df$"60.0%")
#         df$null.10 <- sign(df$"45.0%") == sign(df$"55.0%")
#         df
#     }
# )

names(dfs) <- methods
by <- if (model == "GT") {
    c(
        "test", "gene", "snp", "n_tot", "n_ase", "component",
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
mdf <- mdf %>% mutate(discrepancy = mean.hmc - mean.vb)




dv <- dfs$vb[, c("mean", "component", "snp", "gene")]
dv <- dv[order(dv$component, dv$mean), ]
dv$test <- paste(dv$gene, dv$snp, sep = "_")
levs <- dv[dv$component == "both", "test"]
dv$testn <- factor(dv$test, levels = levs)

g <- ggplot(dv) +
    aes(testn, mean, colour = component) +
    geom_point(alpha = 0.5) +
    scale_colour_brewer(palette = "Set2") +
    labs(x = "Association") +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )
ggsave(
    sprintf("%s/%s/estimates/component-estimates-colour.png", fpath, model),
    width = 7, height = 7
)

ds <- dfs$sampling[, c("mean", "component", "snp", "gene")]
ds <- ds[order(ds$component, ds$mean), ]
ds$test <- paste(ds$gene, ds$snp, sep = "_")
levs <- ds[ds$component == "both", "test"]
ds$testn <- factor(ds$test, levels = levs)

g <- ggplot(ds) +
    aes(testn, mean, colour = component) +
    geom_point(alpha = 0.5) +
    scale_colour_brewer(palette = "Set2") +
    labs(x = "Association") +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )
ggsave(
    sprintf("%s/%s/estimates/component-estimates-colour.png", fpath, model),
    width = 7, height = 7
)

mdf$label <- gsub("inter", "NB", mdf$component)
mdf$label <- gsub("intra", "BB", mdf$label)

lims <- range(c(mdf$mean.hmc, mdf$mean.vb))

g <- ggplot(mdf) +
    aes(mean.hmc, mean.vb) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.5, size = 0.7) +
    facet_wrap(~label) +
    lims(x = lims, y = lims) +
    labs(x = "HMC estimate", y = "VB estimate")
ggsave(
    sprintf("%s/%s/estimates/component-estimates-facet.png", fpath, model),
    width = 10, height = 6
)

g <- ggplot(mdf) +
    aes(mean.hmc, discrepancy) +
    geom_point(alpha = 0.5, size = 0.7) +
    facet_wrap(~label) +
    labs(x = "HMC estimate", y = "Discrepancy")
ggsave(
    sprintf("%s/%s/diag/discrepancy-components.png", fpath, model),
    width = 10, height = 6
)


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

saveRDS(
    x,
    sprintf("rds/%s_component_discrepancies_vb_%1.0e.rds", model, tol)
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
            geom_point(size = 0.8, alpha = 0.7)
        } else if (is.logical(cmdf[[x]]) || is.character(cmdf[[x]]) || is.factor(cmdf[[x]])) {
            geom_boxplot(fill = "grey80")
        }
        g <- ggplot(cmdf) +
            aes_string(x, sprintf("abs(%s)", type), colour = "component") +
            geom +
            geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
            # scale_y_log10() +
            labs(x = x, y = "Discrepancy")

        ggsave(
            sprintf("%s/%s/diag/comp_%s_%s_%s.png", fpath, model, type, x, method),
            width = 7, height = 7
        )
    }
}
