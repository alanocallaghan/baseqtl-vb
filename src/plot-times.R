library("ggplot2")
library("ggdist")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
library("geomtextpath")
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

# fpath <- sprintf("fig_%1.0e", tol)
fpath <- "fig"

# "optimizing",
methods <- c("vb", "sampling")
dfs <- lapply(methods,
    function(method) {
        cat(method, "\n")
        combfile <- sprintf("rds/%s/%s_combined.rds", model, method)
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
# mdf <- merge(dfs[[1]], dfs[[2]], by = by, suffix = c(".vb", ".map"))
# mdf <- merge(mdf, dfs[[3]], by = by, suffix = c("", ".hmc"))
mdf <- merge(dfs[["vb"]], dfs[["sampling"]], by = by, suffix = c(".vb", ".hmc"))

mdf <- mdf %>%
    mutate(disc_vb = mean.vb - mean.hmc
    # , disc_map = mean.map - mean
)



# mtol <- sprintf("vb_%1.0e", tol)
mtol <- sprintf("vb")

cmdf <- mdf
cmdf <- cmdf %>% mutate(discrepancy = mean.hmc - mean.vb)

x <- cmdf %>%
    arrange(-abs(disc_vb)) %>%
    top_n(50, disc_vb) %>%
    select(vb = mean.vb, hmc = mean.hmc, discrepancy = disc_vb, snp, gene)

saveRDS(x,
    sprintf("rds/%s_discrepancies_vb_%1.0e.rds", model, tol)
)

mname <- "ADVI"
method <- "vb"
r <- range(c(cmdf$mean.hmc, cmdf$mean.vb))
## not se mean but other hmc se
diag_vars <- c(
    "gene",
    "Rhat",
    "khat",
    "n_eff.hmc", "time.hmc", "n_ase",
    "se_mean.hmc", "sd.hmc",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)
for (x in diag_vars) {    
    g <- ggplot(cmdf) +
        aes_string(x, "abs(discrepancy)") +
        geom_point(size = 0.8) +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")
        ) +
        labs(x = x, y = "Discrepancy")

    ggsave(
        sprintf("%s/%s/diag/disc_%s_%s.png", fpath, model, x, method),
        width = 7, height = 7
    )
}

g <- ggplot(cmdf) +
    aes(time.hmc, disc_vb) +
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
    geom_smooth(method = "lm", formula = y ~ x) -> g
ggsave(
    sprintf("%s/%s/time/nsnp.png", fpath, model),
    width = 7, height = 7
)


cmdf <- cmdf[cmdf$n_eff.hmc > minEff & cmdf$Rhat < maxRhat, ]
## from PSIS paper, arxiv 1507.02646
# cmdf <- cmdf[cmdf$khat < 0.7, ]

gdf <- group_by(cmdf, gene)

lab_str <- "HMC: %s\nVB: %s\n"
cmdf$nullstr95 <- sprintf(lab_str,
    ifelse(cmdf$null.95.hmc, "no", "yes"),
    ifelse(cmdf$null.95.vb, "no", "yes")
)
cmdf$nullstr99 <- sprintf(lab_str,
    ifelse(cmdf$null.99.hmc, "no", "yes"),
    ifelse(cmdf$null.99.vb, "no", "yes")
)

null_levs <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "no",  "yes", "yes")
)
null_ord <- null_levs[c(1, 3, 2, 4)]
scale <- scale_colour_brewer(palette = "Paired", limits = null_levs, name = NULL)
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
    labs(x = "MCMC estimate", y = "VB estimate") +
    scale
ggsave(
    sprintf("%s/%s/estimates/point-estimates-95.png", fpath, model),
    width = 7, height = 7
)

g <- ggplot(mdfs[order(mdfs$nullstr99), ]) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    lims(x = lim, y = lim) +
    labs(x = "MCMC estimate", y = "VB estimate") +
    scale
ggsave(
    sprintf("%s/%s/estimates/point-estimates-99.png", fpath, model),
    width = 7, height = 7
)


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
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
        ggtitle("Sensitivity vs total time for ADVI")
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
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
        ggtitle("Specificity vs total time for ADVI")
    ggsave(
        sprintf("%s/%s/time/time_vs_spec_vb_%s.png", fpath, model, t),
        width = 5, height = 5
    )

    g <- ggplot() +
        aes(1 - spec_vb, sens_vb) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("VB") +
        lims(x = 0:1, y = 0:1)
    ggsave(
        sprintf("%s/%s/roc/roc_vb_%s.png", fpath, model, t),
        width = 5, height = 5
    )
}



gdf <- group_by(cmdf, gene)

## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t, ".hmc")
    levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
    sens_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            gdf %>%
                summarise(t = any(.data[[bc]]), e = any(.data[[column]])) %>%
                sens(
                    ## TRUE means that the HPD interval is one side of zero (sig)
                    ## FALSE means it overlaps zero (null)
                    truth = factor(t, levels = c(TRUE, FALSE)),
                    estimate = factor(e, levels = c(TRUE, FALSE))
                ) %>%
                pull(.estimate)
        }
    )
    spec_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            gdf %>%
                summarise(t = any(.data[[bc]]), e = any(.data[[column]])) %>%
                spec(
                    ## TRUE means that the HPD interval is one side of zero (sig)
                    ## FALSE means it overlaps zero (null)
                    truth = factor(t, levels = c(TRUE, FALSE)),
                    estimate = factor(e, levels = c(TRUE, FALSE))
                ) %>%
                pull(.estimate)
        }
    )
    time_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            gd <- gdf %>%
                summarise(e = any(.data[[column]]))
            sum(cmdf$time.hmc[cmdf$gene %in% gd$gene[gd$e]]) + sum(cmdf$time.vb)
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
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
        ggtitle("Sensitivity vs total time for ADVI")
    ggsave(
        sprintf("%s/%s/time/time_vs_sens_gene_vb_%s.png", fpath, model, t),
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
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
        ggtitle("Specificity vs total time for ADVI")
    ggsave(
        sprintf("%s/%s/time/time_vs_spec_gene_vb_%s.png", fpath, model, t),
        width = 5, height = 5
    )

    g <- ggplot() +
        aes(1 - spec_vb, sens_vb) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("VB") +
        lims(x = 0:1, y = 0:1)
    ggsave(
        sprintf("%s/%s/roc/roc_vb_gene_%s.png", fpath, model, t),
        width = 5, height = 5
    )
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
