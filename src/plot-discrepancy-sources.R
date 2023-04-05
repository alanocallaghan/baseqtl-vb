library("ggplot2")
library("ggdist")
library("ggrastr")
library("ggpointdensity")
library("dplyr")
library("viridis")
library("argparse")
library("latex2exp")
library("cowplot")

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "GT",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-2,
    type = "double"
)

options(mc.cores = 8)

args <- parser$parse_args()
source("src/functions.R")
theme_set(theme_bw())

tol <- args[["tolerance"]]
model <- args[["model"]]
mname <- "ADVI"
method <- "vb"

tol_str <- sprintf("%s_%1.0e", method, tol)

fpath <- sprintf("fig_%1.0e", tol)
mkfigdir(fpath, model)


infile <- sprintf("rds/%s_discrepancies_vb_%1.0e.rds", model, tol)
# if (!file.exists(infile)) {

maxRhat <- 1.1 ## from baseqtl-paper repo
minEff <- 500 ## from stan docs (-ish)


# "optimizing",
methods <- c("vb", "sampling")
dfs <- lapply(
    methods,
    function(method) {
        cat(method, "\n")
        mtol <- mtol(method, tol)
        combfile <- sprintf("rds/%s/%s_combined.rds", model, mtol)
        if (file.exists(combfile)) {
            return(readRDS(combfile))
        }
    }
)

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
mdf <- merge(dfs[["vb"]], dfs[["sampling"]], by = by, suffix = c(".vb", ".hmc"))

cmdf <- mdf
cmdf <- cmdf %>% mutate(
    discrepancy = mean.hmc - mean.vb,
    abs_discrepancy = abs(discrepancy),
    relative_discrepancy = discrepancy / ((mean.hmc + mean.vb) / 2),
    abs_relative_discrepancy = abs(relative_discrepancy),
    hpd.width.95.vb = abs(`2.5%.vb` - `97.5%.vb`),
    hpd.width.95.hmc = abs(`2.5%.hmc` - `97.5%.hmc`),
    hpd.width.99.vb = abs(`0.5%.vb` - `99.5%.vb`),
    hpd.width.99.hmc = abs(`0.5%.hmc` - `99.5%.hmc`),
    discrepancy_95hpdi_width = hpd.width.95.hmc - hpd.width.95.vb,
    discrepancy_99hpdi_width = hpd.width.99.hmc - hpd.width.99.vb,
    abs_discrepancy_99hpdi_width = abs(discrepancy_99hpdi_width),
    relative_discrepancy_99hpdi_width = discrepancy_99hpdi_width / ((hpd.width.99.hmc + hpd.width.99.vb) / 2)
)
cmdf <- cmdf[cmdf$discrepancy < 10, ]
cmdf <- cmdf[cmdf$n_eff.hmc > minEff & cmdf$Rhat < maxRhat, ]
## from PSIS paper, arxiv 1507.02646
# cmdf <- cmdf[cmdf$khat < 0.7, ]

################################################################################
## supplementary plots of discrepancy versus basic diagnostic vars
################################################################################

diagname <- function(x) {
    c(
        "gene" = "Gene",
        "Rhat" = TeX("\\hat{R}"),
        "khat" = TeX("\\hat{k}"),
        "converged" = "ADVI converged",
        "niter" = "ADVI iterations",
        "n_eff.hmc" = "Effective sample size",
        "time.hmc" = "Time taken",
        "n_tot" = "Sample size",
        "p_het" = "Proportion heterozygous",
        "se_mean.hmc" = "SE(mean)",
        "sd.hmc" = "Posterior SD",
        "mean_count" = "mean(RNAseq counts)",
        "sd_count" = "SD(RNAseq counts)",
        "n_wt" = "Number of mut individuals",
        "n_het" = "Number of het individuals",
        "n_hom" = "Number of hom ref individuals"
    )[[x]]
}
typename <- function(x) {
    c(
        "abs_discrepancy" = "Absolute discrepancy in mean",
        "abs_discrepancy_99hpdi_width" = "Absolute discrepancy in\nHPD interval width"
    )[[x]]
}

cmdf <- cmdf %>%
    arrange(-abs(discrepancy))
cmdf$top50 <- c(rep(TRUE, 50), rep(FALSE, nrow(cmdf) - 50))

cmdf <- cmdf %>%
    mutate(
        disc_sc_sdh = discrepancy / sd.hmc,
        disc_sc_sdv = discrepancy / sd.vb,
        disc_sc_meanh = discrepancy / mean.hmc,
        disc_sc_meanv = discrepancy / mean.vb
    )

r <- range(c(cmdf$mean.hmc, cmdf$mean.vb))
## not se mean but other hmc se
diag_vars <- c(
    "gene",
    "Rhat",
    "khat",
    # "converged",
    "niter",
    "n_eff.hmc", "time.hmc", "n_ase",
    "se_mean.hmc", "sd.hmc",
    "n_tot",
    "p_het",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)
# cmdf$converged <- factor(ifelse(as.logical(cmdf$converged), TRUE, FALSE))
for (type in c("abs_discrepancy", "abs_discrepancy_99hpdi_width")) {
    # for (type in c("discrepancy", "disc_sc_sdh", "disc_sc_sdh", "disc_sc_meanh", "disc_sc_meanv")) {
    for (x in diag_vars) {
        geom <- if (is.numeric(cmdf[[x]])) {
            geom_point(shape = 16, size = 0.8, alpha = 0.6
                # , aes(colour = top50)
            )
        } else if (is.logical(cmdf[[x]]) || is.character(cmdf[[x]]) || is.factor(cmdf[[x]])) {
            geom_boxplot(fill = "grey80")
        }
        g <- ggplot(arrange(cmdf, abs(discrepancy))) +
            aes(.data[[x]], .data[[type]]) +
            geom +
            # geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
            # scale_colour_brewer(palette = "Set1", name = "Top 50", direction = -1) +
            # scale_y_log10() +
            theme(legend.position = "bottom") +
            labs(x = diagname(x), y = typename(type))

        ggsave(
            sprintf("%s/%s/diag/%s_%s_%s.png", fpath, model, type, gsub("\\.", "_", x), method),
            width = 5, height = 3
        )
    }
}

sub_diag_vars <- c("niter", "n_eff.hmc")
g1 <- ggplot(arrange(cmdf, abs(discrepancy))) +
    aes(.data[["niter"]], .data[["abs_discrepancy"]]) +
    geom_pointdensity(shape = 16, size = 0.8, alpha = 0.6) +
    scale_colour_viridis(guide = "none") +
    theme(legend.position = "bottom") +
    labs(x = diagname("niter"), y = "Absolute discrepancy in mean")
g2 <- ggplot(arrange(cmdf, abs(discrepancy))) +
    aes(.data[["n_eff.hmc"]], .data[["abs_discrepancy"]]) +
    geom_pointdensity(shape = 16, size = 0.8, alpha = 0.6) +
    theme(legend.position = "bottom") +
    scale_colour_viridis(guide = "none") +
    labs(x = diagname("n_eff.hmc"), y = "Absolute discrepancy in mean")
## replace this with estimate of the number of mixture components, boxplot
g3 <- ggplot(arrange(cmdf, abs(discrepancy))) +
    aes(.data[["sd.hmc"]], .data[["abs_discrepancy"]]) +
    geom_pointdensity(shape = 16, size = 0.8, alpha = 0.6) +
    theme(legend.position = "bottom") +
    scale_colour_viridis(guide = "none") +
    labs(x = "Posterior SD (HMC)", y = "Absolute discrepancy in mean")
g4 <- ggplot(arrange(cmdf, abs(discrepancy))) +
    aes(.data[["khat"]], .data[["abs_discrepancy"]]) +
    geom_pointdensity(shape = 16, size = 0.8, alpha = 0.6) +
    theme(legend.position = "bottom") +
    scale_colour_viridis(guide = "none") +
    labs(x = diagname("khat"), y = "Absolute discrepancy in mean")
gg <- cowplot::plot_grid(g1, g2, labels = if (model == "GT") c("A", "B") else c("C", "D"))
ggsave(
    sprintf("%s/%s/diag/disc_niter_eff.pdf", fpath, model),
    width = 5.5, height = 3.5
)
