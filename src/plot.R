library("argparse")
library("ggplot2")
library("baseqtl")
library("ggpointdensity")
library("ggExtra")
library("viridis")
library("dplyr")
theme_set(theme_bw())

parser <- ArgumentParser()
parser$add_argument(
    "-i", "--inference",
    default = "sampling",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-2,
    type = "double"
)
parser$add_argument(
    "-m", "--model",
    default = "GT",
    type = "character"
)

args <- parser$parse_args()
method <- args[["inference"]]
tol <- args[["tolerance"]]
model <- args[["model"]]

mname <- switch(method,
    "vb" = "ADVI",
    "optimizing" = "MAP",
    "sampling" = "HMC"
)

# mtol <- if (method == "vb") sprintf("vb_%1.0e", tol) else method
mtol <- method

maxRhat <- 1.1
minEff <- 500

sfile <- sprintf("rds/%s/sfile.rds", model)
combfile <- sprintf("rds/%s/%s_combined.rds", model, method)

approx_res_df <- readRDS(combfile)
sample_res_df <- readRDS(sfile)

if (model == "GT") {
    by <- c("gene", "snp", "test")
} else {
    by <- c("gene", "snp", "test", "condition")
}

sample_res_df <- sample_res_df[sample_res_df$n_eff > minEff & sample_res_df$Rhat < maxRhat, ]
if (method == "sampling") {
    approx_res_df <- approx_res_df[approx_res_df$n_eff > minEff & approx_res_df$Rhat < maxRhat, ]
}

approx_res_df$test <- paste(approx_res_df$gene, approx_res_df$snp, sep = "_")
sample_res_df$test <- paste(sample_res_df$Gene_id, sample_res_df$tag , sep = "_")
sample_res_df$gene <- sample_res_df$Gene_id
sample_res_df$snp <- sample_res_df$tag
sample_res_df$Gene_id <- NULL
sample_res_df$tag <- NULL


approx_res_df$null.95 <- sign(approx_res_df$"2.5%") == sign(approx_res_df$"97.5%")
approx_res_df$null.50 <- sign(approx_res_df$"25.0%") == sign(approx_res_df$"75.0%")
sample_res_df$null.95 <- sign(sample_res_df$"log2_aFC_2.5.") == sign(sample_res_df$"log2_aFC_97.5.")
sample_res_df$null.50 <- sign(sample_res_df$"log2_aFC_25.") == sign(sample_res_df$"log2_aFC_75.")

sample_res_df <- sample_res_df[sample_res_df$test %in% approx_res_df$test, ]
im <- match(sample_res_df$test, approx_res_df$test)
approx_res_df <- approx_res_df[im, ]
stopifnot(all(sample_res_df$test == approx_res_df$test))

r <- range(c(sample_res_df$log2_aFC_mean, approx_res_df$mean))

mdf <- merge(approx_res_df, sample_res_df, by = by, suffix = c(".VB", ".HMC"))
mdf$null.99.VB <- ifelse(mdf$null.99.VB, "no", "yes")
mdf$null.95.VB <- ifelse(mdf$null.95.VB, "no", "yes")
mdf$null.50.VB <- ifelse(mdf$null.50.VB, "no", "yes")
mdf$null.95.HMC <- ifelse(mdf$null.95.HMC, "no", "yes")
mdf$null.50.HMC <- ifelse(mdf$null.50.HMC, "no", "yes")

lab_str <- paste0(mname, ": %s\nHMC: %s\n")
levs <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "no",  "yes", "yes")
)
# levs <- rev(levs)
flev <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "yes",  "no", "yes")
)

mdf$Null.99 <- sprintf(lab_str, mdf$null.99.VB, mdf$null.99.HMC)
mdf$Null.95 <- sprintf(lab_str, mdf$null.95.VB, mdf$null.95.HMC)
mdf$Null.50 <- sprintf(lab_str, mdf$null.50.VB, mdf$null.50.HMC)
mdf$Null.95.50 <- sprintf(lab_str, mdf$null.50.VB, mdf$null.95.HMC)

mdf$Null.99 <- factor(mdf$Null.99, levels = flev)
mdf$Null.95 <- factor(mdf$Null.95, levels = flev)
mdf$Null.50 <- factor(mdf$Null.50, levels = flev)
mdf$Null.95.50 <- factor(mdf$Null.95.50, levels = flev)

mdf$discrepancy <- mdf$mean - mdf$log2_aFC_mean

## not se mean but other hmc se
diag_vars <- c(
    "gene",
    # "Rhat",
    # "n_eff",
    "time", "n_ase",
    "log2_aFC_se_mean", "log2_aFC_sd",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)

mdf2 <- mdf[mdf$Rhat.VB < maxRhat & mdf$Rhat.HMC < maxRhat & mdf$n_eff.VB > minEff & mdf$n_eff.HMC > minEff, ]
mdf2 <- mdf2[mdf2$mean < 2, ]

g <- ggplot(mdf2) +
    aes_string("log2_aFC_mean", "mean", colour = "log2_aFC_se_mean") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_viridis(name = "SE (x)") +
    lims(x = r, y = r) +
    labs(x = "Old estimate", y = "New estimate")
ggsave(sprintf("fig/%s/estimates/se_x.png", model), width = 7, height = 7)
g <- ggplot(mdf2) +
    aes_string("log2_aFC_mean", "mean", colour = "se_mean") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_viridis(name = "SE (y)") +
    lims(x = r, y = r) +
    labs(x = "Old estimate", y = "New estimate")
ggsave(sprintf("fig/%s/estimates/se_y.png", model), width = 7, height = 7)

for (x in diag_vars) {
    scale <- if (x == "gene") scale_colour_discrete(guide="none") else scale_colour_viridis()

    g <- ggplot(mdf) +
        aes_string(x, "abs(discrepancy)") +
        geom_point(size = 0.8) +
        # geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
        labs(x = x, y = "Discrepancy")

    ggsave(sprintf("fig/%s/diag/disc_%s_%s.png", model, x, mtol), width = 7, height = 7)
}

mdf <- mdf[order(mdf$Null.99), ]

g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.99) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/%s/estimates/%s_mcmc_all_categorical_99.png", model, mtol), width = 7, height = 7)

mdf <- mdf[order(mdf$Null.95), ]
g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.95) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/%s/estimates/%s_mcmc_all_categorical_95.png", model, mtol), width = 7, height = 7)

mdf <- mdf[order(mdf$Null.95.50), ]
g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.95.50) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/%s/estimates/%s_mcmc_all_categorical_50.png", model, mtol), width = 7, height = 7)
