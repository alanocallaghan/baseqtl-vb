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

maxRhat <- 1.1 ## from baseqtl-paper repo
minEff <- 500 ## from stan docs (-ish)

args <- parser$parse_args()
tol <- args[["tolerance"]]
model <- args[["model"]]

fpath <- sprintf("fig_%1.0e", tol)
source("src/functions.R")
mkfigdir(fpath, model)

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
df_vb_hmc <- merge(dfs[["vb"]], dfs[["sampling"]], by = by, suffix = c(".vb", ".hmc"))

aucs <- df_vb_hmc |>
    group_by(seed) |>
    summarise(
        auc = roc_auc_vec(
            truth = factor(null.99.hmc),
            estimate = PEP,
            event_level = "second"
        )
    )

if (model == "GT") {
    mean_sd_df <- df_vb_hmc |>
        group_by(gene, snp)
} else {
    mean_sd_df <- df_vb_hmc |>
        group_by(gene, snp, condition)
}
mean_sd_df <- mean_sd_df  |>
    summarise(
        mean_hmc = mean(mean.hmc),
        mean_vb = mean(mean.vb),
        sd_vb = sd(mean.vb)
    )

g <- ggplot(mean_sd_df) +
    aes(x = mean_hmc, y = mean_vb, colour = sd_vb) +
    geom_point() +
    scale_colour_viridis()
# ggsave("tmp.png")
ggsave(sprintf("%s/%s/variability/mean-sd-vb.pdf", fpath, model), width = 5, height = 5)


# g <- ggplot(df_vb_hmc[df_vb_hmc$seed == 42, ]) +
#     aes(x = mean.hmc, y = mean.vb) +
#     geom_point() +
#     scale_colour_viridis()
# ggsave("tmp2.png")
# ggsave(sprintf("%s/%s/variability/mean-sd-disc.png", fpath, model))

disc_df <- df_vb_hmc |>
    mutate(disc = abs(mean.vb - mean.hmc))
if (model == "GT") {
    disc_df <- disc_df |>
        group_by(gene, snp)
} else {
    disc_df <- disc_df |>
        group_by(gene, snp, condition)
}
disc_df <- disc_df |>
    summarise(
        mean_disc = mean(disc),
        sd_disc = sd(disc)
    )

g <- ggplot(disc_df) +
    aes(x = mean_disc, y = sd_disc) +
    geom_point() +
    labs(x = "Mean discrepancy", y = "SD of discrepancy")
ggsave(sprintf("%s/%s/variability/mean-sd-disc.pdf", fpath, model), width = 5, height = 5)
