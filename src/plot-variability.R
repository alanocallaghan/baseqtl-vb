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



lab_str <- "HMC: %s\nADVI: %s"
df_vb_hmc$nullstr99 <- sprintf(lab_str,
    ifelse(df_vb_hmc$null.99.hmc, "yes", "no"),
    ifelse(df_vb_hmc$null.99.vb, "yes", "no")
)

null_levs <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "no",  "yes", "yes")
)
null_ord <- null_levs[c(1, 3, 2, 4)]


df_summ <- df_vb_hmc |>
    group_by(gene, snp) |>
    summarise(
        minvb = min(mean.vb), maxvb = max(mean.vb), meanvb = mean(mean.vb),
        meanhmc = mean(mean.hmc), minhmc = min(`2.5%.hmc`), maxhmc = min(`97.5%.hmc`),
        sig = unique(nullstr99),
        .groups = "drop_last"
    )
lims <- range(as.matrix(df_summ[, 3:8]))

g <- ggplot(df_summ) +
    aes(x = meanhmc, y = meanvb) +
    geom_pointrange(aes(ymin = minvb, ymax = maxvb), fatten = 1, alpha = 0.25) +
    # geom_pointrange(aes(xmin = minhmc, xmax = maxhmc), fatten = 1, alpha = 0.25) +
    facet_wrap(~sig) +
    scale_x_continuous(name = "HMC estimate", limits = lims) +
    scale_y_continuous(name = "ADVI estimate", limits = lims)
ggsave(sprintf("%s/%s/variability/mean-var-pointrange.pdf", fpath, model), width = 12, height = 10)
# ggsave("tmp.png")


if (model == "GT") {
    n_vb <- dfs$vb |>
        group_by(seed) |>
        summarise(
            n = n_distinct(test),
            .groups = "drop_last"
        )
    n_samp <- dfs$sampling |>
        summarise(
            n = n_distinct(test),
            .groups = "drop_last"
        )
} else {
    n_vb <- dfs$vb |>
        group_by(seed, condition) |>
        summarise(
            n = n_distinct(test),
            .groups = "drop_last"
        )
    n_samp <- dfs$sampling |>
        group_by(condition) |>
        summarise(
            n = n_distinct(test),
            .groups = "drop_last"
        )

}


aucs <- df_vb_hmc |>
    group_by(seed) |>
    summarise(
        auc = roc_auc_vec(
            truth = factor(null.99.hmc),
            estimate = PEP,
            event_level = "second"
        ),
        .groups = "drop_last"
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
        sd_vb = sd(mean.vb),
        .groups = "drop_last"
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
        sd_disc = sd(disc),
        .groups = "drop_last"
    )

g <- ggplot(disc_df) +
    aes(x = mean_disc, y = sd_disc) +
    geom_point() +
    labs(x = "Mean discrepancy", y = "SD of discrepancy")
ggsave(sprintf("%s/%s/variability/mean-sd-disc.pdf", fpath, model), width = 5, height = 5)




## mean/sd of discrepancy across seeds
## mabs???/sd of discrepancy across seeds
## to decide which seed to use

# disc_seed_df <- df_vb_hmc |>
#     mutate(disc = abs(mean.vb - mean.hmc))
# if (model == "GT") {
#     disc_seed_df <- disc_seed_df |>
#         group_by(seed)
# } else {
#     disc_seed_df <- disc_seed_df |>
#         group_by(seed)
# }
# disc_seed_df <- disc_seed_df |> summarise(
#         sabs = sum(abs(disc)),
#         mabs = mean(abs(disc)),
#         medabs = median(abs(disc)),
#         minabs = min(abs(disc)),
#         maxabs = max(abs(disc)),
#         sdabs = sd(abs(disc)),
#         sd = sd(disc),
        # .groups = "drop_last"
#     )
