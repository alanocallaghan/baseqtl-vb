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
    default = "GT",
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
dfs <- lapply(
    methods,
    function(method) {
        cat(method, "\n")
        mtol <- sprintf("%s_%1.0e", method, tol)
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
    group_by(seed.vb) |>
    summarise(
        auc = roc_auc_vec(
            truth = factor(null.99.hmc),
            estimate = PEP.vb,
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
mean_sd_df <- mean_sd_df |>
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


# g <- ggplot(df_vb_hmc[df_vb_hmc$seed.vb == 42, ]) +
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
    geom_pointdensity() +
    scale_colour_viridis(guide = "none") +
    labs(x = "Mean absolute difference", y = "SD of absolute difference")
g <- cowplot::plot_grid(g, labels = if (model == "GT") "A" else "B")
ggsave(sprintf("%s/%s/variability/mean-sd-disc.pdf", fpath, model), width = 3.5, height = 3.5)



perf_dfs <- lapply(
    unique(df_vb_hmc$seed.vb),
    function(seed) {
        df <- df_vb_hmc[df_vb_hmc$seed.vb == seed, ]
        prob <- df[["PEP.vb"]] %||% df[["PEP"]]
        pred_vb <- prediction(
            predictions = prob,
            labels = factor(df[["null.99.hmc"]], levels = c(TRUE, FALSE))
        )

        perf_auroc_vb <- performance(pred_vb, "auc")@y.values[[1]]
        perf_aupr_vb <- performance(pred_vb, "aucpr")@y.values[[1]]
        perf_pr_vb <- performance(pred_vb, "prec", "rec")
        perf_roc_vb <- performance(pred_vb, "tpr", "fpr")
        data.frame(
            seed = seed,
            auroc = perf_auroc_vb,
            aupr = perf_aupr_vb,
            pr_x = perf_pr_vb@x.values[[1]],
            pr_y = perf_pr_vb@y.values[[1]],
            roc_x = perf_roc_vb@x.values[[1]],
            roc_y = perf_roc_vb@y.values[[1]]
        )
    }
)
perf_df <- do.call(rbind, perf_dfs)

g1 <- ggplot(perf_df) +
    aes(roc_x, roc_y, colour = factor(seed)) +
    geom_path() +
    lims(x = 0:1, y = 0:1) +
    labs(x = "True positive rate", y = "False positive rate") +
    scale_colour_brewer(palette = "Accent", name = NULL, guide = "none") +
    theme(legend.position = "bottom") +
    annotate(
        label = paste0("AUC = ", paste(round(range(perf_df$auroc), digits = 3), collapse = "-")),
        geom = "text",
        x = 0.6,
        y = 0.5
    )
ggsave(sprintf("%s/%s/variability/roc.pdf", fpath, model), width = 3.5, height = 3.5)


time_dfs <- lapply(
    unique(df_vb_hmc$seed.vb),
    function(seed) {
        df <- df_vb_hmc[df_vb_hmc$seed.vb == seed, ]
        df <- df[, !grepl("%", colnames(df))]
        bc <- "null.99.hmc"
        peps <- seq(0, 1, by = 0.01)
        peps <- peps[peps != 1]
        sens_vb <- sapply(
            peps,
            function(x) {
                sens_vec(
                    truth = factor(df[[bc]], levels = c(TRUE, FALSE)),
                    estimate = factor(df$PEP.vb > x, levels = c(TRUE, FALSE))
                )
            }
        )
        tabs <- sapply(
            peps,
            function(x) {
                mean(df$PEP.vb > x)
            }
        )

        spec_vb <- sapply(
            peps,
            function(x) {
                spec_vec(
                    truth = factor(df[[bc]], levels = c(TRUE, FALSE)),
                    estimate = factor(df$PEP.vb > x, levels = c(TRUE, FALSE))
                )
            }
        )
        time_vb <- sapply(
            peps,
            function(x) {
                sum(df$time.hmc[df$PEP.vb > x]) + sum(df$time.vb)
            }
        )
        data.frame(seed = seed, pep = peps, sens = sens_vb, spec = spec_vb, time = time_vb)
    }
)
time_df <- do.call(rbind, time_dfs)

hmctimes <- dfs$sampling$time


g2 <- ggplot(time_df) +
    aes(sens, time / 3600, colour = factor(seed)) +
    geom_path() +
    scale_colour_brewer(palette = "Accent", name = NULL, guide = "none") +
    geom_hline(
        yintercept = sum(hmctimes) / 3600,
        linetype = "dashed"
    ) +
    annotate(
        geom = "text",
        x = mean(range(time_df$sens)),
        y = sum(hmctimes) / 3600,
        label = "Total time\nwithout screening",
        vjust = -0.2,
    ) +
    ylim(0, max(sum(hmctimes), time_df$time)) +
    # scale_y_log10() +
    labs(x = "Sensitivity", y = "Total time (hr)")
ggsave(sprintf("%s/%s/variability/time-pep.pdf", fpath, model), width = 3.5, height = 3.5)
# ggsave("tmp.pdf")

plot_grid(g1, g2, labels = "AUTO")
ggsave(sprintf("%s/%s/variability/time-roc-pep.pdf", fpath, model), width = 5.5, height = 3)

time_dfs_levs <- lapply(
    unique(df_vb_hmc$seed.vb),
    function(seed) {
        df <- df_vb_hmc[df_vb_hmc$seed.vb == seed, ]
        levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
        df <- add_nulls(df, suffix = ".hmc")
        df <- add_nulls(df, suffix = ".vb")
        bc <- "null.99.hmc"
        sens_vb <- sapply(
            levs,
            function(x) {
                column <- sprintf("null.%s.vb", x)
                sens_vec(
                    truth = factor(df[[bc]], levels = c(TRUE, FALSE)),
                    estimate = factor(df[[column]], levels = c(TRUE, FALSE))
                )
            }
        )
        spec_vb <- sapply(
            levs,
            function(x) {
                column <- sprintf("null.%s.vb", x)
                spec_vec(
                    truth = factor(df[[bc]], levels = c(TRUE, FALSE)),
                    estimate = factor(df[[column]], levels = c(TRUE, FALSE))
                )
            }
        )
        time_vb <- sapply(
            levs,
            function(x) {
                column <- sprintf("null.%s.vb", x)
                sum(df$time.hmc[df[[column]]]) + sum(df$time.vb)
            }
        )

        data.frame(seed = seed, sens = sens_vb, spec = spec_vb, time = time_vb)
    }
)
time_df_levs <- do.call(rbind, time_dfs_levs)


g <- ggplot(time_df_levs) +
    aes(sens, time / 3600, colour = factor(seed)) +
    geom_path() +
    scale_colour_brewer(palette = "Accent", name = NULL, guide = "none") +
    geom_hline(
        yintercept = sum(hmctimes) / 3600,
        linetype = "dashed"
    ) +
    annotate(
        geom = "text",
        x = mean(range(time_df_levs$sens)),
        y = sum(hmctimes) / 3600,
        label = "Total time\nwithout screening",
        vjust = -0.2,
    ) +
    ylim(0, max(sum(hmctimes), time_df$time) / 3600) +
    # scale_y_log10() +
    labs(x = "Sensitivity", y = "Total time (hr)")
ggsave(sprintf("%s/%s/variability/time-levs.pdf", fpath, model), width = 3.5, height = 3.5)

# ggsave("tmp2.png")


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
