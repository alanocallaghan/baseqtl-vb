library("rstan")
library("ggplot2")
library("ggdist")
library("ggrastr")
library("distributional")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
library("dplyr")
library("argparse")
library("latex2exp")
library("cowplot")
library("BiocParallel")
rstan_options(auto_write = TRUE)

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
df <- readRDS(infile)

init_types <- c("random", "fixed", "old")
rerun_files <- sprintf(
    "rds/%s_discrepancy_%s_rerun_%s.rds",
    model, tol_str, init_types
)


dfs <- lapply(rerun_files, readRDS)
df_all <- do.call(rbind, dfs)
rerun_df <- readRDS(sprintf("rds/%s_rerun_discrepant_full.rds", model))

df_all$rn <- rownames(df_all)
df_all$rowidx <- gsub(
    "(Psoriasis_skin|normal_skin)?\\.?(phi|theta|bj|betas\\[1\\]|betas\\[2\\]|rai0\\[\\d+\\]|lp__)(\\d*$)",
    "\\3",
    df_all$rn
)
df_all$rowidx <- as.numeric(df_all$rowidx)
df_all <- df_all[df_all$rowidx != "", ]
df_all$rowidx[is.na(df_all$rowidx)] <- 0
df_all <- df_all[!is.na(df_all$param), ]

df_new <- df_all[grep("bj", df_all$param), ]

by <- if (model == "GT") {
    c("gene", "snp")
} else {
    c("gene", "snp", "condition")
}

cols <- c("mean", "khat", "time", "gene", "snp", "init", "rowidx")
if (model == "noGT") cols <- c(cols, "condition")

dfn <- df_new[, cols]
dfn[, c("mean", "khat", "time")] <- lapply(
    dfn[, c("mean", "khat", "time")],
    as.numeric
)
dfm <- merge(
    dfn, df,
    by = by,
    suffixes = c(".new", ".old")
)

dfm <- dfm %>% mutate(
    new_discrepancy = hmc - mean,
    discrepancy = discrepancy
)

df_sum <- dfm %>%
    group_by(gene, snp, init) %>%
    summarise(
        sd_mean = sd(mean),
        discrepancy = unique(discrepancy),
        sd_disc = sd(new_discrepancy),
        mean_disc = mean(new_discrepancy),
        mean_abs_disc = mean(abs(new_discrepancy)),
        sd_abs_disc = sd(abs(new_discrepancy)),
        .groups = "drop"
    )
g <- ggplot(df_sum) +
    aes(init, sd_abs_disc) +
    geom_violin() +
    geom_jitter(size = 0.5, width = 0.25, height = 0) +
    labs(x = "Initialisation type", y = "sd(abs(discrepancy))")
ggsave(
    sprintf("%s/%s/diag/init_sd_disc.pdf", fpath, model),
    width = 2.5, height = 2.5
)

g <- ggplot(df_sum) +
    aes(init, mean_disc) +
    geom_violin() +
    geom_jitter(size = 0.5, width = 0.25, height = 0) +
    labs(x = "Initialisation type", y = "mean(discrepancy)")
ggsave(
    sprintf("%s/%s/diag/init_mean_disc.png", fpath, model),
    width = 5, height = 5
)

g <- ggplot(df_sum) +
    aes(init, mean_abs_disc) +
    geom_violin() +
    geom_jitter(size = 0.5, width = 0.25, height = 0) +
    labs(x = "Initialisation type", y = "mean(abs(discrepancy))")
ggsave(
    sprintf("%s/%s/diag/init_mean_abs_disc.pdf", fpath, model),
    width = 2.5, height = 2.5
)

g <- ggplot(df_sum) +
    aes(init, sd_mean) +
    geom_violin() +
    geom_jitter(size = 0.5, width = 0.25, height = 0) +
    labs(x = "Initialisation type", y = "sd(estimate)")
ggsave(sprintf("%s/%s/diag/init_se.png", fpath, model), width = 5, height = 5)


## posterior samples HMC vs VB, split by chain
full_posteriors_vb <- readRDS(
    sprintf("rds/%s_full_posterior_discrepant_vb.rds", model)
)
full_posteriors_hmc <- readRDS(
    sprintf("rds/%s_full_posterior_discrepant_hmc.rds", model)
)

## if noGT, full_posteriors_hmc[[i]] and full_posteriors_vb[[i]] are lists
if (model == "noGT") {
    dfs_slabplot <- lapply(
        seq_len(length(full_posteriors_hmc)),
        function(i) {
            dfh <- data.frame(
                method = "hmc",
                association = as.character(i),
                chain = as.character(rep(1:4, each = 1000))
            )
            dfh[, names(full_posteriors_hmc[[i]])] <- sapply(
                full_posteriors_hmc[[i]], function(x) extract(x, par = "bj")$bj
            )
            dfv <- data.frame(
                method = "vb",
                association = as.character(i),
                chain = as.character(rep(1:4, each = 1000))
            )
            dfv[, names(full_posteriors_vb[[i]])] <- sapply(
                full_posteriors_vb[[i]], function(x) extract(x, par = "bj")$bj
            )
            df <- rbind(dfh, dfv)
            df <- reshape2::melt(df, id.vars = c("method", "association", "chain"))
            colnames(df)[4:5] <- c("condition", "bj")
            df
        }
    )
} else {
    dfs_slabplot <- lapply(1:length(full_posteriors_hmc), function(i) {
        rbind(
            data.frame(
                method = "hmc",
                bj = extract(full_posteriors_hmc[[i]], par = "bj")$bj,
                association = as.character(i),
                chain = rep(1:4, each = 1000)
            ),
            data.frame(
                method = "vb",
                bj = extract(full_posteriors_vb[[i]], par = "bj")$bj,
                association = as.character(i),
                chain = "vb"
            )
        )
    })
}

df_slabplot <- do.call(rbind, dfs_slabplot)
df_slabplot <- df_slabplot %>% filter(association %in% c(1:15))

df_slabplot$association_condition <- paste(df_slabplot$association, df_slabplot$condition)
plot_aes <- if (model == "noGT") {
    aes(x = association_condition, y = bj, group = method, colour = chain, order = chain)
} else {
    aes(x = association, y = bj, group = method, colour = chain, order = chain)
}

g <- ggplot(df_slabplot) +
    plot_aes +
    geom_dots(
        position = position_dodge(),
        # height = 2, dotsize = 0.5, stackratio = 1.1, scale = 1.5
    ) +
    labs(x = "Association", y = TeX("$\\beta_{aFC}$")) +
    scale_fill_manual(
        name = "Chain",
        aesthetics = c("fill", "colour"),
        values = setNames(
            c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "grey80"),
            c(1:4, "vb")
        )
    ) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_blank()
    )
# ggsave(sprintf("%s/%s/diag/full-posteriors-dots.png", fpath, model), width = 20, height = 5)

plot_aes <- if (model == "noGT") {
    aes(x = association_condition, y = bj, colour = chain, fill = chain)
} else {
    aes(x = association, y = bj, colour = chain, fill = chain)
}
g <- ggplot(df_slabplot) +
    plot_aes +
    stat_slab(position = position_dodge(width = 0.75), slab_size = 0.3, fill = "grey90") +
    labs(x = "Association", y = TeX("$\\beta_{aFC}$")) +
    scale_fill_manual(
        name = "Chain",
        aesthetics = c("fill", "colour"),
        values = setNames(
            c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "grey50"),
            c(1:4, "vb")
        )
    ) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_blank()
    )
ggsave(sprintf("%s/%s/diag/full-posteriors-slab.pdf", fpath, model), width = 20, height = 5)





################################################################################
## Unused
################################################################################


# # maybe split into pointrange by init type?
# g <- ggplot(dfm) +
#     aes(abs(discrepancy), abs(new_discrepancy), colour = init) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_jitter(size = 0.5, width = 0.001, height = 0) +
#     labs(x = "Original discrepancy", y = "New discrepancy") +
#     scale_colour_brewer(palette = "Set2", name = "Initialisation")
# ggsave(
#     sprintf("%s/%s/diag/rerun_discrepant.png", fpath, model),
#     width = 5, height = 5
# )


# g <- ggplot(dfm) +
#     aes(init, abs(new_discrepancy)) +
#     geom_violin(alpha = 0) +
#     geom_segment(
#         aes(
#             x = 0.5, xend = 5,
#             y = abs(discrepancy), yend = abs(discrepancy),
#             colour = gene
#         ),
#         linetype = "dashed", alpha = 0.7, linewidth = 0.4
#     ) +
#     geom_violin() +
#     geom_boxplot(outlier.colour = NA, width = 0.25) +
#     geom_jitter(aes(colour = gene), width = 0.25, alpha = 0.3, size = 0.3) +
#     scale_colour_discrete(name = NULL, guide = "none") +
#     labs(x = "Initialisation method", y = "abs(discrepancy)")
# ggsave(
#     sprintf("%s/%s/diag/init_discrepant.png", fpath, model),
#     width = 5, height = 5
# )


# g <- ggplot(dfm) +
#     aes(init, abs(new_discrepancy)) +
#     geom_violin(alpha = 0) +
#     geom_segment(
#         aes(
#             x = 0.5, xend = 5,
#             y = abs(discrepancy), yend = abs(discrepancy),
#             colour = gene
#         ),
#         linetype = "dashed", alpha = 0.2, linewidth = 0.4
#     ) +
#     geom_violin() +
#     geom_boxplot(outlier.colour = NA, width = 0.25) +
#     geom_jitter(aes(colour = gene), width = 0.25, , alpha = 0.3, size = 0.3) +
#     scale_y_log10() +
#     scale_colour_discrete(name = NULL, guide = "none") +
#     labs(x = "Initialisation method", y = "abs(discrepancy)")
# ggsave(
#     sprintf("%s/%s/diag/init_discrepant_log.png", fpath, model),
#     width = 5, height = 5
# )



# md_df <- dfm %>%
#     group_by(gene, snp) %>%
#     summarise(
#         new_mabs_disc = mean(abs(new_discrepancy)),
#         mabs_disc = mean(abs(discrepancy)),
#         .groups = "drop"
#     )

# # todo: make this y pointrange
# g <- ggplot(md_df) +
#     aes(mabs_disc, new_mabs_disc) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_point() +
#     labs(x = "New mean absolute discrepancy", y = "Original absolute discrepancy")
# ggsave(
#     sprintf("%s/%s/diag/comp-discrepancy.png", fpath, model),
#     width = 5, height = 5
# )

