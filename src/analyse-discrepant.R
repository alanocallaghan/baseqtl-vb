library("rstan")
library("ggplot2")
library("ggdist")
library("distributional")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
library("dplyr")
library("argparse")

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
parser$add_argument(
    "-i", "--inference",
    default = "vb",
    type = "character"
)

args <- parser$parse_args()
source("src/functions.R")
theme_set(theme_bw())

tol <- args[["tolerance"]]
model <- args[["model"]]
method <- args[["inference"]]

tol_str <- if (method == "vb") sprintf("%s_%1.0e", method, tol) else method
fun <- match.fun(method)

file <- sprintf("rds/%s_discrepancies_%s.rds", model, tol_str)
df <- readRDS(file)

fpath <- sprintf("fig_%1.0e", tol)
dir.create(fpath, showWarnings = FALSE, recursive = TRUE)
# fpath <- "fig"
n_replicates <- 20


init_types <- c("random", "fixed", "old")
# init_types <- c("random", "opt", "fixed", "old")
rerun_files <- sprintf("rds/%s_discrepancy_%s_rerun_%s.rds", model, tol_str, init_types)

fit_fun <- match.fun(paste("fit_stan", model, sep = "_"))



if (!all(file.exists(rerun_files))) {
        
    fr_res <- parallel::mclapply(1:nrow(df), function(i) {
        gene_data <- get_gene_data(df[i, "gene"], model)
        out <- replicate(
            n_replicates,
            {
                d <- fit_fun(
                    gene_data,
                    df[i, "snp"],
                    gene = df[i, "gene"],
                    init = "random",
                    algorithm = "fullrank",
                    vars = NULL
                )
                # ds <- d[d$condition == df[i, "condition"], ]
                # abs(df[i, "hmc"]) - ds[grep("bj", rownames(ds)), "mean"]
                d$param <- gsub("1$", "", rownames(d))
                d
            },
            simplify = FALSE
        )
        do.call(rbind, out)
    }, mc.cores = 8)
    fr_df <- do.call(rbind, fr_res)
    fr_df$init <- "random"
    fr_df$algorithm <- "fullrank"
    saveRDS(fr_df, sprintf("rds/%s_discrepant_%s_fr.rds", model, tol_str))

    rand_res <- parallel::mclapply(1:nrow(df), function(i) {
        gene_data <- get_gene_data(df[i, "gene"], model)
        out <- replicate(
            n_replicates,
            {
                d <- fit_fun(
                    gene_data,
                    df[i, "snp"],
                    gene = df[i, "gene"],
                    init = "random",
                    vars = NULL
                )
                # ds <- d[d$condition == df[i, "condition"], ]
                # abs(df[i, "hmc"]) - ds[grep("bj", rownames(ds)), "mean"]
                d$param <- gsub("1$", "", rownames(d))
                d
            },
            simplify = FALSE
        )
        do.call(rbind, out)
    }, mc.cores = 8)
    rand_df <- do.call(rbind, rand_res)
    rand_df$init <- "random"
    saveRDS(rand_df, rerun_files[[1]])

    fixed_res <- parallel::mclapply(1:nrow(df), function(i) {
        gene_data <- get_gene_data(df[i, "gene"], model)
        out <- replicate(
            n_replicates,
            {
                d <- fit_fun(gene_data, df[i, "snp"], gene = df[i, "gene"], init = "fixed", vars = NULL)
                d$param <- gsub("1$", "", rownames(d))
                d
            },
            simplify = FALSE
        )
        do.call(rbind, out)
    }, mc.cores = 8)
    fixed_df <- do.call(rbind, fixed_res)
    fixed_df$init <- "fixed"
    saveRDS(fixed_df, rerun_files[[3]])

    old_res <- parallel::mclapply(1:nrow(df), function(i) {
        gene_data <- get_gene_data(df[i, "gene"], model)
        init <- list(
            betas = c(5, 0),
            bj = df[i, "vb"],
            phi = 10,
            theta = 10
        )
        out <- replicate(
            n_replicates,
            {
                d <- fit_fun(gene_data, df[i, "snp"], gene = df[i, "gene"], init = init, vars = NULL)
                d$param <- gsub("1$", "", rownames(d))
                d
            },
            simplify = FALSE
        )
        do.call(rbind, out)
    }, mc.cores = 8)
    old_df <- do.call(rbind, old_res)
    old_df$init <- "old"
    saveRDS(old_df, rerun_files[[4]])

    dfs <- list(
        fixed_df = fixed_df,
        # opt_df = opt_df,
        rand_df = rand_df,
        old_df = old_df
    )

    rerun_hmc <- parallel::mclapply(
        1:nrow(df),
        function(i) {
            gene_data <- get_gene_data(df[i, "gene"], model)
            d <- fit_fun(
                gene_data,
                df[i, "snp"],
                gene = df[i, "gene"],
                method = "sampling",
                cores = 1,
                vars = NULL
            )
            d$param <- gsub("1$", "", rownames(d))
            # ds <- d[d$condition == df[i, "condition"], ]
            # abs(df[i, "hmc"]) - ds[grep("bj", rownames(ds)), "mean"]
            d
        },
        mc.cores = 8
    )
    rerun_df <- do.call(rbind, rerun_hmc)
    saveRDS(rerun_df, sprintf("rds/%s_rerun_discrepant_full.rds", model))

}

dfs <- lapply(rerun_files, readRDS)
df_all <- do.call(rbind, dfs)
rerun_df <- readRDS(sprintf("rds/%s_rerun_discrepant_full.rds", model))

df_all$rn <- rownames(df_all)
df_all$rowidx <- gsub("(phi|theta|bj|betas\\[1\\]|betas\\[2\\]|lp__)(\\d*$)", "\\2", df_all$rn)
df_all$rowidx <- as.numeric(df_all$rowidx)
df_all$rowidx[is.na(df_all$rowidx)] <- 0

df_new <- df_all[df_all$param == "bj", ]

by <- if (model == "GT") {
    c("gene", "snp")
} else c("gene", "snp", "condition")

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
# ggsave(sprintf("%s/%s/diag/init_sd_disc.png", fpath, model), width = 5, height = 5)
ggsave(sprintf("%s/%s/diag/init_sd_disc.pdf", fpath, model), width = 2.5, height = 3.5)

g <- ggplot(df_sum) +
    aes(init, mean_disc) +
    geom_violin() +
    geom_jitter(size = 0.5, width = 0.25, height = 0) +
    labs(x = "Initialisation type", y = "mean(discrepancy)")
ggsave(sprintf("%s/%s/diag/init_mean_disc.png", fpath, model), width = 5, height = 5)

g <- ggplot(df_sum) +
    aes(init, mean_abs_disc) +
    geom_violin() +
    geom_jitter(size = 0.5, width = 0.25, height = 0) +
    labs(x = "Initialisation type", y = "mean(abs(discrepancy))")
# ggsave(sprintf("%s/%s/diag/init_mean_abs_disc.png", fpath, model), width = 3, height = 3)
ggsave(sprintf("%s/%s/diag/init_mean_abs_disc.pdf", fpath, model), width = 2.5, height = 3.5)

g <- ggplot(df_sum) +
    aes(init, sd_mean) +
    geom_violin() +
    geom_jitter(size = 0.5, width = 0.25, height = 0) +
    labs(x = "Initialisation type", y = "sd(estimate)")
ggsave(sprintf("%s/%s/diag/init_se.png", fpath, model), width = 5, height = 5)


# maybe split into pointrange by init type?
g <- ggplot(dfm) +
    aes(abs(discrepancy), abs(new_discrepancy), colour = init) +
    geom_abline(slope = 1, intercept = 0, linetype =  "dashed") +
    geom_jitter(size = 0.5, width = 0.001, height = 0) +
    labs(x = "Original discrepancy", y = "New discrepancy") +
    scale_colour_brewer(palette = "Set2", name = "Initialisation")
ggsave(sprintf("%s/%s/diag/rerun_discrepant.png", fpath, model), width = 5, height = 5)


g <- ggplot(dfm) +
    aes(init, abs(new_discrepancy)) +
    geom_violin(alpha = 0) +
    geom_segment(
        aes(x = 0.5, xend = 5,
            y = abs(discrepancy), yend = abs(discrepancy),
            colour = gene
        ),
        linetype = "dashed", alpha = 0.7, size = 0.4
    ) +
    geom_violin() +
    geom_boxplot(outlier.colour = NA, width = 0.25) +
    geom_jitter(aes(colour = gene), width = 0.25, alpha = 0.3, size = 0.3) +
    scale_colour_discrete(name = NULL, guide = "none") +
    labs(x = "Initialisation method",  y = "abs(discrepancy)")
ggsave(sprintf("%s/%s/diag/init_discrepant.png", fpath, model), width = 5, height = 5)


g <- ggplot(dfm) +
    aes(init, abs(new_discrepancy)) +
    geom_violin(alpha = 0) +
    geom_segment(
        aes(x = 0.5, xend = 5,
            y = abs(discrepancy), yend = abs(discrepancy),
            colour = gene
        ),
        linetype = "dashed", alpha = 0.2, size = 0.4
    ) +
    geom_violin() +
    geom_boxplot(outlier.colour = NA, width = 0.25) +
    geom_jitter(aes(colour = gene), width = 0.25, , alpha = 0.3, size = 0.3) +
    scale_y_log10() +
    scale_colour_discrete(name = NULL, guide = "none") +
    labs(x = "Initialisation method",  y = "abs(discrepancy)")
ggsave(sprintf("%s/%s/diag/init_discrepant_log.png", fpath, model), width = 5, height = 5)



md_df <- dfm %>%
    group_by(gene, snp) %>%
    summarise(
        new_mabs_disc = mean(abs(new_discrepancy)),
        mabs_disc = mean(abs(discrepancy)),
        .groups = "drop"
    )

# todo: make this y pointrange
g <- ggplot(md_df) +
    aes(mabs_disc, new_mabs_disc) +
    geom_abline(slope = 1, intercept = 0, linetype =  "dashed") +
    geom_point() +
    labs(x = "New mean absolute discrepancy", y = "Original absolute discrepancy")
ggsave(sprintf("%s/%s/diag/comp-discrepancy.png", fpath, model), width = 5, height = 5)





fr_df <- readRDS(sprintf("rds/%s_discrepant_%s_fr.rds", model, tol_str))
# fr_df <- do.call(rbind, fr_df) 

cols <- c("mean", "khat", "time", "gene", "snp")
if (model == "noGT") cols <- c(cols, "condition")

dfn <- fr_df[, cols]
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





if (model == "GT") {
    stop()
}


top_newd <- dfm[order(-abs(dfm$new_discrepancy)), ]
top_newd <- top_newd[1:10, ]

cols <- c("gene", "snp", "condition")

inds <- rerun_df[, cols] == top_newd[1, cols, drop=TRUE]
top_newd_hmc <- rerun_df[apply(inds, 1, all), ]
top_newd_hmc <- top_newd_hmc[, !grepl("%", colnames(top_newd_hmc)), ]

retry_advi <- df_all[df_all$rowidx == top_newd$rowidx[[1]], ]
retry_hmc <- head(top_newd_hmc)

# retry_advi$param <- rownames(retry_advi)
retry_advi <- retry_advi[retry_advi$param != "lp__", ]
# retry_hmc$param <- rownames(retry_hmc)
retry_hmc <- retry_hmc[retry_hmc$param != "lp__", ]

g <- ggplot() +
    stat_halfeye(
        data = retry_advi,
        mapping = aes(
            xdist = dist_normal(mean, sd),
            y = "advi", colour = "advi", fill = "advi"
        )
    ) +
    stat_halfeye(
        data = retry_hmc,
        mapping = aes(
            xdist = dist_normal(mean, sd),
            y = "hmc", colour = "hmc", fill = "hmc"
        )
    ) +
    facet_wrap(~param, scales = "free_x")

ggsave("tmp.pdf")
system("convert tmp.pdf tmp.png")



fit <- fit_mod(
    "noGT",
    top_newd[1, "gene"],
    top_newd[1, "snp"],
    top_newd[["condition"]][[1]],
    chains = 1,
    iter = 10
)


bj <- retry_hmc["bj", "mean"]
b1 <- retry_hmc["betas[1]", "mean"]
b2 <- retry_hmc["betas[2]", "mean"]
phi <- retry_hmc["phi", "mean"]
theta <- retry_hmc["theta", "mean"]

rbj <- range(c(retry_hmc["bj", "mean"], retry_advi["bj", "mean"]))
rb1 <- range(c(retry_hmc["betas[1]", "mean"], retry_advi["betas[1]", "mean"]))
rb2 <- range(c(retry_hmc["betas[2]", "mean"], retry_advi["betas[2]", "mean"]))
rp <- range(c(retry_hmc["phi", "mean"], retry_advi["phi", "mean"]))
rt <- range(c(retry_hmc["theta", "mean"], retry_advi["theta", "mean"]))

n_per_param <- 200
grid <- expand.grid(
    bj = seq(rbj[[1]], rbj[[2]], length.out = n_per_param),
    beta1 = seq(rb1[[1]], rb1[[2]], length.out = n_per_param),
    beta2 = mean(rb2),
    phi = mean(rp),
    theta = mean(rt)
)
rai <- NULL
grid$lp <- sapply(1:nrow(grid), function(i) {
    pars <- list(
        bj = grid[i, "bj"],
        betas = c(grid[i, "beta1"], grid[i, "beta2"]),
        phi = grid[i, "phi"],
        rai0 = rai,
        theta = summary(fit)[[1]]["theta", "mean"]
    )
    log_prob(
        fit,
        upars = unconstrain_pars(fit, pars)
    )
})
g <- ggplot(grid) +
    aes(x = bj, y = beta1, fill = log(abs(lp))) +
    geom_tile() +
    geom_point(
        aes(
            x = retry_hmc["bj", "mean"],
            y = retry_hmc["betas[1]", "mean"],
            colour = "HMC"
        )
    ) +
    geom_point(
        aes(
            x = retry_advi["bj", "mean"],
            y = retry_advi["betas[1]", "mean"],
            colour = "VB"
        )
    ) +
    scale_fill_viridis()

ggsave("tmp2.pdf")
system("convert tmp2.pdf tmp2.png")



grid <- expand.grid(
    bj = seq(rbj[[1]], rbj[[2]], length.out = n_per_param),
    beta1 = mean(rb1),
    beta2 = seq(rb2[[1]], rb2[[2]], length.out = n_per_param),
    phi = mean(rp),
    theta = mean(rt)
)

grid$lp <- sapply(1:nrow(grid), function(i) {
    pars <- list(
        bj = grid[i, "bj"],
        betas = c(grid[i, "beta1"], grid[i, "beta2"]),
        phi = grid[i, "phi"],
        rai0 = rai,
        theta = summary(fit)[[1]]["theta", "mean"]
    )
    log_prob(
        fit,
        upars = unconstrain_pars(fit, pars)
    )
})

g <- ggplot(grid) +
    aes(x = bj, y = beta2, fill = log(abs(lp))) +
    geom_tile() +
    geom_point(
        aes(
            x = retry_hmc["bj", "mean"],
            y = retry_hmc["betas[2]", "mean"],
            colour = "HMC"
        )
    ) +
    geom_point(
        aes(
            x = retry_advi["bj", "mean"],
            y = retry_advi["betas[2]", "mean"],
            colour = "VB"
        )
    ) +
    scale_fill_viridis()
# ggsave(sprintf("%s/%s/diag/%s_%s_%s_grid-beta2.png", fpath, model, mygene, mysnp, mycondition))
ggsave("tmp3.pdf")
system("convert tmp3.pdf tmp3.png")


if (model == "noGT" & !file.exists("rds/noGT_worst_newmod.rds")) {
    newmod <- stan_model("src/stan/noGT_nb_ase.stan")
    gene_data <- get_gene_data(top_newd$gene[[1]], model)

    retry_newmod <- replicate(n_replicates,
        {
            d <- fit_fun(
                gene_data,
                top_newd$snp[[1]],
                gene = top_newd$gene[[1]],
                init = "random",
                model = newmod,
                vars = NULL
            )
            d <- d[d$condition == top_newd$condition[[1]], ]
            d <- d[, !grepl("%", colnames(d))]
            d$param <- gsub("1$", "", rownames(d))
            d
        },
        simplify = FALSE
    )
    newmod_df <- do.call(rbind, retry_newmod)
    saveRDS(newmod_df, "rds/noGT_worst_newmod.rds")
}
newmod_df <- readRDS("rds/noGT_worst_newmod.rds")




# dfma <- merge(df_new, df, suffix = c("new", "old"))
# dfma <- dfma[order(dfma$hmc), ]
# dfma$ind <- 1:nrow(dfma)

# g <- ggplot(dfma) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = vb_low,
#             y = vb,
#             ymax = vb_high,
#             colour = "vb"
#         ),
#         position = position_nudge(x = 0.25)
#     ) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = `5.0%`,
#             y = mean,
#             ymax = `95.0%`,
#             colour = "vb_rerun"
#         ),
#         position = position_nudge(x = -0.25)
#     ) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = hmc_low,
#             y = hmc,
#             ymax = hmc_high,
#             colour = "hmc"
#         )
#     )
# ggsave(sprintf("%s/%s/diag/discrepant_posteriors_out.png", fpath, model), width = 20, height = 7)


# dfmma <- dfma[abs(dfma$discrepancy) != max(abs(dfma$discrepancy)), ]

# g <- ggplot(dfmma) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = vb_low,
#             y = vb,
#             ymax = vb_high,
#             colour = "vb"
#         ),
#         position = position_nudge(x = 0.25)
#     ) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = `5.0%`,
#             y = mean,
#             ymax = `95.0%`,
#             colour = "vb_rerun"
#         ),
#         position = position_nudge(x = -0.25)
#     ) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = hmc_low,
#             y = hmc,
#             ymax = hmc_high,
#             colour = "hmc"
#         )
#     )
# ggsave(sprintf("%s/%s/diag/discrepant_posteriors.png", fpath, model), width = 20, height = 7)



# df_lots <- lapply(1:5, function(j) {
#     new_res <- parallel::mclapply(1:nrow(df), function(i) {
#         fit_fun(df[i, "gene"], df[i, "snp"])
#     }, mc.cores = 5)
#     df_new <- do.call(rbind, new_res)
#     dfm <- merge(df_new[, c("mean", "khat", "time", "gene", "snp")], df)
#     dfm <- dfm %>% mutate(new_discrepancy = hmc - mean)
#     dfm
# })
