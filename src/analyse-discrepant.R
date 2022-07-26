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
library("latex2exp")
library("cowplot")
rstan_options(auto_write = TRUE)

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
    full_posteriors_hmc <- parallel::mclapply(
        1:nrow(df),
        function(i) {
            gene_data <- get_gene_data(df[i, "gene"], model)
            d <- fit_fun(
                gene_data,
                df[i, "snp"],
                gene = df[i, "gene"],
                method = "sampling",
                cores = 1,
                summarise_posterior = FALSE,
                vars = NULL
            )
        },
        mc.cores = 8
    )
    saveRDS(full_posteriors_hmc, sprintf("rds/%s_full_posterior_discrepant_hmc.rds", model))


    full_posteriors_vb <- parallel::mclapply(
        1:nrow(df),
        function(i) {
            gene_data <- get_gene_data(df[i, "gene"], model)
            d <- fit_fun(
                gene_data,
                df[i, "snp"],
                gene = df[i, "gene"],
                method = "vb",
                output_samples = 4000,
                summarise_posterior = FALSE,
                vars = NULL
            )
        },
        mc.cores = 8
    )
    saveRDS(full_posteriors_vb, sprintf("rds/%s_full_posterior_discrepant_vb.rds", model))

    df_slabplot <- lapply(1:length(full_posteriors_hmc), function(i) {
        rbind(
            data.frame(
                method = "hmc",
                bj = extract(full_posteriors_hmc[[i]], par="bj")$bj,
                association = as.character(i),
                chain = rep(1:4, each = 1000)
            ),
            data.frame(
                method = "vb",
                bj = extract(full_posteriors_vb[[i]], par="bj")$bj,
                association = as.character(i),
                chain = "vb"
            )
        )
    })
    df_slabplot <- do.call(rbind, df_slabplot)
    g <- ggplot(df_slabplot) +
        aes(x = association, y = bj, group = method, fill = chain, colour = chain, order = chain) +
        geom_dots(position = position_dodge(), scale = 0.8) +
        labs(x = "Association", y = TeX("$\\beta_j$")) +
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

    ggsave("tmp.png", width = 20, height = 5)

    ## filter to one association
    ## plot two density plots - one for VB, one stacked + split by colour for HMC

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
ggsave(sprintf("%s/%s/diag/init_sd_disc.pdf", fpath, model), width = 2.5, height = 2.5)

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
ggsave(sprintf("%s/%s/diag/init_mean_abs_disc.pdf", fpath, model), width = 2.5, height = 2.5)

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





# fr_df <- readRDS(sprintf("rds/%s_discrepant_%s_fr.rds", model, tol_str))
# # fr_df <- do.call(rbind, fr_df) 

# cols <- c("mean", "khat", "time", "gene", "snp")
# if (model == "noGT") cols <- c(cols, "condition")

# dfn <- fr_df[, cols]
# dfn[, c("mean", "khat", "time")] <- lapply(
#     dfn[, c("mean", "khat", "time")],
#     as.numeric
# )
# dfm <- merge(
#     dfn, df,
#     by = by,
#     suffixes = c(".new", ".old")
# )

# dfm <- dfm %>% mutate(
#     new_discrepancy = hmc - mean,
#     discrepancy = discrepancy
# )





# if (model == "GT") {
#     stop()
# }


top_newd <- dfm[order(-abs(dfm$new_discrepancy)), ]
top_newd <- top_newd[1:10, ]

cols <- if (model == "GT") c("gene", "snp") else c("gene", "snp", "condition")

inds <- rerun_df[, cols] == top_newd[1, cols, drop=TRUE]
top_newd_hmc <- rerun_df[apply(inds, 1, all), ]
top_newd_hmc <- top_newd_hmc[, !grepl("%", colnames(top_newd_hmc)), ]

retry_advi <- df_all[df_all$rowidx == top_newd$rowidx[[1]], ]
retry_hmc <- top_newd_hmc

# retry_advi$param <- rownames(retry_advi)
retry_advi <- retry_advi[retry_advi$param != "lp__", ]
# retry_hmc$param <- rownames(retry_hmc)
retry_hmc <- retry_hmc[retry_hmc$param != "lp__", ]

# g <- ggplot() +
#     stat_halfeye(
#         data = retry_advi,
#         mapping = aes(
#             xdist = dist_normal(mean, sd),
#             y = "advi", colour = "advi", fill = "advi"
#         )
#     ) +
#     stat_halfeye(
#         data = retry_hmc,
#         mapping = aes(
#             xdist = dist_normal(mean, sd),
#             y = "hmc", colour = "hmc", fill = "hmc"
#         )
#     ) +
#     facet_wrap(~param, scales = "free_x")

# ggsave("tmp.pdf")
# system("convert tmp.pdf tmp.png")



mysnp <- top_newd[1, "snp"]
mygene <- top_newd[1, "gene"]
mycondition <- top_newd[1, "condition"]

fit <- if (model == "noGT") {
    fit_mod(
        "noGT",
        mygene,
        mysnp,
        mycondition,
        chains = 1,
        iter = 10
    )
} else {
    fit_mod(
        "GT",
        mygene,
        mysnp,
        chains = 1,
        iter = 10
    )
}


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

rbj <- expanded_range(rbj, maxabs=20)
rb1 <- expanded_range(rb1)
rb2 <- expanded_range(rb2)

n_per_param <- 200
grid <- expand.grid(
    bj = seq(rbj[[1]], rbj[[2]], length.out = n_per_param),
    beta1 = seq(rb1[[1]], rb1[[2]], length.out = n_per_param),
    beta2 = mean(rb2),
    phi = mean(rp),
    theta = mean(rt)
)
rai <- if (model == "noGT") NULL else retry_hmc[grep("rai0", rownames(retry_hmc)), "mean"]
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
        upars = unconstrain_pars(fit, pars),
        adjust_transform = TRUE
    )
})
lpr <- range(grid$lp, finite = TRUE)
g1 <- ggplot(grid) +
    # aes(x = bj, y = beta1, fill = log(abs(lp))) +
    aes(x = bj, y = beta1, fill = lp) +
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
    scale_fill_viridis() +
    scale_colour_brewer(palette = "Set1", name = "Point estimate") +
    labs(x = TeX("$\\beta_j$"), y = TeX("$\\gamma_0$"))

fname <- if (model == "GT") {
        sprintf(
            "%s/%s/diag/%s_%s_grid-beta1.pdf",
            fpath, model, mygene, mysnp
        )
    } else {
        sprintf(
            "%s/%s/diag/%s_%s_%s_grid-beta1.pdf",
            fpath, model, mygene, mysnp, mycondition
        )
}
ggsave(fname, width = 5, height = 3)

# ggsave("tmp2.pdf")
# system("convert tmp2.pdf tmp2.png")



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
        upars = unconstrain_pars(fit, pars),
        adjust_transform = TRUE
    )
})
lpr <- range(c(lpr, grid$lp), finite = TRUE)

g2 <- ggplot(grid) +
    # aes(x = bj, y = beta2, fill = log(abs(lp))) +
    aes(x = bj, y = beta2, fill = lp) +
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
    scale_fill_viridis() +
    scale_colour_brewer(palette = "Set1", name = "Point estimate") +
    labs(x = TeX("$\\beta_j$"), y = TeX("$\\gamma_1$"))

fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/%s_%s_grid-beta2.pdf",
        fpath, model, mygene, mysnp
    )
} else {
    sprintf(
        "%s/%s/diag/%s_%s_%s_grid-beta2.pdf",
        fpath, model, mygene, mysnp, mycondition
    )
}

ggsave(fname, width = 5, height = 3)



grid <- expand.grid(
    bj = mean(rbj),
    beta1 = seq(rb1[[1]], rb1[[2]], length.out = n_per_param),
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
        upars = unconstrain_pars(fit, pars),
        adjust_transform = TRUE
    )
})
lpr <- range(c(lpr, grid$lp), finite = TRUE)

g3 <- ggplot(grid) +
    # aes(x = beta1, y = beta2, fill = log(abs(lp))) +
    aes(x = beta1, y = beta2, fill = lp) +
    geom_tile() +
    geom_point(
        aes(
            x = retry_hmc["betas[1]", "mean"],
            y = retry_hmc["betas[2]", "mean"],
            colour = "HMC"
        )
    ) +
    geom_point(
        aes(
            x = retry_advi["betas[1]", "mean"],
            y = retry_advi["betas[2]", "mean"],
            colour = "VB"
        )
    ) +
    scale_fill_viridis() +
    scale_colour_brewer(palette = "Set1", name = "Point estimate") +
    labs(x = TeX("$\\gamma_0$"), y = TeX("$\\gamma_1$"))


fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/%s_%s_grid-beta1-beta2.pdf",
        fpath, model, mygene, mysnp
    )
} else {
    sprintf(
        "%s/%s/diag/%s_%s_%s_grid-beta1-beta2.pdf",
        fpath, model, mygene, mysnp, mycondition
    )
}
ggsave(fname, width = 5, height = 3)

legend <- get_legend(
    g1 + 
        scale_fill_viridis(limits = lpr, name = "Log posterior density") +
        guides(color = guide_legend(nrow = 1)) +
        theme(
            legend.position = "bottom",
            legend.text = element_text(hjust = 1, angle = 45)
        )
)

t <- list(
    theme(legend.position = "none"),
    scale_fill_viridis(limits = lpr, name = "Log posterior density")
)
g_all <- plot_grid(
    plot_grid(g1 + t, g2 + t, g3 + t, labels = "AUTO", nrow = 1),
    legend,
    rel_heights = c(0.85, 0.15),
    nrow = 2
)
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/%s_%s_grid-all.pdf",
        fpath, model, mygene, mysnp
    )
} else {
    sprintf(
        "%s/%s/diag/%s_%s_%s_grid-all.pdf",
        fpath, model, mygene, mysnp, mycondition
    )
}
ggsave(fname, width = 8, height = 4)



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
# newmod_df <- readRDS("rds/noGT_worst_newmod.rds")

if (model == "noGT") set.seed(35)

fit_real_hmc <- if (model == "noGT") {
    fit_mod("noGT", mygene, mysnp, mycondition)
} else {
    fit_mod("GT", mygene, mysnp)
}
fit_real_vb <- if (model == "noGT") {
    fit_mod("noGT", mygene, mysnp, mycondition, method = "vb")
} else {
    fit_mod("GT", mygene, mysnp, method = "vb")
}
summary(fit_real_hmc, pars = "bj")[[1]]
summary(fit_real_vb, pars = "bj")[[1]]


dhmc <- extract(fit_real_hmc)
dvb <- extract(fit_real_vb)

g1 <- ggplot() +
    stat_halfeye(aes(x = dhmc$bj, y = "HMC")) +
    stat_halfeye(aes(x = dvb$bj, y = "VB")) +
    labs(x = TeX("$\\beta_j$"), y = "Inference\nmethod")
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/%s_%s_method.pdf",
        fpath, model, mygene, mysnp
    )
} else {
    sprintf(
        "%s/%s/diag/%s_%s_%s_method.pdf",
        fpath, model, mygene, mysnp, mycondition
    )
}
ggsave(fname, width = 2, height = 3)

inter <- stan_model("src/stan/GT_nb_ase_refbias_inter.stan")
intra <- stan_model("src/stan/GT_nb_ase_refbias_intra.stan")
between <- fit_mod("GT", mygene, mysnp, stanmodel = inter)
within <- fit_mod("GT", mygene, mysnp, stanmodel = intra)

dbetween <- extract(between)
dwithin <- extract(within)

g2 <- ggplot() +
    stat_halfeye(aes(x = dbetween$bj, y = "Between\n(normal mixture prior)")) +
    stat_halfeye(aes(x = dwithin$bj, y = "Within\n(normal mixture prior)")) +
    labs(x = TeX("$\\beta_j$"), y = "Model")
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/%s_%s_component.pdf",
        fpath, model, mygene, mysnp
    )
} else {
    sprintf(
        "%s/%s/diag/%s_%s_%s_component.pdf",
        fpath, model, mygene, mysnp, mycondition
    )
}
ggsave(fname, width = 2, height = 3)

inter_p <- stan_model("src/stan/GT_nb_ase_refbias_inter_prior.stan")
intra_p <- stan_model("src/stan/GT_nb_ase_refbias_intra_prior.stan")

betweenp <- fit_mod("GT", mygene, mysnp, stanmodel = inter_p)
withinp <- fit_mod("GT", mygene, mysnp, stanmodel = intra_p)

dbetweenp <- extract(betweenp)
dwithinp <- extract(withinp)

g3 <- ggplot() +
    stat_halfeye(aes(x = dbetweenp$bj, y = "Between\n(normal prior)")) +
    stat_halfeye(aes(x = dwithinp$bj, y = "Within\n(normal prior)")) +
    labs(x = TeX("$\\beta_j$"), y = "Model")
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/%s_%s_component_prior.pdf",
        fpath, model, mygene, mysnp
    )
} else {
    sprintf(
        "%s/%s/diag/%s_%s_%s_component_prior.pdf",
        fpath, model, mygene, mysnp, mycondition
    )
}
ggsave(fname, width = 2, height = 3)

lims <- xlim(-1.2, 1.2)
p <- plot_grid(g1 + lims, g2 + lims, g3 + lims, ncol = 1, labels = "AUTO", align = "v")

fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/%s_%s_method_comp_prior.pdf",
        fpath, model, mygene, mysnp
    )
} else {
    sprintf(
        "%s/%s/diag/%s_%s_%s_method_comp_prior.pdf",
        fpath, model, mygene, mysnp, mycondition
    )
}
ggsave(fname, width = 4, height = 5)

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
