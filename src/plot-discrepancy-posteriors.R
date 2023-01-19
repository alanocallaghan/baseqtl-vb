library("ggplot2")
library("ggdist")
library("ggrastr")
library("distributional")
library("ggpointdensity")
library("viridis")
library("dplyr")
library("argparse")
library("latex2exp")
library("cowplot")
library("baseqtl")
library("rstan")

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


top_newd <- dfm[order(-abs(dfm$new_discrepancy)), ]
top_newd <- top_newd[1:10, ]

cols <- if (model == "GT") c("gene", "snp") else c("gene", "snp", "condition")

inds <- rerun_df[, cols] == top_newd[1, cols, drop = TRUE]
top_newd_hmc <- rerun_df[apply(inds, 1, all), ]
top_newd_hmc <- top_newd_hmc[, !grepl("%", colnames(top_newd_hmc)), ]

retry_advi <- df_all[df_all$rowidx == top_newd$rowidx[[1]], ]
retry_hmc <- top_newd_hmc

# retry_advi$param <- rownames(retry_advi)
retry_advi <- retry_advi[retry_advi$param != "lp__", ]
# retry_hmc$param <- rownames(retry_hmc)
retry_hmc <- retry_hmc[retry_hmc$param != "lp__", ]


mysnp <- top_newd[1, "snp"]
mysnp_clean <- make.names(mysnp)
mygene <- top_newd[1, "gene"]
mycondition <- top_newd[1, "condition"]

if (model == "noGT") {
    retry_advi <- retry_advi[grep(mycondition, rownames(retry_advi)), ]
}


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


bj <- retry_hmc[retry_hmc$param == "bj", "mean"]
b1 <- retry_hmc[retry_hmc$param == "betas[1]", "mean"]
b2 <- retry_hmc[retry_hmc$param == "betas[2]", "mean"]
phi <- retry_hmc[retry_hmc$param == "phi", "mean"]
theta <- retry_hmc[retry_hmc$param == "theta", "mean"]

rbj <- range(c(retry_hmc[retry_hmc$param == "bj", "mean"], retry_advi[retry_advi$param == "bj", "mean"]))
rb1 <- range(c(retry_hmc[retry_hmc$param == "betas[1]", "mean"], retry_advi[retry_advi$param == "betas[1]", "mean"]))
rb2 <- range(c(retry_hmc[retry_hmc$param == "betas[2]", "mean"], retry_advi[retry_advi$param == "betas[2]", "mean"]))
rp <- range(c(retry_hmc[retry_hmc$param == "phi", "mean"], retry_advi[retry_advi$param == "phi", "mean"]))
rt <- range(c(retry_hmc[retry_hmc$param == "theta", "mean"], retry_advi[retry_advi$param == "theta", "mean"]))

# error is here
rbj <- expanded_range(rbj, maxabs = 20)
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
            x = retry_hmc[retry_hmc$param == "bj", "mean"],
            y = retry_hmc[retry_hmc$param == "betas[1]", "mean"],
            colour = "HMC"
        )
    ) +
    geom_point(
        aes(
            x = retry_advi[retry_advi$param == "bj", "mean"],
            y = retry_advi[retry_advi$param == "betas[1]", "mean"],
            colour = "VB"
        )
    ) +
    scale_fill_viridis() +
    scale_colour_brewer(palette = "Set1", name = "Point estimate") +
    labs(x = TeX("$\\beta_{aFC}$"), y = TeX("$\\gamma_0$"))
gt1 <- rasterise(g1, dpi = 300)

fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-beta1.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-beta1.pdf",
        fpath, model
    )
}
# ggsave(fname, width = 5, height = 3)





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
            x = retry_hmc[retry_hmc$param == "bj", "mean"],
            y = retry_hmc[retry_hmc$param == "betas[2]", "mean"],
            colour = "HMC"
        )
    ) +
    geom_point(
        aes(
            x = retry_advi[retry_advi$param == "bj", "mean"],
            y = retry_advi[retry_advi$param == "betas[2]", "mean"],
            colour = "VB"
        )
    ) +
    scale_fill_viridis() +
    scale_colour_brewer(palette = "Set1", name = "Point estimate") +
    labs(x = TeX("$\\beta_{aFC}$"), y = TeX("$\\gamma_1$"))
gt2 <- rasterise(g2, dpi = 300)

fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-beta2.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-beta2.pdf",
        fpath, model
    )
}

# ggsave(fname, width = 5, height = 3)



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
            x = retry_hmc[retry_hmc$param == "betas[1]", "mean"],
            y = retry_hmc[retry_hmc$param == "betas[2]", "mean"],
            colour = "HMC"
        )
    ) +
    geom_point(
        aes(
            x = retry_advi[retry_advi$param == "betas[1]", "mean"],
            y = retry_advi[retry_advi$param == "betas[2]", "mean"],
            colour = "VB"
        )
    ) +
    scale_fill_viridis() +
    scale_colour_brewer(palette = "Set1", name = "Point estimate") +
    labs(x = TeX("$\\gamma_0$"), y = TeX("$\\gamma_1$"))
gt3 <- rasterise(g3, dpi = 300)


fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-beta1-beta2.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-beta1-beta2.pdf",
        fpath, model
    )
}
# ggsave(file = fname, width = 5, height = 3)



legend <- get_legend(
    g1 +
        scale_fill_viridis(
            limits = lpr, name = "Log posterior density",
            guide = guide_colourbar(
                label.hjust = 1,
                label.theme = element_text(angle = 45, size = 10)
            )
        ) +
        guides(color = guide_legend(nrow = 1)) +
        theme(
            legend.position = "bottom"
            # , legend.text = element_text(hjust = 1, angle = 45)
        )
)

t <- list(
    theme(legend.position = "none"),
    scale_fill_viridis(
        limits = lpr, name = "Log posterior density"
    )
)
labels <- if (model == "GT") c(LETTERS[1:3]) else LETTERS[4:6]
g_all <- plot_grid(
    plot_grid(g1 + t, g2 + t, g3 + t, labels = labels, nrow = 1),
    legend,
    rel_heights = c(0.8, 0.2),
    nrow = 2
)
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-all.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-all.pdf",
        fpath, model
    )
}
print(fname)
ggsave(fname, width = 8, height = 5)

t <- list(
    theme(legend.position = "none"),
    scale_fill_viridis(
        limits = lpr, name = "Log posterior density"
    )
)
labels <- if (model == "GT") c(LETTERS[1:3]) else LETTERS[4:6]
g_all <- plot_grid(
    plot_grid(gt1 + t, gt2 + t, gt3 + t, labels = labels, nrow = 1),
    legend,
    rel_heights = c(0.8, 0.2),
    nrow = 2
)
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-all-rast.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_grid-all-rast.pdf",
        fpath, model
    )
}
print(fname)
ggsave(fname, width = 8, height = 5)


if (model == "noGT") set.seed(35)

fit_real_hmc <- fit_mod(model, mygene, mysnp, mycondition)
fit_real_vb <- fit_mod(model, mygene, mysnp, mycondition, method = "vb")

summary(fit_real_hmc, pars = "bj")[[1]]
summary(fit_real_vb, pars = "bj")[[1]]


dhmc <- extract(fit_real_hmc)
dvb <- extract(fit_real_vb)

g1 <- ggplot() +
    stat_halfeye(aes(x = dhmc$bj, y = "HMC")) +
    stat_halfeye(aes(x = dvb$bj, y = "VB")) +
    labs(x = TeX("$\\beta_{aFC}$"), y = "Inference\nmethod")
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_method.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_method.pdf",
        fpath, model
    )
}
# ggsave(fname, width = 2, height = 3)

rbias <- if (model == "noGT") "" else "_refbias"

inter <- stan_model(sprintf("src/stan/%s_nb_ase%s_inter.stan", model, rbias))
intra <- stan_model(sprintf("src/stan/%s_nb_ase%s_intra.stan", model, rbias))
between <- fit_mod(model, mygene, mysnp, condition = mycondition, stanmodel = inter)
within <- fit_mod(model, mygene, mysnp, condition = mycondition, stanmodel = intra)

dbetween <- extract(between)
dwithin <- extract(within)

g2 <- ggplot() +
    stat_halfeye(aes(x = dbetween$bj, y = "Between\n(mixture prior)")) +
    stat_halfeye(aes(x = dwithin$bj, y = "Within\n(mixture prior)")) +
    labs(x = TeX("$\\beta_{aFC}$"), y = "Model")
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_component.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_component.pdf",
        fpath, model
    )
}
# ggsave(fname, width = 2, height = 3)

set.seed(42)

inter_p <- stan_model(sprintf("src/stan/%s_nb_ase%s_inter_prior.stan", model, rbias))
intra_p <- stan_model(sprintf("src/stan/%s_nb_ase%s_intra_prior.stan", model, rbias))

betweenp <- fit_mod(model, mygene, mysnp, condition = mycondition, stanmodel = inter_p)
withinp <- fit_mod(model, mygene, mysnp, condition = mycondition, stanmodel = intra_p)

dbetweenp <- extract(betweenp)
dwithinp <- extract(withinp)

g3 <- ggplot() +
    stat_halfeye(aes(x = dbetweenp$bj, y = "Between\n(normal prior)")) +
    stat_halfeye(aes(x = dwithinp$bj, y = "Within\n(normal prior)")) +
    labs(x = TeX("$\\beta_{aFC}$"), y = "Model")
fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_component_prior.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_component_prior.pdf",
        fpath, model
    )
}
# ggsave(fname, width = 2, height = 3)

if (model == "GT") {
    lims <- xlim(-1.2, 1.2)
} else {
    lims <- NULL
}
p <- plot_grid(g1 + lims, g2 + lims, g3 + lims, ncol = 1, labels = "AUTO", align = "v")

fname <- if (model == "GT") {
    sprintf(
        "%s/%s/diag/most_discrepant_method_comp_prior.pdf",
        fpath, model
    )
} else {
    sprintf(
        "%s/%s/diag/most_discrepant_method_comp_prior.pdf",
        fpath, model
    )
}
print(fname)
cat(mygene, mysnp, mycondition, file = sprintf("%s/%s/diag/worst_gene.txt", fpath, model))
ggsave(fname, width = 4, height = 5)
