library("rstan")
library("ggplot2")
library("ggdist")
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
    "-i", "--index",
    default = 1,
    type = "integer"
)

args <- parser$parse_args()
source("src/functions.R")
theme_set(theme_bw())

tol <- args[["tolerance"]]
model <- args[["model"]]
method <- "vb"
index <- args[["index"]]


tol_str <- if (method == "vb") sprintf("%s_%1.0e", method, tol) else method

file <- sprintf("rds/%s_discrepancies_%s.rds", model, tol_str)
df <- readRDS(file)

n_per_param <- 100
fpath <- sprintf("fig_%1.0e", tol)
dir.create(fpath, showWarnings = FALSE, recursive = TRUE)
# fpath <- "fig"

top_df <- df[index, ]
mysnp <- top_df$snp
mygene <- top_df$gene
mycondition <- top_df$condition

fit <- fit_mod(model, mygene, mysnp, mycondition)


b1 <- summary(fit)[[1]]["betas[1]", "mean"]
b2 <- summary(fit)[[1]]["betas[2]", "mean"]
phi <- summary(fit)[[1]]["phi", "mean"]
theta <- summary(fit)[[1]]["theta", "mean"]


sp <- summary(fit)[[1]]
rai <- sp[grep("rai0", rownames(sp), value = TRUE), "mean"]


grid <- expand.grid(
    bj = seq(-abs(top_df$vb), abs(top_df$vb), length.out = n_per_param) * 1.5,
    beta1 = seq(2, b1 * 2, length.out = n_per_param),
    beta2 = b2,
    phi = phi,
    theta = theta
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
    aes(x = bj, y = beta1, fill = lp) +
    geom_tile() +
    geom_segment(
        x = top_df$hmc,
        y = min(grid$beta1),
        yend = max(grid$beta1),
        xend = top_df$hmc,
        linetype = "dashed",
        aes(colour = "HMC")
    ) +
    geom_segment(
        x = top_df$vb,
        y = min(grid$beta1),
        yend = max(grid$beta1),
        xend = top_df$vb,
        linetype = "dashed",
        aes(colour = "VB")
    ) +
    scale_fill_viridis() +
    labs(x = Tex("\\beta_j"), y = Tex("\\gamma_0"))
ggsave(sprintf("%s/%s/diag/%s_%s_%s_grid-beta1.png", fpath, model, mygene, mysnp, mycondition))



grid <- expand.grid(
    bj = seq(-abs(top_df$vb), abs(top_df$vb), length.out = n_per_param) * 2,
    beta1 = b1,
    beta2 = seq(-abs(b2), b2, length.out = n_per_param) * 2,
    phi = phi,
    theta = theta
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
    aes(x = bj, y = beta2, fill = lp) +
    geom_tile() +
    geom_segment(
        x = top_df$hmc,
        y = min(grid$beta2),
        yend = max(grid$beta2),
        xend = top_df$hmc,
        linetype = "dashed",
        aes(colour = "HMC")
    ) +
    geom_segment(
        x = top_df$vb,
        y = min(grid$beta2),
        yend = max(grid$beta2),
        xend = top_df$vb,
        linetype = "dashed",
        aes(colour = "VB")
    ) +
    scale_fill_viridis() +
    labs(x = Tex("\\beta_j"), y = Tex("\\gamma_1"))
ggsave(sprintf("%s/%s/diag/%s_%s_%s_grid-beta2.png", fpath, model, mygene, mysnp, mycondition))





# grid <- expand.grid(
#     bj = top_df$hmc,
#     beta1 = seq(-b1, b1, length.out = n_per_param) * 2,
#     beta2 = seq(-abs(b2), b2, length.out = n_per_param) * 2,
#     phi = phi,
#     theta = theta
# )
# grid$lp <- sapply(1:nrow(grid), function(i) {
#     pars <- list(
#         bj = grid[i, "bj"],
#         betas = c(grid[i, "beta1"], grid[i, "beta2"]),
#         phi = grid[i, "phi"],
#         theta = summary(fit)[[1]]["theta", "mean"]
#     )
#     log_prob(
#         fit,
#         upars = unconstrain_pars(fit, pars)
#     )
# })

# g <- ggplot(grid) +
#     aes(x = beta1, y = beta2, fill = lp) +
#     geom_tile() +
#     scale_fill_viridis()
# ggsave("grid-betas.png")





# grid <- expand.grid(
#     bj = seq(-abs(top_df$vb), abs(top_df$vb), length.out = n_per_param) * 2,
#     beta1 = b1,
#     beta2 = b2,
#     phi = seq(phi * 0.5, phi * 2, length.out = n_per_param),
#     theta = theta
# )

# grid$lp <- sapply(1:nrow(grid), function(i) {
#     pars <- list(
#         bj = grid[i, "bj"],
#         betas = c(grid[i, "beta1"], grid[i, "beta2"]),
#         phi = grid[i, "phi"],
#         theta = summary(post)[[1]]["theta", "mean"]
#     )
#     log_prob(
#         post,
#         upars = unconstrain_pars(post, pars)
#     )
# })

# g <- ggplot(grid) +
#     aes(x = bj, y = lp) +
#     geom_point()
# ggsave("grid-phi.png")




# grid <- expand.grid(
#     bj = seq(-abs(top_df$vb), abs(top_df$vb), length.out = n_per_param) * 2,
#     beta1 = b1,
#     beta2 = b2,
#     phi = phi,
#     theta = seq(theta * 0.5, theta * 2, length.out = n_per_param)
# )

# grid$lp <- sapply(1:nrow(grid), function(i) {
#     pars <- list(
#         bj = grid[i, "bj"],
#         betas = c(grid[i, "beta1"], grid[i, "beta2"]),
#         phi = grid[i, "phi"],
#         theta = summary(post)[[1]]["theta", "mean"]
#     )
#     log_prob(
#         post,
#         upars = unconstrain_pars(post, pars)
#     )
# })

# g <- ggplot(grid) +
#     aes(x = bj, y = lp) +
#     geom_point()
# ggsave("grid-theta.png")
