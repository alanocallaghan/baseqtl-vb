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

fit_fun <- match.fun(paste("fit_stan", model, sep = "_"))


old_res <- lapply(1:2, function(i) {
    gene_data <- get_gene_data(df[i, "gene"], model)
    init <- list(
        betas = c(5, 0),
        bj = df[i, "hmc"],
        phi = 10,
        theta = 10
    )
    out <- replicate(
        n_replicates,
        {
            fit_fun(
                gene_data, df[i, "snp"],
                gene = df[i, "gene"],
                init = init,
                method = "sampling"
            )
        },
        simplify = FALSE
    )
    do.call(rbind, out)
})
old_df <- do.call(rbind, old_res)
old_df$init <- "old"


rand_res <- lapply(1:2, function(i) {
    gene_data <- get_gene_data(df[i, "gene"], model)
    out <- replicate(
        n_replicates,
        fit_fun(
            gene_data,
            df[i, "snp"],
            gene = df[i, "gene"],
            init = "random",
            method = "sampling"
        ),
        simplify = FALSE
    )
    do.call(rbind, out)
})
rand_df <- do.call(rbind, rand_res)
rand_df$init <- "random"

dfb <- rbind(rand_df, old_df)

g <- ggplot(dfb) +
    aes(x = init, y = time) +
    geom_violin()

ggsave(paste0(model, ".png"),
    width = 5, height = 5
)
