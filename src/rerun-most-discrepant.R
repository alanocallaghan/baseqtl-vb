
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
    default = "GT",
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

init_types <- c("random", "fixed", "old")
rerun_files <- sprintf(
    "rds/%s_discrepancy_%s_rerun_%s.rds",
    model, tol_str, init_types
)


df <- readRDS(infile)


n_replicates <- 20


init_types <- c("random", "fixed", "old")
rerun_files <- sprintf(
    "rds/%s_discrepancy_%s_rerun_%s.rds",
    model, tol_str, init_types
)
fit_fun <- match.fun(paste("fit_stan", model, sep = "_"))
set.seed(42)

if (!file.exists(f <- sprintf("rds/%s_full_posterior_discrepant_hmc.rds", model))) {
    full_posteriors_hmc <- parallel::mclapply(
        seq_len(nrow(df)),
        function(i) {
            gene_data <- get_gene_data(df[i, "gene"], model)
            d <- fit_fun(
                gene_data,
                df[i, "snp"],
                seed = i,
                gene = df[i, "gene"],
                method = "sampling",
                cores = 1,
                summarise_posterior = FALSE,
                vars = NULL
            )
        }
    )
    saveRDS(
        full_posteriors_hmc,
        f
    )
}
if (!file.exists(f <- sprintf("rds/%s_full_posterior_discrepant_vb.rds", model))) {
    full_posteriors_vb <- parallel::mclapply(
        seq_len(nrow(df)),
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
        }
    )
    saveRDS(
        full_posteriors_vb,
        f
    )
}

if (!file.exists(f <- rerun_files[[1]])) {
    print("rerunning with random init")
    rand_res <- parallel::mclapply(
        seq_len(nrow(df)),
        function(i) {
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
                    d$param <- gsub("(normal_skin|Psoriasis_skin)?\\.?", "", rownames(d))
                    d
                },
                simplify = FALSE
            )
            do.call(rbind, out)
        }
    )
    rand_df <- do.call(rbind, rand_res)
    rand_df$init <- "random"
    saveRDS(rand_df, f)
}

if (!file.exists(f <- rerun_files[[2]])) {
    print("running with fixed init")
    fixed_res <- parallel::mclapply(
        seq_len(nrow(df)),
        function(i) {
            gene_data <- get_gene_data(df[i, "gene"], model)
            out <- replicate(
                n_replicates,
                {
                    d <- fit_fun(
                        gene_data,
                        df[i, "snp"],
                        gene = df[i, "gene"],
                        init = "fixed",
                        vars = NULL
                    )
                    d$param <- gsub("(normal_skin|Psoriasis_skin)?\\.?", "", rownames(d))
                    d
                },
                simplify = FALSE
            )
            do.call(rbind, out)
        }
    )
    fixed_df <- do.call(rbind, fixed_res)
    fixed_df$init <- "fixed"
    saveRDS(fixed_df, f)
}

if (!file.exists(f <- rerun_files[[3]])) {
    print("running with old init")
    old_res <- parallel::mclapply(
        seq_len(nrow(df)),
        function(i) {
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
                    d <- fit_fun(
                        gene_data,
                        df[i, "snp"],
                        gene = df[i, "gene"],
                        init = init,
                        vars = NULL
                    )
                    d$param <- gsub("(normal_skin|Psoriasis_skin)?\\.?", "", rownames(d))
                    d
                },
                simplify = FALSE
            )
            do.call(rbind, out)
        }
    )
    old_df <- do.call(rbind, old_res)
    old_df$init <- "old"
    saveRDS(old_df, f)
}

if (!file.exists(f <- sprintf("rds/%s_rerun_discrepant_full.rds", model))) {
    print("rerunning HMC")
    rerun_hmc <- parallel::mclapply(
        seq_len(nrow(df)),
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
            d$param <- gsub("(normal_skin|Psoriasis_skin)?\\.?", "", rownames(d))
            d
        }
    )
    rerun_df <- do.call(rbind, rerun_hmc)
    saveRDS(rerun_df, f)
}



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
