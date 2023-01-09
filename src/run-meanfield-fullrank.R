### comparing fullrank to meanfield approximations
### pick 1 association per gene
### compare posterior means and variances

library("baseqtl")
library("argparse")
library("rstan")
library("BiocParallel")

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    # default = "GT",
    default = "noGT",
    type = "character"
)
parser$add_argument(
    "-n", "--n_iterations",
    default = 50000,
    type = "integer"
)
parser$add_argument(
    "-c", "--cores",
    default = 1,
    type = "integer"
)
parser$add_argument(
    "-s", "--seed",
    default = 42,
    type = "integer"
)
parser$add_argument(
    "-g", "--gene",
    default = "ENSG00000130363", ## nogt
    # default = "ENSG00000025708", ## gt
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-3,
    type = "double"
)

args <- parser$parse_args()


source("src/functions.R")
tol <- args[["tolerance"]]
n_iterations <- args[["n_iterations"]]
gene <- args[["gene"]]
seed <- args[["seed"]]
method <- "vb"
model <- args[["model"]]
fun <- match.fun(method)

register(MulticoreParam(workers = args[["cores"]]))
# register(SerialParam())

mtol <- if (method == "vb") sprintf("%s_%1.0e", method, tol) else method

fit_fun <- if (model == "GT") fit_stan_GT else fit_stan_noGT

covariates <- get_covariates(model)
gene_data <- get_gene_data(gene, model)
snps <- get_snps(gene_data)


use <- seq_len(min(10, length(snps)))
res <- bplapply(snps[use],
    function(snp) {
        types <- c("meanfield", "fullrank")
        set.seed(seed)
        out <- lapply(
            types,
            function(type) {
                cat("Running", snp, type, "\n")
                tab <- fit_fun(
                    gene_data = gene_data,
                    gene = gene,
                    snp = snp,
                    covariates = covariates,
                    method = method,
                    algorithm = type,
                    summarise_posterior = FALSE,
                    tol = tol
                )
                tab
            }
        )
        names(out) <- types
        out[["hmc"]] <- fit_fun(
            gene_data = gene_data,
            gene = gene,
            snp = snp,
            covariates = covariates,
            method = "sampling",
            summarise_posterior = FALSE,
            cores = 1,
            tol = tol
        )
        out
    }
)
names(res) <- snps[use]
f <- sprintf("rds/%s/%s/meanfield/%s_s%d.rds", model, mtol, gene, seed)
saveRDS(res, f)
print(f)

cat("Done!\n")
