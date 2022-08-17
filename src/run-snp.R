library("baseqtl")
library("argparse")
library("rstan")

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "GT",
    # default = "noGT",
    type = "character"
)
parser$add_argument(
    "-i", "--inference",
    # default = "vb",
    default = "pathfinder",
    type = "character"
)
parser$add_argument(
    "-n", "--n_iterations",
    default = 50000,
    type = "integer"
)
parser$add_argument(
    "-c", "--cores",
    default = 50000,
    type = "integer"
)
parser$add_argument(
    "-g", "--gene",
    # default = "ENSG00000002330", ## nogt
    default = "ENSG00000025708", ## gt
    type = "character"
)
# parser$add_argument(
#     "-s", "--snp",
#     default = "17718119:G:T",
#     type = "character"
# )
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-3,
    type = "double"
)

args <- parser$parse_args()


source("src/functions.R")
source("src/pathfinder.R")
tol <- args[["tolerance"]]
n_iterations <- args[["n_iterations"]]
gene <- args[["gene"]]
method <- args[["inference"]]
model <- args[["model"]]
fun <- match.fun(method)


mtol <- if (method == "vb") sprintf("%s_%1.0e", method, tol) else method

# mtol <- method

fit_fun <- if (model == "GT") fit_stan_GT else fit_stan_noGT

covariates <- get_covariates(model)
gene_data <- get_gene_data(gene, model)
snps <- get_snps(gene_data)

tmp <- parallel::mclapply(snps,
    function(snp) {
        # cat("Running", snp, "\n")
        if (model == "GT") {
            file <- sprintf("rds/GT/%s/%s_%s.rds", mtol, gene, snp)
        } else {
            file <- sprintf("rds/noGT/%s/%s_%s.rds", mtol, gene, snp)
        }
        if (file.exists(file)) return()
        tab <- fit_fun(
            gene_data = gene_data,
            gene = gene,
            snp = snp,
            covariates = covariates,
            method = method,
            tol = tol
        )
        saveRDS(tab, file)
    }, mc.cores = args[["cores"]]
)

print(sprintf("rds/%s/%s/%s_done", model, mtol, gene))
file.create(sprintf("rds/%s/%s/%s_done", model, mtol, gene))
cat("Done!\n")
