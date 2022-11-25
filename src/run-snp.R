library("baseqtl")
library("argparse")
library("rstan")
library("BiocParallel")

parser <- ArgumentParser()
parser$add_argument(
  "-m", "--model",
  default = "GT",
  # default = "noGT",
  type = "character"
)
parser$add_argument(
  "-i", "--inference",
  default = "sampling",
  # default = "pathfinder",
  type = "character"
)
parser$add_argument(
  "-n", "--n_iterations",
  default = 50000,
  type = "integer"
)
parser$add_argument(
  "-s", "--seed",
  # default = 42,
  type = "integer"
)
parser$add_argument(
  "-c", "--cores",
  default = 1,
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
  default = 1e-2,
  type = "double"
)

args <- parser$parse_args()

register(MulticoreParam(workers = args[["cores"]]))
source("src/functions.R")
source("src/pathfinder.R")
tol <- args[["tolerance"]]
seed <- args[["seed"]]
n_iterations <- args[["n_iterations"]]
gene <- args[["gene"]]
method <- args[["inference"]]
model <- args[["model"]]
fun <- match.fun(method)


mtol <- sprintf("%s_%1.0e", method, tol)

# mtol <- method

fit_fun <- if (model == "GT") fit_stan_GT else fit_stan_noGT

covariates <- get_covariates(model)
gene_data <- get_gene_data(gene, model)
snps <- get_snps(gene_data)


res <- bplapply(
  snps,
  # res <- lapply(snps,
  function(snp) {
    cat("Running", snp, "\n")
    tab <- fit_fun(
      gene_data = gene_data,
      gene = gene,
      snp = snp,
      covariates = covariates,
      method = method,
      seed = seed,
      tol = tol
    )
    tab <- as.data.frame(tab)
    tab$seed <- seed
    tab
    # }
  }
)

if (is.null(seed)) {
  f <- sprintf("rds/%s/%s/%s_sNULL.rds", model, mtol, gene)
} else {
  f <- sprintf("rds/%s/%s/%s_s%d.rds", model, mtol, gene, seed)
}
print(f)
saveRDS(res, f)
cat("Done!\n")
# ./noGT/vb_1e-03/ENSG00000002330_sNULL.rds
