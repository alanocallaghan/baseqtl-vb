library("baseqtl")
library("argparse")
library("rstan")

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "GT",
    type = "character"
)
parser$add_argument(
    "-i", "--inference",
    default = "vb",
    type = "character"
)
parser$add_argument(
    "-n", "--n_iterations",
    default = 50000,
    type = "double"
)
parser$add_argument(
    "-g", "--gene",
    default = "ENSG00000128311",
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
method <- args[["inference"]]
model <- args[["model"]]
gene <- args[["gene"]]
fun <- match.fun(method)

components <- c("inter", "intra", "both")


mtol <- if (method == "vb") sprintf("%s_%1.0e", method, tol) else method
fit_fun <- if (model == "GT") fit_stan_GT else fit_fun <- fit_stan_noGT

covariates <- get_covariates(model)
gene_data <- get_gene_data(gene, model)
snps <- get_snps(gene_data)


res <- lapply(snps,
    function(snp) {
        file <- sprintf("rds/%s/components/%s/%s_%s.rds", model, mtol, gene, snp)
        # if (file.exists(file)) return()
        tabs <- list()
        for (component in components) {  
            modfile <- sprintf("src/stan/%s_nb_ase_refbias_%s.stan", model, component)
            mod <- stan_model(modfile)

            tabs[[component]] <- fit_fun(
                gene_data = gene_data,
                gene = gene,
                snp = snp,
                model = mod,
                covariates = covariates,
                method = method
            )
        }
        tab <- do.call(rbind, tabs)
        dir.create(
            sprintf("rds/%s/components/%s/", model, mtol),
            showWarnings = FALSE,
            recursive = TRUE
        )
        # saveRDS(tab, file)
    }
)
file <- sprintf("rds/%s/components/%s/%s.rds", model, mtol, gene)
cat(file, "\n")
saveRDS(res, sprintf("rds/%s/components/%s/%s_done", model, mtol, gene))
cat("Done!\n")
