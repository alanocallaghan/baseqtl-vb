library("baseqtl")
library("argparse")
library("rstan")

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "noGT",
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
    default = "ENSG00000010810",
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
tol <- args[["tolerance"]]
n_iterations <- args[["n_iterations"]]
gene <- args[["gene"]]
optimizing <- function(...) {
    rstan::optimizing(..., draws = 1000)
}
vb <- function(...) {
    # rstan::vb(..., grad_samples = 5, elbo_samples = 1000, tol_rel_obj = 1e-3)
    for (i in 1:10) {
        f <- try(rstan::vb(..., tol_rel_obj = tol, iter = n_iterations))
        if (!inherits(f, "try-error")) {
            return (f)
        }
    }
}
sampling <- function(...) {
    rstan::sampling(..., chains = 4, open_progress = FALSE)
}

method <- args[["inference"]]
model <- args[["model"]]
fun <- match.fun(method)

probs <- seq(0.005, 0.995, by = 0.005)

mtol <- if (method == "vb") sprintf("%s_%1.0e", method, tol) else method

# mtol <- method

if (model == "GT") {
    probs <- seq(0.005, 0.995, by = 0.005)
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

    
    ls_mat <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")
    
    gene_datas <- readRDS(
        sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
    )
    snps <- names(gene_datas)
    
    fit_stan <- function(gene, snp) {
        out <- in.neg.beta.prob.eff2(gene_datas[[snp]])
        out$cov <- cbind(out$cov, ls_mat[, gene])
        out$K <- 2
        out$k <- 2
        out$aveP <- c(0, 0)
        out$sdP <- c(0.0309, 0.3479)
        out$mixP <- c(0.97359164, 0.02640836)
        t0 <- proc.time()
        elbo_text <- capture.output(
            post <- fun(
                baseqtl:::stanmodels$GT_nb_ase_refbias, data = out
            )
        )
        time <- proc.time() - t0
        if (method == "optimizing") {
            tab <- posterior::summarise_draws(
                post$theta_tilde,
                mean,
                sd,
                ~quantile(.x, probs = probs)
            )
            tab <- tab[tab$variable == "bj", ]
            prob <- mean(post$theta_tilde[, "bj"] > 0)
        } else {
            tab <- rstan::summary(
                post,
                pars = "bj",
                probs = probs
            )$summary
            prob <- mean(extract(post, pars="bj")$bj > 0)
        }
        tab <- as.data.frame(tab)
        if (prob < 0.5) {
            prob <- 1 - prob
        }
        tab$prob <- prob
        if (method == "vb") {
            elbo <- parse_elbo(elbo_text)
            attr(tab, "elbo") <- elbo
            tab$converged <- elbo[nrow(elbo), "iter"] != n_iterations
            tab$niter <- elbo[nrow(elbo), "iter"]
        }
        tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
        tab$gene <- gene
        tab$time <- time[["elapsed"]]
        tab$snp <- snp
        tab
    }
} else {
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"

    covariates <- list(
        normal_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/normal_skin_scaled_gc_libsize.rds"),
        Psoriasis_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/Psoriasis_skin_scaled_gc_libsize.rds")
    )
    files <- c(
        normal_skin = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, gene),
        Psoriasis_skin = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, gene)
    )
    gene_datas <- lapply(
        files,
        readRDS
    )
    snps_each <- lapply(gene_datas, names)
    snps <- do.call(intersect, unname(snps_each))

    fit_stan <- function(gene, snp) {
        tabs <- lapply(names(files),
            function(condition) {
                data <- in.neg.beta.noGT.eff2(
                    gene_datas[[condition]][[snp]],
                    covar = covariates[[condition]][names(gene_datas[[condition]][[snp]]$NB$counts), gene, drop = FALSE]
                )
                data$k <- 3
                data$aveP <- c(0, 0, 0)
                data$sdP <- c(0.0436991990773286, 0.34926955206545, 0.4920048983496)
                data$mixP <- c(-0.0460439385014068, -3.50655789731998, -4.19970507787993)

                model <- baseqtl:::stanmodels$noGT_nb_ase
                if (!is.null(data$ai0)) {
                    model <- baseqtl:::stanmodels$noGT_nb_ase_refbias
                }
                t0 <- proc.time()
                elbo_text <- capture.output(
                    post <- fun(
                        model, data = data
                    )
                )
                time <- proc.time() - t0
                if (method == "optimizing") {
                    tab <- posterior::summarise_draws(
                        post$theta_tilde,
                        mean,
                        sd,
                        ~quantile(.x, probs = probs)
                    )
                    tab <- tab[tab$variable %in% "bj", ]
                    prob <- mean(post$theta_tilde[, "bj"] > 0)
                } else {
                    tab <- rstan::summary(
                        post,
                        pars = "bj",
                        probs = probs
                    )$summary
                    prob <- mean(extract(post, pars="bj")$bj > 0)
                }
                tab <- as.data.frame(tab)
                tab$prob <- ifelse(prob < 0.5, 1 - prob, prob)
                if (method == "vb") {
                    elbo <- parse_elbo(elbo_text)
                    attr(tab, "elbo") <- elbo
                    tab$converged <- elbo[nrow(elbo), "iter"] != n_iterations
                    tab$niter <- elbo[nrow(elbo), "iter"]
                }
                tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
                tab$gene <- gene
                tab$snp <- snp
                tab$time <- time[["elapsed"]]
                tab$condition <- condition
                tab
            }
        )
        do.call(rbind, tabs)
    }
}

tmp <- lapply(snps,
    function(snp) {
        if (model == "GT") {
            file <- sprintf("rds/GT/%s/%s_%s.rds", mtol, gene, snp)
        } else {
            file <- sprintf("rds/noGT/%s/%s_%s.rds", mtol, gene, snp)
        }
        if (file.exists(file)) return()
        tab <- fit_stan(gene, snp)
        saveRDS(tab, file)
    }
)

print(sprintf("rds/%s/%s/%s_done", model, mtol, gene))
file.create(sprintf("rds/%s/%s/%s_done", model, mtol, gene))
cat("Done!\n")
