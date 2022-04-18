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

components <- c("inter", "intra", "both")

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

# probs <- c(0.005, 0.025, 0.25, 0.50, 0.75, 0.975, 0.995)

mtol <- if (method == "vb") sprintf("%s_%1.0e", method, tol) else method

# mtol <- method


if (model == "GT") {
    probs <- seq(0.005, 0.995, by = 0.005)
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

    files <- list.files(dir)
    files <- grep("rbias", files, value=TRUE)
    stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)

    ls_mat <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")
    gene_datas <- readRDS(
        sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
    )
    snps <- names(gene_datas)
    fit_stan <- function(gene) {
        out <- in.neg.beta.prob.eff2(gene_datas[[snp]])
        out$cov <- cbind(out$cov, ls_mat[, gene])
        out$K <- 2
        out$k <- 2
        out$aveP <- c(0, 0)
        out$sdP <- c(0.0309, 0.3479)
        out$mixP <- c(0.97359164, 0.02640836)
        while (TRUE) {
            t0 <- proc.time()
            elbo_text <- capture.output(
                post <- try(
                    fun(
                        mod, data = out
                    )
                )
            )
            time <- proc.time() - t0
            if (!inherits(post, "try-error")) break
        }
        if (method == "optimizing") {
            tab <- posterior::summarise_draws(
                post$theta_tilde,
                mean,
                sd,
                ~quantile(.x, probs = probs)
            )
            tab <- tab[tab$variable == "bj", ]
        } else {
            tab <- rstan::summary(
                post,
                pars = "bj",
                probs = probs
            )$summary
        }
        tab <- as.data.frame(tab)
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
        tab$component <- component
        tab
    }
} else {
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"
    files <- list.files(dir, full.names = TRUE)

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
    snps <- do.call(intersect, snps_each)

    fit_stan <- function(gene) {
        lapply(names(files),
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
                elbo_text <- capture.output(
                    time <- system.time(
                        post <- fun(
                            model, data = data
                        )
                    )
                )
                if (method == "optimizing") {
                    tab <- posterior::summarise_draws(
                        post$theta_tilde,
                        mean,
                        sd,
                        ~quantile(.x, probs = probs)
                    )
                    tab <- tab[tab$variable %in% "bj", ]
                } else {
                    tab <- rstan::summary(
                        post,
                        pars = "bj",
                        probs = probs
                    )$summary
                }
                tab <- as.data.frame(tab)
                if (method == "vb") {
                    elbo <- parse_elbo(elbo_text)
                    attr(tab, "elbo") <- elbo
                    if (elbo[nrow(elbo), "iter"] == n_iterations) {
                        tab$converged <- FALSE
                    }
                }
                tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
                tab$gene <- gene
                tab$snp <- snp
                tab$time <- time[["elapsed"]]
                tab$condition <- condition
                file <- sprintf("rds/noGT/%s/%s_%s_%s.rds", mtol, gene, snp, condition)
                saveRDS(tab, file)
                tab
            }
        )
    }
}

tmp <- lapply(snps,
    function(snp) {
        file <- sprintf("rds/GT/components/%s/%s_%s.rds", mtol, args[["gene"]], snp)
        if (file.exists(file)) return()
        tabs <- list()
        for (component in components) {  
            modfile <- sprintf("src/stan/GT_nb_ase_refbias_%s.stan", component)
            mod <- stan_model(modfile)

            tabs[[component]] <- fit_stan(args[["gene"]], snp)
        }
        tab <- do.call(rbind, tabs)


        dir.create(
            sprintf("rds/GT/components/%s/", mtol),
            showWarnings = FALSE,
            recursive = TRUE
        )
        saveRDS(tab, file)
    }
)
file.create(sprintf("rds/GT/components/%s/%s_done", mtol, args[["gene"]]))
cat("Done!\n")
