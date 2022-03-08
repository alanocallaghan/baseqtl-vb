library("baseqtl")
library("argparse")
library("rstan")

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--method",
    type = "character"
)
args <- parser$parse_args()

optimizing <- function(...) {
    rstan::optimizing(..., draws = 1000)
}
vb <- function(...) {
    # rstan::vb(..., grad_samples = 5, elbo_samples = 1000, tol_rel_obj = 1e-3)
    while (TRUE) {
        f <- try(rstan::vb(..., tol_rel_obj = 1e-3, iter = 20000))
        if (!inherits(f, "try-error")) {
            return (f)
        }
    }
}
sampling <- function(...) {
    rstan::sampling(..., chains = 1, open_progress = FALSE)
}

if (is.null(method <- args[["method"]])) {
    method <- "vb"
}

fun <- match.fun(method)

# probs <- c(0.005, 0.025, 0.25, 0.50, 0.75, 0.975, 0.995)
probs <- seq(0.005, 0.995, by = 0.005)
dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

files <- list.files(dir)
files <- grep("rbias", files, value=TRUE)
stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))

ls_mat <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")

disc <- readRDS("rds/GT_discrepancies_vb.rds")

fit_stan <- function(i) {
    line <- disc[i, ]
    gene <- line[["gene"]]
    snp <- line[["snp"]]
    file <- sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
    gene_data <- readRDS(file)
    snp_in <- gene_data[[snp]]
    gene_datas <- readRDS(file)
    snps_data <- lapply(
        gene_datas,
        function(x) {
            out <- in.neg.beta.prob.eff2(x)
            out$cov <- cbind(out$cov, ls_mat[, gene])
            out$K <- 2
            out$k <- 2
            out$aveP <- c(0, 0)
            out$sdP <- c(0.0309, 0.3479)
            out$mixP <- c(0.97359164, 0.02640836)
            capture.output(
                time <- system.time(
                    post <- fun(
                        baseqtl:::stanmodels$GT_nb_ase_refbias, data = out
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
                tab <- tab[tab$variable == "bj", ]
            } else {
                tab <- rstan::summary(
                    post,
                    pars = "bj",
                    probs = probs
                )$summary
            }
            tab <- as.data.frame(tab)
            tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
            tab$gene <- gene
            tab$time <- time[["elapsed"]]
            tab
        }
    )
    names(snps_data) <- names(gene_datas)
    snps_data
}
all_stan_res <- lapply(1:nrow(disc), fit_stan)



names(all_stan_res) <- genes



saveRDS(all_stan_res, sprintf("rds/%s/disc.rds", method))

