library("baseqtl")
library("argparse")
library("rstan")
library("parallel")


parser <- ArgumentParser()
parser$add_argument( 
    "-i", "--inference",
    default = "vb",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-2,
    type = "double"
)
args <- parser$parse_args()

tol <- args[["tolerance"]]
method <- args[["inference"]]
n_iterations <- 50000

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

fun <- match.fun(method)

source("generatives.R")

# sim_mod <- cmdstan_model("gt_nb_sim.stan")
# gt_nb_mod <- cmdstan_model("gt_nb.stan")
# gt_nb_ase_mod <- cmdstan_model("gt_nb_ase.stan")

N <- 200
n_sims <- 500

sims_gt_nb_ase <- mclapply(
    seq_len(n_sims),
    function(i) {
        while (TRUE) {
            sim <- gt_nb_ase_gen(N)
            if (
                max(sim$Y, na.rm = TRUE) < .Machine$integer.max &
                all(sim$Y > 0)
            ) break
        }
        data <- list(
            N = N,
            A = N,
            L = N,
            K = 2,
            k = 2,
            Y = sim$Y,
            g = sim$g,
            gase = sim$g,
            m = sim$m,
            n = sim$n,
            pH = rep(1, N),
            ai0 = rep(0, N),
            sdai0 = rep(0.01017496, N),
            s = rep(1, N),
            cov = sim$cov,
            aveP = c(0, 0),
            sdP = c(0.0309, 0.3479),
            mixP = c(0.9736, 0.0264)
        )
        fit <- sampling(
            baseqtl:::stanmodels$GT_nb_ase,
            data = data,
            cores = 1
        )
        summ <- summary(fit, variables = "bj")
        qc_pass <- summ$summary["bj", "Rhat"] < 1.05 &
            summ$summary["bj", "n_eff"] > 500
        data.frame(
            true = sim$bj,
            est = summ$summary["bj", "mean"],
            null = sim$null,
            qc_pass = qc_pass,
            rhat = summ$summary["bj", "Rhat"],
            ess = summ$summary["bj", "n_eff"]
        )
    },
    mc.cores = 6
)

df_ase <- do.call(rbind, sims_gt_nb_ase)

plot(df_ase[1:2][df_ase[[3]], ])
plot(df_ase[1:2])

