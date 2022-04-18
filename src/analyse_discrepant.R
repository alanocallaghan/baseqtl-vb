library("rstan")
library("ggplot2")
library("ggdist")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
library("geomtextpath")
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
    default = 1e-3,
    type = "double"
)
args <- parser$parse_args()

tol <- args[["tolerance"]]
model <- args[["model"]]
tol_str <- sprintf("%1.0e", tol)

file <- sprintf("rds/%s_discrepancies_vb_%s.rds", model, tol_str)
df <- readRDS(file)


fpath <- sprintf("fig_%1.0e", tol)
dir.create(fpath)
# fpath <- "fig"


method <- "vb"
vb <- function(...) {
    rstan::vb(..., tol_rel_obj = tol, iter = 20000)
}
source("src/fit_funs.R")

new_res <- parallel::mclapply(1:nrow(df), function(i) {
    fit_stan(df[i, "gene"], df[i, "snp"])
}, mc.cores = 5)
df_new <- do.call(rbind, new_res)

dfm <- merge(df_new[, c("mean", "khat", "time", "gene", "snp")], df)

dfm <- dfm %>% mutate(
    new_discrepancy = hmc - mean,
    discrepancy = discrepancy
)
g <- ggplot(dfm) +
    aes(abs(discrepancy), abs(new_discrepancy)) +
    geom_point()
ggsave(sprintf("%s/%s/diag/rerun_discrepant.png", fpath, model), width = 5, height = 5)

dfmm <- dfm[dfm$discrepancy != max(dfm$discrepancy), ]
g <- ggplot(dfmm) +
    aes(abs(discrepancy), abs(new_discrepancy)) +
    geom_point()
ggsave(sprintf("%s/%s/diag/rerun_discrepant1.png", fpath, model), width = 5, height = 5)


dfma <- merge(df_new, df, suffix = c("new", "old"))
dfma <- dfma[order(dfma$hmc), ]
dfma$ind <- 1:nrow(dfma)

g <- ggplot(dfma) +
    geom_pointrange(
        aes(x = factor(ind),
            ymin = vb_low,
            y = vb,
            ymax = vb_high,
            colour = "vb"
        ),
        position = position_nudge(x = 0.25)
    ) +
    geom_pointrange(
        aes(x = factor(ind),
            ymin = `5.0%`,
            y = mean,
            ymax = `95.0%`,
            colour = "vb_rerun"
        ),
        position = position_nudge(x = -0.25)
    ) +
    geom_pointrange(
        aes(x = factor(ind),
            ymin = hmc_low,
            y = hmc,
            ymax = hmc_high,
            colour = "hmc"
        )
    )
ggsave(sprintf("%s/%s/diag/discrepant_posteriors_out.png", fpath, model), width = 20, height = 7)


dfmma <- dfma[abs(dfma$discrepancy) != max(abs(dfma$discrepancy)), ]

g <- ggplot(dfmma) +
    geom_pointrange(
        aes(x = factor(ind),
            ymin = vb_low,
            y = vb,
            ymax = vb_high,
            colour = "vb"
        ),
        position = position_nudge(x = 0.25)
    ) +
    geom_pointrange(
        aes(x = factor(ind),
            ymin = `5.0%`,
            y = mean,
            ymax = `95.0%`,
            colour = "vb_rerun"
        ),
        position = position_nudge(x = -0.25)
    ) +
    geom_pointrange(
        aes(x = factor(ind),
            ymin = hmc_low,
            y = hmc,
            ymax = hmc_high,
            colour = "hmc"
        )
    )
ggsave(sprintf("%s/%s/diag/discrepant_posteriors.png", fpath, model), width = 20, height = 7)



df_lots <- lapply(1:5, function(j) {
    new_res <- parallel::mclapply(1:nrow(df), function(i) {
        fit_stan(df[i, "gene"], df[i, "snp"])
    }, mc.cores = 5)
    df_new <- do.call(rbind, new_res)
    dfm <- merge(df_new[, c("mean", "khat", "time", "gene", "snp")], df)
    dfm <- dfm %>% mutate(new_discrepancy = hmc - mean)
    dfm
})
