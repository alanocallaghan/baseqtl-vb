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
    default = "noGT",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-2,
    type = "double"
)
args <- parser$parse_args()

theme_set(theme_bw())
tol <- args[["tolerance"]]
model <- args[["model"]]
tol_str <- sprintf("%1.0e", tol)

file <- sprintf("rds/%s_discrepancies_vb_%s.rds", model, tol_str)
df <- readRDS(file)

fpath <- sprintf("fig_%1.0e", tol)
dir.create(fpath, showWarnings = FALSE, recursive = TRUE)
# fpath <- "fig"
n_replicates <- 10

method <- "vb"
vb <- function(...) {
    rstan::vb(..., tol_rel_obj = tol, iter = 50000)
}
source("src/fit_funs.R")

rerun_file <- sprintf("rds/%s_discrepancy_%s_%s_rerun.rds", model, method, tol_str)

if (!file.exists(rerun_file)) {
    rand_res <- parallel::mclapply(1:nrow(df), function(i) {
        out <- replicate(
            n_replicates,
            fit_stan(df[i, "gene"], df[i, "snp"], init = "random"),
            simplify = FALSE
        )
        do.call(rbind, out)
    }, mc.cores = 8)
    rand_df <- do.call(rbind, rand_res)
    rand_df$init <- "random"

    opt_res <- parallel::mclapply(1:nrow(df), function(i) {
        out <- replicate(
            n_replicates,
            fit_stan(df[i, "gene"], df[i, "snp"], init = "opt"),
            simplify = FALSE
        )
        do.call(rbind, out)
    }, mc.cores = 8)
    opt_df <- do.call(rbind, opt_res)
    opt_df$init <- "opt"

    fixed_res <- parallel::mclapply(1:nrow(df), function(i) {
        out <- replicate(
            n_replicates,
            fit_stan(df[i, "gene"], df[i, "snp"], init = "fixed"),
            simplify = FALSE
        )
        do.call(rbind, out)
    }, mc.cores = 8)
    fixed_df <- do.call(rbind, fixed_res)
    fixed_df$init <- "fixed"


    df_new <- rbind(fixed_df, opt_df, rand_df)
    saveRDS(df_new, rerun_file)
}
df_new <- readRDS(rerun_file)

dfm <- merge(
    df_new[, c("mean", "khat", "time", "gene", "snp", "init")], df,
    by = c("gene", "snp"),
    suffixes = c(".new", ".old")
)

dfm <- dfm %>% mutate(
    new_discrepancy = hmc - mean,
    discrepancy = discrepancy
)

g <- ggplot(dfm) +
    aes(abs(discrepancy), abs(new_discrepancy), colour = init) +
    geom_abline(slope = 1, intercept = 0, linetype =  "dashed") +
    geom_jitter(size = 0.5, width = 0.001, height = 0) +
    labs(x = "Original discrepancy", y = "New discrepancy") +
    scale_colour_brewer(palette = "Set2", name = "Initialisation")
ggsave(sprintf("%s/%s/diag/rerun_discrepant.png", fpath, model), width = 5, height = 5)


g <- ggplot(dfm) +
    aes(init, abs(new_discrepancy)) +
    geom_boxplot() +
    geom_jitter(width = 0.25)
ggsave(sprintf("%s/%s/diag/init_discrepant.png", fpath, model), width = 5, height = 5)


md_df <- dfm %>%
    group_by(gene, snp) %>%
    summarise(
        new_mabs_disc = mean(abs(new_discrepancy)),
        mabs_disc = mean(abs(discrepancy)),
        .groups = "drop"
    )

g <- ggplot(md_df) +
    aes(mabs_disc, new_mabs_disc) +
    geom_abline(slope = 1, intercept = 0, linetype =  "dashed") +
    geom_point() +
    labs(x = "New mean absolute discrepancy", y = "Original absolute discrepancy")
ggsave(sprintf("%s/%s/diag/comp-discrepancy.png", fpath, model), width = 5, height = 5)





# dfma <- merge(df_new, df, suffix = c("new", "old"))
# dfma <- dfma[order(dfma$hmc), ]
# dfma$ind <- 1:nrow(dfma)

# g <- ggplot(dfma) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = vb_low,
#             y = vb,
#             ymax = vb_high,
#             colour = "vb"
#         ),
#         position = position_nudge(x = 0.25)
#     ) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = `5.0%`,
#             y = mean,
#             ymax = `95.0%`,
#             colour = "vb_rerun"
#         ),
#         position = position_nudge(x = -0.25)
#     ) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = hmc_low,
#             y = hmc,
#             ymax = hmc_high,
#             colour = "hmc"
#         )
#     )
# ggsave(sprintf("%s/%s/diag/discrepant_posteriors_out.png", fpath, model), width = 20, height = 7)


# dfmma <- dfma[abs(dfma$discrepancy) != max(abs(dfma$discrepancy)), ]

# g <- ggplot(dfmma) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = vb_low,
#             y = vb,
#             ymax = vb_high,
#             colour = "vb"
#         ),
#         position = position_nudge(x = 0.25)
#     ) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = `5.0%`,
#             y = mean,
#             ymax = `95.0%`,
#             colour = "vb_rerun"
#         ),
#         position = position_nudge(x = -0.25)
#     ) +
#     geom_pointrange(
#         aes(x = factor(ind),
#             ymin = hmc_low,
#             y = hmc,
#             ymax = hmc_high,
#             colour = "hmc"
#         )
#     )
# ggsave(sprintf("%s/%s/diag/discrepant_posteriors.png", fpath, model), width = 20, height = 7)



# df_lots <- lapply(1:5, function(j) {
#     new_res <- parallel::mclapply(1:nrow(df), function(i) {
#         fit_stan(df[i, "gene"], df[i, "snp"])
#     }, mc.cores = 5)
#     df_new <- do.call(rbind, new_res)
#     dfm <- merge(df_new[, c("mean", "khat", "time", "gene", "snp")], df)
#     dfm <- dfm %>% mutate(new_discrepancy = hmc - mean)
#     dfm
# })
