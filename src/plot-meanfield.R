### comparing fullrank to meanfield approximations
### pick 1 association per gene
### compare posterior means and variances

library("baseqtl")
library("argparse")
library("ggdist")
library("dplyr")
library("tidyr")
library("transport")
library("rstan")
library("ggpointdensity")
library("viridis")

parser <- ArgumentParser()
parser$add_argument(
  "-m", "--model",
  default = "GT",
  type = "character"
)
parser$add_argument(
  "-t", "--tolerance",
  default = 1e-2,
  type = "double"
)

args <- parser$parse_args()

source("src/functions.R")
tol <- args[["tolerance"]]
model <- args[["model"]]

theme_set(theme_bw())

fpath <- sprintf("fig_%1.0e", tol)
source("src/functions.R")
mkfigdir(fpath, model)

files <- list.files(sprintf("rds/%s/vb_%1.0e/meanfield/", model, tol), full.names = TRUE)
res <- lapply(files, readRDS)
genes <- gsub(".*(ENSG\\d+)_.*.rds", "\\1", files)
names(res) <- genes

draws <- lapply(names(res), function(gene) {
  lapply(names(res[[gene]]), function(snp) {
    tab <- lapply(
      names(res[[gene]][[snp]]),
      function(method) {
        x <- res[[gene]][[snp]][[method]]
        if (model == "noGT") {
          out <- lapply(x, function(y) extract(y, pars = "bj")$bj)
          names(out) <- names(x)
        } else {
          out <- extract(x, pars = "bj")$bj
        }
        out
      }
    )
    tab <- as.data.frame(tab, check.names = FALSE)
    if (model == "noGT") {
      colnames(tab) <- paste(names(tab), rep(names(res[[gene]][[snp]]), each = 2))
    } else {
      colnames(tab) <- names(res[[gene]][[snp]])
    }
    tab$snp <- snp
    tab$gene <- gene
    tab
  })
})

if (model == "noGT") {
  draw_dfs <- lapply(draws, function(x) do.call(rbind, x))
  draws_df <- do.call(rbind, draw_dfs)
  mdf <- pivot_longer(draws_df,
    `normal_skin meanfield`:`Psoriasis_skin hmc`,
    names_sep = " ", names_to = c("condition", "method")
  )
  g <- ggplot(mdf[mdf$gene %in% genes[1:5], ]) +
    aes(x = snp, colour = Method, fill = Method, y = value) +
    stat_pointinterval(position = position_dodge()) +
    facet_grid(cols = vars(gene), rows = vars(condition), scale = "free") +
    scale_colour_brewer(palette = "Paired", aesthetics = c("fill", "colour")) +
    theme_bw()
  ggsave("tmp.png", width = 12, height = 5)

  em_df <- mdf |>
    pivot_wider(
      names_from = method,
      values_from = value
    ) |>
    group_by(gene, snp, condition) |>
    summarise(
      # kl_mf = KL.empirical(meanfield[[1]], hmc[[1]]),
      # kl_fr = KL.empirical(fullrank[[1]], hmc[[1]]),
      ## KL is messy, would need to compute approximate densities
      ks_mf = ks.test(meanfield[[1]], hmc[[1]], exact = FALSE)$statistic,
      ks_fr = ks.test(fullrank[[1]], hmc[[1]], exact = FALSE)$statistic,
      wass_mf = wasserstein1d(meanfield[[1]], hmc[[1]]),
      wass_fr = wasserstein1d(fullrank[[1]], hmc[[1]]),
      .groups = "drop_last"
    )
  g <- ggplot(em_df) +
    aes(wass_mf, wass_fr) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_pointdensity() +
    # geom_point() +
    scale_colour_viridis(guide = "none") +
    labs(x = "Wassertein(mean field VB | HMC)", y = "Wassertein(full rank VB | HMC)")
  # ggsave("was2.png")
  ggsave(
    sprintf("%s/%s/diag/wasserstein-meanfield-fullrank.pdf", fpath, model),
    width = 5, height = 5
  )

  g <- ggplot(em_df) +
    aes(ks_mf, ks_fr) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_pointdensity() +
    # geom_point() +
    scale_colour_viridis(guide = "none") +
    labs(x = "KS(mean field VB || HMC)", y = "KS(full rank VB || HMC)")
  # ggsave("was2.png")
  ggsave(
    sprintf("%s/%s/diag/KS-meanfield-fullrank.pdf", fpath, model),
    width = 5, height = 5
  )
} else {
  draw_dfs <- lapply(draws, function(x) do.call(rbind, x))
  draws_df <- do.call(rbind, draw_dfs)
  mdf <- tidyr::pivot_longer(draws_df, `meanfield`:`hmc`)

  # g <- ggplot(mdf[mdf$gene %in% genes[1:6], ]) +
  #     aes(x = snp, colour = name, fill = name, y = value) +
  #     stat_pointinterval(position = position_dodge()) +
  #     facet_wrap(~gene, scale = "free") +
  #     scale_colour_brewer(palette = "Paired", aesthetics = c("fill", "colour")) +
  #     theme_bw()
  # ggsave("tmp.png",
  #     width = 12, height = 5
  # )

  em_df <- draws_df |>
    group_by(gene, snp) |>
    summarise(
      ks_mf = ks.test(meanfield, hmc, exact = FALSE)$statistic,
      ks_fr = ks.test(fullrank, hmc, exact = FALSE)$statistic,
      wass_mf = wasserstein1d(meanfield, hmc),
      wass_fr = wasserstein1d(fullrank, hmc),
      .groups = "drop_last"
    )
  g1 <- ggplot(em_df) +
    aes(wass_mf, wass_fr) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_pointdensity() +
    # geom_point() +
    scale_colour_viridis(guide = "none") +
    labs(x = "Wassertein(mean field VB | HMC)", y = "Wassertein(full rank VB | HMC)")
  ggsave(
    sprintf("%s/%s/diag/wasserstein-meanfield-fullrank.pdf", fpath, model),
    width = 4, height = 4
  )

  g2 <- ggplot(em_df) +
    aes(ks_mf, ks_fr) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_pointdensity() +
    # geom_point() +
    scale_colour_viridis(guide = "none") +
    labs(x = "KS(mean field VB | HMC)", y = "KS(full rank VB | HMC)")
  ggsave(
    sprintf("%s/%s/diag/KS-meanfield-fullrank.pdf", fpath, model),
    width = 4, height = 4
  )
  gg <- cowplot::plot_grid(g1, g2, labels = "AUTO")
  ggsave(
    sprintf("%s/%s/diag/meanfield-fullrank.pdf", fpath, model),
    width = 5.5, height = 3
  )
}
