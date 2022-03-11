library("argparse")
library("ggplot2")
library("baseqtl")
library("ggpointdensity")
library("ggExtra")
library("viridis")
library("dplyr")
theme_set(theme_bw())

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--method",
    default = "sampling",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-2,
    type = "double"
)

args <- parser$parse_args()
method <- args[["method"]]
tol <- args[["tolerance"]]

mname <- switch(method,
    "vb" = "ADVI",
    "optimizing" = "MAP",
    "sampling" = "HMC"
)

mtol <- if (method == "vb") sprintf("vb_%1.0e", tol) else method

dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

files <- list.files(dir)
files <- grep("rbias", files, value=TRUE)
stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))

dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

outfiles <- list.files(sprintf("rds/GT/%s/", mtol), pattern = "ENSG*", full.names = TRUE)
genes <- unique(gsub(".*(ENSG\\d+).*", "\\1", outfiles))
infiles <- sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, genes)

dfs <- parallel::mclapply(
    1:length(genes),
    function(i) {
        cat(i, "/", length(genes), "\n")
        gene <- genes[[i]]
        infile <- infiles[[i]]
        outfiles <- list.files(
            sprintf("rds/GT/%s/", mtol),
            pattern = paste0(gene, ".*"),
            full.names = TRUE
        )
        # outfile <- outfiles[[i]]
        out <- do.call(rbind, lapply(outfiles, readRDS))
        inp <- readRDS(infile)
        if (!length(out)) {
            return(list())
        }
        snp <- out$snp
        covars <- lapply(inp[snp],
            function(x) {
                inp1 <- in.neg.beta.prob.eff2(x)
                data.frame(
                    n_tot = inp1$N,
                    n_ase = inp1$A,
                    mean_count = mean(log1p(inp1$Y)),
                    sd_count = sd(log1p(inp1$Y)),
                    n_wt = sum(inp1$g == 0),
                    n_het = sum(abs(inp1$g) == 1),
                    n_hom = sum(abs(inp1$g) == 2)
                )
            }
        )
        covars <- do.call(rbind, covars)
        # df <- do.call(rbind, out)
        df <- out
        df <- cbind(covars, df)
        df$snp <- snp
        df
    }, mc.cores = 8
)

approx_res_df <- do.call(rbind, dfs)

sample_res_df <- do.call(
    rbind,
    lapply(genes,
        function(gene) {
            read.table(sprintf("%s/rbias.%s.stan.summary.txt", dir, gene), header=TRUE)
        }
    )
)

sample_res_df <- sample_res_df[sample_res_df$n_eff > 500 & sample_res_df$Rhat < 1.05, ]
if (method == "sampling") {
    approx_res_df <- approx_res_df[approx_res_df$n_eff > 500 & approx_res_df$Rhat < 1.05, ]
}

approx_res_df$test <- paste(approx_res_df$gene, approx_res_df$snp, sep = "_")
# approx_res_df$test <- paste(approx_res_df$gene, approx_res_df$rSNP, sep = "_")
sample_res_df$test <- paste(sample_res_df$Gene_id, sample_res_df$tag , sep = "_")


approx_res_df$null.95 <- sign(approx_res_df$"2.5%") == sign(approx_res_df$"97.5%")
approx_res_df$null.50 <- sign(approx_res_df$"25.0%") == sign(approx_res_df$"75.0%")
sample_res_df$null.95 <- sign(sample_res_df$"log2_aFC_2.5.") == sign(sample_res_df$"log2_aFC_97.5.")
sample_res_df$null.50 <- sign(sample_res_df$"log2_aFC_25.") == sign(sample_res_df$"log2_aFC_75.")


sample_res_df <- sample_res_df[sample_res_df$test %in% approx_res_df$test, ]
im <- match(sample_res_df$test, approx_res_df$test)
approx_res_df <- approx_res_df[im, ]
stopifnot(all(sample_res_df$test == approx_res_df$test))

r <- range(c(sample_res_df$log2_aFC_mean, approx_res_df$mean))
mdf <- merge(approx_res_df, sample_res_df, by = "test", suffix = c(".VB", ".HMC"))
mdf$null.99.VB <- ifelse(mdf$null.99.VB, "no", "yes")
mdf$null.95.VB <- ifelse(mdf$null.95.VB, "no", "yes")
mdf$null.50.VB <- ifelse(mdf$null.50.VB, "no", "yes")
mdf$null.95.HMC <- ifelse(mdf$null.95.HMC, "no", "yes")
mdf$null.50.HMC <- ifelse(mdf$null.50.HMC, "no", "yes")

lab_str <- paste0(mname, ": %s\nHMC: %s\n")
levs <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "no",  "yes", "yes")
)
# levs <- rev(levs)
flev <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "yes",  "no", "yes")
)

mdf$Null.99 <- sprintf(lab_str, mdf$null.99.VB, mdf$null.99.HMC)
mdf$Null.95 <- sprintf(lab_str, mdf$null.95.VB, mdf$null.95.HMC)
mdf$Null.50 <- sprintf(lab_str, mdf$null.50.VB, mdf$null.50.HMC)
mdf$Null.95.50 <- sprintf(lab_str, mdf$null.50.VB, mdf$null.95.HMC)

mdf$Null.99 <- factor(mdf$Null.99, levels = flev)
mdf$Null.95 <- factor(mdf$Null.95, levels = flev)
mdf$Null.50 <- factor(mdf$Null.50, levels = flev)
mdf$Null.95.50 <- factor(mdf$Null.95.50, levels = flev)

mdf$discrepancy <- mdf$mean - mdf$log2_aFC_mean

# g <- ggplot(mdf) +
#     aes_string("mean", "mean_count") +
#     geom_point(size = 0.8)
# ggsave("tmp1.png", width = 7, height = 7)

## not se mean but other hmc se
diag_vars <- c(
    "gene",
    # "Rhat",
    # "n_eff",
    "time", "n_ase",
    "log2_aFC_se_mean", "log2_aFC_sd",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)

mdf2 <- mdf[mdf$Rhat.VB < 1.05 & mdf$Rhat.HMC < 1.05 & mdf$n_eff.VB > 500 & mdf$n_eff.HMC > 500, ]

g <- ggplot(mdf2) +
    aes_string("log2_aFC_mean", "mean", colour = "log2_aFC_se_mean") +
    geom_point(size = 0.8) +
    scale_colour_viridis(name = "SE (x)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "Old estimate", y = "New estimate")
ggsave("fig/GT/estimates/se1.png", width = 7, height = 7)
g <- ggplot(mdf2) +
    aes_string("log2_aFC_mean", "mean", colour = "se_mean") +
    geom_point(size = 0.8) +
    scale_colour_viridis(name = "SE (y)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "Old estimate", y = "New estimate")
ggsave("fig/GT/estimates/se2.png", width = 7, height = 7)
stop()


for (x in diag_vars) {
    scale <- if (x == "gene") scale_colour_discrete(guide="none") else scale_colour_viridis()
    
    # g <- ggplot(mdf[order(mdf[[x]]), ]) +
    #     aes_string("log2_aFC_mean", "mean", colour = x) +
    #     geom_point(size = 0.8) +
    #     scale +
    #     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    #     lims(x = r, y = r) +
    #     labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
    # ggsave(sprintf("fig/GT/diag/%s_%s.png", x, mtol), width = 7, height = 7)
    
    g <- ggplot(mdf) +
        aes_string(x, "abs(discrepancy)") +
        geom_point(size = 0.8) +
        # scale_colour_viridis() +
        # geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        # lims(x = r, y = r) +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")
        ) +
        labs(x = x, y = "Discrepancy")

    ggsave(sprintf("fig/GT/diag/disc_%s_%s.png", x, mtol), width = 7, height = 7)
}


mdf$gf <- as.numeric(factor(mdf$gene))

tmp <- mdf %>%
    group_by(gene) %>%
    summarise(mean_disc = mean(discrepancy), sum_disc = sum(abs(discrepancy) > 0.3))
tmp %>% arrange(-sum_disc)
gg <- tmp %>% arrange(-sum_disc) %>% top_n(5) %>% pull(gene)

g <- ggplot(mdf) +
    aes(gene, discrepancy, colour = gene %in% gg) +
    geom_violin() +
    # geom_jitter(size = 0.8) +
    # scale_colour_viridis() +
    # geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    # lims(x = r, y = r) +
    geom_smooth(method = "loess", formula = y ~ x) +
    theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Gene", y = "Discrepancy")

ggsave(sprintf("fig/GT/diag/disc_gene_box_%s.png", mtol), width=7, height = 7)


# g <- ggplot() +
#     aes(sample_res_df$log2_aFC_mean, approx_res_df$mean) +
#     geom_pointdensity() +
#     scale_colour_viridis() +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     lims(x = r, y = r) +
#     labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
# ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_all.png", mtol), width = 7, height = 7)


mdf <- mdf[order(mdf$Null.99), ]

g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.99) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_all_categorical_99.png", mtol), width = 7, height = 7)


mdf <- mdf[order(mdf$Null.95), ]
g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.95) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_all_categorical_95.png", mtol), width = 7, height = 7)

mdf <- mdf[order(mdf$Null.95.50), ]
g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.95.50) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_all_categorical_50.png", mtol), width = 7, height = 7)


ggplot() +
    aes(mdf$time) +
    geom_histogram(
        bins = nclass.FD(mdf$time),
        colour = "grey60",
        fill = "grey90",
        boundary = 0
    ) +
    geom_vline(xintercept = median(mdf$time), linetype = "dashed") +
    labs(x = "Time (s)", y = "Frequency")
ggsave(sprintf("fig/GT/time/time_dist_all_%s.png", mtol), width = 6, height = 6)



sample_data <- readRDS(
    "/rds/project/cew54/rds-cew54-wallace-share/Projects/baseqtl/data/btrecase.GT.nGT.all.rds"
)
ss <- as.data.frame(sample_data$RNA)
ss$test <- paste(ss$Gene_id, ss$tag , sep = "_")
ss$null.50 <- sign(ss$"log2_aFC_75%") == sign(ss$"log2_aFC_25%")

ss <- ss[ss$test %in% approx_res_df$test, ]
im <- match(ss$test, approx_res_df$test)
approx_res_df_sub <- approx_res_df[im, ]
approx_res_df_sub$null.50 <- sign(approx_res_df_sub$"75.0%") == sign(approx_res_df_sub$"25.0%")
stopifnot(all(ss$test == approx_res_df_sub$test))
r2 <- range(c(ss$log2_aFC_mean, approx_res_df_sub$mean))

# g <- ggplot() +
#     aes(ss$log2_aFC_mean, approx_res_df_sub$mean) +
#     geom_pointdensity() +
#     scale_colour_viridis() +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     lims(x = r2, y = r2) +
#     labs(x = "MCMC estimate", y = sprintf("%s estimate", method))
# ggsave(file = sprintf("fig/GT/estimates/%s_mcmc.png", mtol), width = 7, height = 7)

mdf2 <- merge(approx_res_df_sub, ss, by = "test", suffix = c(".VB", ".HMC"))
mdf2$null.99.VB <- ifelse(mdf2$null.99.VB, "no", "yes")
mdf2$null.95.VB <- ifelse(mdf2$null.95.VB, "no", "yes")
mdf2$null.50.VB <- ifelse(mdf2$null.50.VB, "no", "yes")
mdf2$null.50.HMC <- ifelse(mdf2$null.50.HMC, "no", "yes")


mdf2$Null.99 <- sprintf(lab_str, mdf2$null.99.VB, mdf2$null.99.HMC)
mdf2$Null.95 <- sprintf(lab_str, mdf2$null.95.VB, mdf2$null.95.HMC)
mdf2$Null.50 <- sprintf(lab_str, mdf2$null.50.VB, mdf2$null.50.HMC)
mdf2$Null.95.50 <- sprintf(lab_str, mdf2$null.50.VB, mdf2$null.95.HMC)

mdf2$Null.99 <- factor(mdf2$Null.99, levels = flev)
mdf2$Null.95 <- factor(mdf2$Null.95, levels = flev)
mdf2$Null.50 <- factor(mdf2$Null.50, levels = flev)
mdf2$Null.95.50 <- factor(mdf2$Null.95.50, levels = flev)


mdf2 <- mdf2[order(mdf2$Null.99), ]
g <- ggplot(mdf2) +
    aes(log2_aFC_mean, mean, colour = Null.99) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r2, y = r2) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", method))

# ggMarginal(g)
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_categorical_99.png", mtol), width = 7, height = 7)


mdf2 <- mdf2[order(mdf2$Null.95), ]
g <- ggplot(mdf2) +
    aes(log2_aFC_mean, mean, colour = Null.95) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r2, y = r2) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", method))

# ggMarginal(g)
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_categorical_95.png", mtol), width = 7, height = 7)


mdf2 <- mdf2[order(mdf2$Null.95.50), ]
g <- ggplot(mdf2) +
    aes(log2_aFC_mean, mean, colour = Null.95.50) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r2, y = r2) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", method))

# ggMarginal(g)
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_categorical_50.png", mtol), width = 7, height = 7)


ggplot() +
    aes(mdf2$time) +
    geom_histogram(
        bins = nclass.FD(mdf2$time),
        colour = "grey60",
        fill = "grey90",
        boundary = 0
    ) +
    geom_vline(xintercept = median(mdf2$time), linetype = "dashed") +
    labs(x = "Time (s)", y = "Frequency")
ggsave(sprintf("fig/GT/time/time_dist_all_%s.png", mtol), width = 6, height = 6)
