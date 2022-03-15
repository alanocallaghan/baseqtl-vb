library("argparse")
library("baseqtl")
library("ggplot2")
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
tol <- args[["tolerance"]]
method <- args[["method"]]

mname <- switch(method,
    "vb" = "ADVI",
    "optimizing" = "MAP",
    "sampling" = "HMC"
)
mtol <- if (method == "vb") sprintf("vb_%1.0e", tol) else method

dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"

# files <- list.files(dir)
# files <- grep("refbias", files, value = TRUE)
outfiles <- list.files(sprintf("rds/noGT/%s/", mtol), pattern = "ENSG*", full.names=TRUE)
genes <- unique(gsub(".*(ENSG\\d+)..*", "\\1", outfiles))

infiles <- list(
    normal_skin = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, genes),
    Psoriasis_skin = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, genes)
)
i <- 1

process_nogt <- function(x) {
    inp1 <- in.neg.beta.noGT.eff2(x)
    data.frame(
        n_tot = inp1$N,
        n_ase = inp1$A,
        mean_ase = mean(log1p(inp1$m)),
        sd_ase = sd(log1p(inp1$m)),
        mean_count = mean(log1p(inp1$Y)),
        sd_count = sd(log1p(inp1$Y)),
        n_wt = sum(inp1$gase == 0),
        n_het = sum(abs(inp1$gase) == 1),
        n_hom = sum(abs(inp1$gase) == 2)
    )
}
dfs <- parallel::mclapply(
    1:length(genes),
    function(i) {
        cat(i, "/", length(genes), "\n")
        infile_norm <- infiles[[1]][[i]]
        infile_pso <- infiles[[2]][[i]]
        
        outfiles_norm <- list.files(
            sprintf("rds/noGT/%s/", mtol),
            pattern = sprintf("%s_.*_normal_skin.rds", genes[[i]]),
            full.names = TRUE
        )
        outfiles_pso <- list.files(
            sprintf("rds/noGT/%s/", mtol),
            pattern = sprintf("%s_.*_Psoriasis_skin.rds", genes[[i]]),
            full.names = TRUE
        )
        outs_norm <- do.call(rbind, lapply(outfiles_norm, readRDS))
        outs_pso <- do.call(rbind, lapply(outfiles_pso, readRDS))
        out <- rbind(outs_norm, outs_pso)
        cn <- setdiff(
            colnames(out),
            c("n_eff", "Rhat", "null.99", "gene", "time", "snp", "condition")
        )
        out[, cn] <- out[, cn] / log(2)

        inp_norm <- readRDS(infile_norm)
        inp_pso <- readRDS(infile_pso)
        if (!length(out)) {
            return(list())
        }
        covars_norm <- lapply(inp_norm, process_nogt)
        covars_pso <- lapply(inp_pso, process_nogt)
        covars_norm <- do.call(rbind, covars_norm)
        covars_pso <- do.call(rbind, covars_pso)
        covars_norm$condition <- "normal_skin"
        covars_pso$condition <- "Psoriasis_skin"
        covars <- rbind(covars_norm, covars_pso)
        covars$gene <- genes[[i]]
        covars$snp <- rownames(covars)
        # df_norm <- do.call(rbind, out[[1]])
        # df_pso <- do.call(rbind, out[[2]])
        # df_all <- rbind(df_norm, df_pso)
        df_all <- merge(covars, out)
        df_all
    },
    mc.cores = 8
)



# all_stan_res <- lapply(
#     list.files(sprintf("rds/noGT/%s/", method), pattern = "ENSG*", full.names=TRUE),
#     readRDS
# )
# all_stan_res <- lapply(all_stan_res,
#     function(x) {
#         if (!length(x)) return(NULL)
#         c(x[[1]], x[[2]])
#     }
# )
# all_stan_res <- Reduce(c, all_stan_res)
# all_stan_res <- do.call(rbind, all_stan_res)
# approx_res_df <- as.data.frame(all_stan_res)


approx_res_df <- do.call(rbind, dfs)
approx_res_df$test <- paste(approx_res_df$gene, approx_res_df$snp , sep = "_")

summary_files <- list.files(dir, pattern = "stan.summary", full.names = TRUE)
sdata <- data.frame(
    tissue = gsub(".*(Psoriasis_skin|normal_skin).*", "\\1", summary_files)
)
ll <- lapply(seq_along(summary_files),
    function(i) {
        x <- read.table(summary_files[[i]], header = TRUE)
        x$condition <- sdata[i, ]
        x
    }
)
cn <- Reduce(intersect, lapply(ll, colnames))
ll <- lapply(ll, function(x) x[, cn])
sample_res_df <- do.call(rbind, ll)


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


mdf <- merge(approx_res_df, sample_res_df, by = c("test", "condition"), suffix = c(".VB", ".HMC"))
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

mdf <- mdf %>% mutate(allele_freq = (n_het + (n_hom * 2)) / ((n_wt + n_het + n_hom) * 2))


mdf2 <- mdf[mdf$Rhat.VB < 1.05 & mdf$Rhat.HMC < 1.05 & mdf$n_eff.VB > 500 & mdf$n_eff.HMC > 500, ]

g <- ggplot(mdf2) +
    aes_string("log2_aFC_mean", "mean", colour = "log2_aFC_se_mean") +
    geom_point(size = 0.8) +
    scale_colour_viridis(name = "SE (x)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "Old estimate", y = "New estimate")
ggsave("fig/noGT/estimates/se1.png", width = 7, height = 7)
g <- ggplot(mdf2) +
    aes_string("log2_aFC_mean", "mean", colour = "se_mean") +
    geom_point(size = 0.8) +
    scale_colour_viridis(name = "SE (y)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "Old estimate", y = "New estimate")
ggsave("fig/noGT/estimates/se2.png", width = 7, height = 7)

for (x in c("gene", "time", "allele_freq", "mean_ase", "sd_ase", "n_ase", "mean_count", "sd_count", "n_wt", "n_het", "n_hom")) {
    scale <- if (x == "gene") scale_colour_discrete(guide="none") else scale_colour_viridis()
    # g <- ggplot(mdf) +
    #     aes_string("mean", "log2_aFC_mean", colour = x) +
    #     geom_point(size = 0.8) +
    #     scale +
    #     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    #     lims(x = r, y = r) +
    #     labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
    # ggsave(sprintf("fig/noGT/diag/xy_%s_%s.png", x, mtol), width = 7, height = 7)
    
    g <- ggplot(mdf) +
        aes_string(x, "discrepancy") +
        geom_point(size = 0.8) +
        # scale_colour_viridis() +
        # geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        # lims(x = r, y = r) +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
        labs(x = x, y = "Discrepancy")

    ggsave(sprintf("fig/noGT/diag/disc_%s_%s.png", x, mtol), width = 7, height = 7)
}


tmp <- mdf %>%
    group_by(gene) %>%
    summarise(mean_disc = mean(discrepancy), sum_disc = sum(abs(discrepancy) > 0.3))
tmp %>% arrange(-sum_disc)
gg <- tmp %>% arrange(-sum_disc) %>% top_n(5, wt = sum_disc) %>% pull(gene)

g <- ggplot(mdf) +
    aes(gene, discrepancy, colour = gene %in% gg) +
    geom_violin() +
    # geom_jitter(size = 0.8) +
    # scale_colour_viridis() +
    # geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    # lims(x = r, y = r) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Gene", y = "Discrepancy")

ggsave(sprintf("fig/noGT/diag/disc_gene_box_%s.png", mtol), width = 7, height = 7)


# g <- ggplot() +
#     aes(approx_res_df$mean, sample_res_df$log2_aFC_mean) +
#     geom_pointdensity() +
#     scale_colour_viridis() +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     lims(x = r, y = r) +
#     labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
# ggsave(file = sprintf("fig/noGT/estimates/%s_mcmc_all.png", mtol), width = 7, height = 7)


mdf <- mdf[order(mdf$Null.99), ]

g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.99) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/noGT/estimates/%s_mcmc_all_categorical_99.png", mtol), width = 7, height = 7)



mdf <- mdf[order(mdf$Null.95), ]
g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.95) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/noGT/estimates/%s_mcmc_all_categorical_95.png", mtol), width = 7, height = 7)

mdf <- mdf[order(mdf$Null.95.50), ]
g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.95.50) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/noGT/estimates/%s_mcmc_all_categorical_50.png", mtol), width = 7, height = 7)


summ <- mdf %>%
    group_by(gene) %>%
    summarise(
        anynull99vb = any(null.99.VB == "no"),
        anynull99hmc = any(null.99.HMC == "no"),
        anynull95vb = any(null.95.VB == "no"),
        anynull95hmc = any(null.95.HMC == "no"),
        anynull50vb = any(null.50.VB == "no"),
        anynull50hmc = any(null.50.HMC == "no")
    )
summ[-1] <- lapply(summ[-1], factor, levels = c("FALSE", "TRUE"))

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
ggsave(sprintf("fig/noGT/time/time_dist_all_%s.png", mtol), width = 6, height = 6)

