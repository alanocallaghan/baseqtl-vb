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
    type = "character"
)
args <- parser$parse_args()

if (is.null(method <- args[["method"]])) {
    method <- "vb"
}
mname <- switch(method,
    "vb" = "ADVI",
    "optimizing" = "MAP",
    "sampling" = "HMC"
)

dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

files <- list.files(dir)
files <- grep("rbias", files, value=TRUE)
stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))

dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"


outfiles <- list.files(sprintf("rds/GT/%s/", method), pattern = "ENSG*", full.names = TRUE)
infiles <- sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, genes)


dfs <- lapply(
    1:length(infiles),
    function(i) {
        cat(i, "/", length(infiles), "\n")
        infile <- infiles[[i]]
        outfile <- outfiles[[i]]
        out <- readRDS(outfile)
        inp <- readRDS(infile)
        if (!length(out)) {
            return(list())
        }
        covars <- lapply(inp,
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
        df <- do.call(rbind, out)
        df <- cbind(covars, df)
        df$snp <- names(out)
        df
    }
    # ,
    # infiles,
    # outfiles,
    # SIMPLIFY = FALSE
)


for (i in 1:length(dfs)) {
    if (length(dfs[[i]])) {
        dfs[[i]]$gene <- gsub(".*(ENSG\\d+).*", "\\1", outfiles[[i]])
    }
}
approx_res_df <- do.call(rbind, dfs)

# all_stan_res <- readRDS(sprintf("rds/GT/%s/all.rds", method))
# table(as.logical(sapply(all_stan_res, length)))
# FALSE  TRUE 
#     6   100 

# approx_res_df <- do.call(
#     rbind,
#     lapply(
#         all_stan_res,
#         function(x) {
#             if (is.null(x)) x <- list()
#             df <- do.call(rbind, x)
#             df$rSNP <- names(x)
#             df
#         }
#     )
# )

# rownames(approx_res_df) <- NULL
# if (is.null(approx_res_df$gene)) approx_res_df$gene <- approx_res_df$feature

sample_res_df <- do.call(
    rbind,
    lapply(genes,
        function(gene) {
            read.table(sprintf("%s/rbias.%s.stan.summary.txt", dir, gene), header=TRUE)
        }
    )
)
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


for (x in c("gene", "time", "n_ase", "mean_count", "sd_count", "n_wt", "n_het", "n_hom")) {
    scale <- if (x == "gene") scale_colour_discrete(guide="none") else scale_colour_viridis()
    g <- ggplot(mdf) +
        aes_string("mean", "log2_aFC_mean", colour = x) +
        geom_point(size = 0.8) +
        scale +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        lims(x = r, y = r) +
        labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
    ggsave(sprintf("fig/GT/diag/%s_%s.png", x, method))
    
    g <- ggplot(mdf) +
        aes_string(x, "discrepancy") +
        geom_point(size = 0.8) +
        # scale_colour_viridis() +
        # geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        # lims(x = r, y = r) +
        geom_smooth() +
        labs(x = x, y = "Discrepancy")

    ggsave(sprintf("fig/GT/diag/disc_%s_%s.png", x, method))
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
    geom_smooth() +
    theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Gene", y = "Discrepancy")

ggsave(sprintf("fig/GT/diag/disc_gene_box_%s.png", method))


# stop()



g <- ggplot() +
    aes(approx_res_df$mean, sample_res_df$log2_aFC_mean) +
    geom_pointdensity() +
    scale_colour_viridis() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_all.png", method), width = 7, height = 7)


mdf <- mdf[order(mdf$Null.99), ]

g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.99) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_all_categorical_99.png", method), width = 7, height = 7)

table(ADVI = mdf$null.99.VB, HMC = mdf$null.99.HMC)
## map
#      HMC
# ADVI     no   yes
#   no     56   304
#   yes    23 49677
## vb
#      HMC
# ADVI     no   yes
#   no     48   609
#   yes    13 39546

yardstick::sens_vec(factor(mdf$null.99.HMC), factor(mdf$null.99.VB))
## vb
# [1] 0.7868852
## map
# [1] 0.7088608
yardstick::spec_vec(factor(mdf$null.99.HMC), factor(mdf$null.99.VB))
## vb
# [1] 0.9848338
## map
# [1] 0.9939177


mdf <- mdf[order(mdf$Null.95), ]
g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.95) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_all_categorical_95.png", method), width = 7, height = 7)

table(ADVI = mdf$null.95.VB, HMC = mdf$null.95.HMC)
## map
#      HMC
# ADVI     no   yes
#   no    110   732
#   yes    61 49157
## vb
#      HMC
# ADVI     no   yes
#   no    113  1293
#   yes    21 38789


yardstick::sens_vec(factor(mdf$null.95.HMC), factor(mdf$null.95.VB))
## map
# [1] 0.6432749
## vb
# [1] 0.8432836
yardstick::spec_vec(factor(mdf$null.95.HMC), factor(mdf$null.95.VB))


yardstick::sens_vec(factor(mdf$null.99.HMC), factor(mdf$null.95.VB))
## map
# [1] 0.7594937
## vb
# [1] 0.9672131
yardstick::spec_vec(factor(mdf$null.99.HMC), factor(mdf$null.95.VB))


mdf <- mdf[order(mdf$Null.95.50), ]
g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.95.50) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_all_categorical_50.png", method), width = 7, height = 7)

table(ADVI = mdf$null.50.VB, HMC = mdf$null.50.HMC)
## vb
#      HMC
# ADVI     no   yes
#   no   2175 13805
#   yes  1859 22377
## map
#      HMC
# ADVI     no   yes
#   no   1375  1935
#   yes  3375 43375
yardstick::sens_vec(factor(mdf$null.50.HMC), factor(mdf$null.50.VB))
## vb
# [1] 0.5391671
## map
# [1] 0.2894737
yardstick::spec_vec(factor(mdf$null.50.HMC), factor(mdf$null.50.VB))
## map
# [1] 0.9572942
## vb
# [1] 0.6184567



yardstick::sens_vec(factor(mdf$null.95.HMC), factor(mdf$null.50.VB))
## map
# [1] 0.9239766
## vb
# [1] 0.9850746

yardstick::spec_vec(factor(mdf$null.95.HMC), factor(mdf$null.50.VB))
## map
# [1] 0.9368197
## vb
# [1] 0.6046105


yardstick::sens_vec(factor(mdf$null.99.HMC), factor(mdf$null.50.VB))
## map
# [1] 0.9367089
## vb
# [1] 1
yardstick::spec_vec(factor(mdf$null.99.HMC), factor(mdf$null.50.VB))
## map
# [1] 0.9352554
## vb
# [1] 0.6035612





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
yardstick::sens(summ, anynull99hmc, anynull99vb)
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.293
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.576
yardstick::spec(summ, anynull99hmc, anynull99vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary         0.929
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1

yardstick::sens(summ, anynull95hmc, anynull95vb)
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary        0.0822
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary        0.0822
yardstick::spec(summ, anynull95hmc, anynull95vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1


yardstick::sens(summ, anynull99hmc, anynull95vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.207
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.0732
yardstick::spec(summ, anynull99hmc, anynull95vb)
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1



yardstick::sens(summ, anynull50hmc, anynull50vb)
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary             0
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.167
yardstick::spec(summ, anynull50hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary          0.98
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary          0.98



yardstick::sens(summ, anynull99hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary        0.0326
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary        0.0326
yardstick::spec(summ, anynull99hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1



yardstick::sens(summ, anynull95hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary        0.0370
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary             0
yardstick::spec(summ, anynull95hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1




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
ggsave(sprintf("fig/GT/time/time_dist_all_%s.png", method), width = 6, height = 6)



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

g <- ggplot() +
    aes(approx_res_df_sub$mean, ss$log2_aFC_mean) +
    geom_pointdensity() +
    scale_colour_viridis() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r2, y = r2) +
    labs(x = sprintf("%s estimate", method), y = "MCMC estimate")
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc.png", method), width = 7, height = 7)

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
    aes(mean, log2_aFC_mean, colour = Null.99) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r2, y = r2) +
    labs(x = sprintf("%s estimate", method), y = "MCMC estimate")

# ggMarginal(g)
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_categorical_99.png", method), width = 7, height = 7)

yardstick::sens_vec(factor(mdf2$null.99.HMC), factor(mdf2$null.99.VB))
## vb
# [1] 0.8
## map
# [1] 0.5294118
yardstick::spec_vec(factor(mdf2$null.99.HMC), factor(mdf2$null.99.VB))
## map
# [1] 0.9553957


table(ADVI = mdf2$null.99.VB, HMC = mdf2$null.99.HMC)
## MAP
#        HMC
# ADVI     no yes
#   FALSE  11 421
#   TRUE    6 202
## vb
#      HMC
# ADVI    no  yes
#   no     8   80
#   yes    2 1002



mdf2 <- mdf2[order(mdf2$Null.95), ]
g <- ggplot(mdf2) +
    aes(mean, log2_aFC_mean, colour = Null.95) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r2, y = r2) +
    labs(x = sprintf("%s estimate", method), y = "MCMC estimate")

# ggMarginal(g)
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_categorical_95.png", method), width = 7, height = 7)

yardstick::sens_vec(factor(mdf2$null.95.HMC), factor(mdf2$null.95.VB))
## vb
# [1] 0.7931034
# map
# [1] 0.325
yardstick::spec_vec(factor(mdf2$null.95.HMC), factor(mdf2$null.95.VB))
## map
# [1] 0.9556962

yardstick::sens_vec(factor(mdf2$null.99.HMC), factor(mdf2$null.95.VB))
## map
# [1] 0.5882353
yardstick::spec_vec(factor(mdf2$null.99.HMC), factor(mdf2$null.95.VB))
## map
# [1] 0.9366906


table(ADVI = mdf2$null.95.VB, HMC = mdf2$null.95.HMC)
## MAP
#      HMC
# ADVI    no  yes
#   no    52   56
#   yes  108 1208
## vb
#      HMC
# ADVI   no yes
#   no   46 102
#   yes  12 932


mdf2 <- mdf2[order(mdf2$Null.95.50), ]
g <- ggplot(mdf2) +
    aes(mean, log2_aFC_mean, colour = Null.95.50) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r2, y = r2) +
    labs(x = sprintf("%s estimate", method), y = "MCMC estimate")

# ggMarginal(g)
ggsave(file = sprintf("fig/GT/estimates/%s_mcmc_categorical_50.png", method), width = 7, height = 7)


table(ADVI = mdf2$null.50.VB, HMC = mdf2$null.50.HMC)
## vb
#      HMC
# ADVI   no yes
#   no  254 298
#   yes 148 392
## map
#      HMC
# ADVI   no yes
#   no  390 122
#   yes 252 660

yardstick::sens_vec(factor(mdf2$null.50.HMC), factor(mdf2$null.50.VB))
## vb
# [1] 0.6318408
## map
# [1] 0.5294118
yardstick::spec_vec(factor(mdf2$null.50.HMC), factor(mdf2$null.50.VB))
## map
# [1] 0.8439898

yardstick::sens_vec(factor(mdf2$null.99.HMC), factor(mdf2$null.50.VB))
## map
# [1] 1
yardstick::spec_vec(factor(mdf2$null.99.HMC), factor(mdf2$null.50.VB))
## map
# [1] 0.6561151
yardstick::sens_vec(factor(mdf2$null.95.HMC), factor(mdf2$null.50.VB))
## map
# [1] 0.925
yardstick::spec_vec(factor(mdf2$null.95.HMC), factor(mdf2$null.50.VB))
## map
# [1] 0.7120253



summ2 <- mdf2 %>%
    group_by(gene) %>%
    summarise(
        anynull99vb = any(null.99.VB == "no"),
        anynull99hmc = any(null.99.HMC == "no"),
        anynull95vb = any(null.95.VB == "no"),
        anynull95hmc = any(null.95.HMC == "no"),
        anynull50vb = any(null.50.VB == "no"),
        anynull50hmc = any(null.50.HMC == "no")
    )
summ2[-1] <- lapply(summ2[-1], factor)
yardstick::sens(summ2, anynull99hmc, anynull99vb)
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.826
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.870
yardstick::spec(summ2, anynull99hmc, anynull99vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary          0.75


yardstick::sens(summ2, anynull95hmc, anynull95vb)
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.762
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.878
yardstick::spec(summ2, anynull95hmc, anynull95vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary         0.778

yardstick::sens(summ2, anynull99hmc, anynull95vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.833
yardstick::spec(summ2, anynull99hmc, anynull95vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>

yardstick::sens(summ2, anynull50hmc, anynull50vb)
## vb
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary           0.1
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.619
yardstick::spec(summ2, anynull50hmc, anynull50vb)


yardstick::sens(summ2, anynull99hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.333
yardstick::spec(summ2, anynull99hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1


yardstick::sens(summ2, anynull95hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 sens    binary         0.367
yardstick::spec(summ2, anynull95hmc, anynull50vb)
## map
# # A tibble: 1 × 3
#   .metric .estimator .estimate
#   <chr>   <chr>          <dbl>
# 1 spec    binary             1



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
ggsave(sprintf("fig/GT/time/time_dist_all_%s.png", method), width = 6, height = 6)
