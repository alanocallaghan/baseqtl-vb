library("ggplot2")
library("ggdist")
library("ggpointdensity")
library("baseqtl")
library("viridis")
library("yardstick")
library("geomtextpath")
library("dplyr")
library("argparse")

theme_set(theme_bw())

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "noGT",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-3,
    type = "double"
)

maxRhat <- 1.1 ## from baseqtl-paper repo
minEff <- 500 ## from stan docs (-ish)

args <- parser$parse_args()
tol <- args[["tolerance"]]

model <- args[["model"]]

# fpath <- sprintf("fig_%1.0e", tol)
fpath <- "fig"

# "optimizing",
methods <- c("vb", "sampling")
dfs <- lapply(methods,
    function(method) {
        cat(method, "\n")
        # mtol <- if (method == "vb") sprintf("%s_%1.0e", method, tol) else method
        mtol <- method
        path <- file.path("rds", model, mtol)
        f <- list.files(path = path, pattern = "ENSG.*.rds", full.names = TRUE)
        input <- lapply(f, readRDS)
        if (model == "noGT") {
            df <- do.call(rbind, input)
            df <- as.data.frame(df)
            df$test <- paste(df$gene, df$snp , sep = "_")
        } else {
            df <- do.call(rbind, input)
            df <- as.data.frame(df)
            df$test <- paste(df$gene, df$snp , sep = "_")
        }
        # if (method == "sampling") {
        #     df <- df[df$n_eff > 500 & df$Rhat < 1.05, ]
        # } else if (method == "vb") {
        #     ## from PSIS paper, arxiv 1507.02646
        #     # df <- df[df$khat < 0.7, ]
        # }
        
        df$method <- gsub("old/", "", method)
        ## TRUE means that the HPD interval is one side of zero (sig)
        ## FALSE means it overlaps zero (null)
        df$null.95 <- sign(df$"2.5%") == sign(df$"97.5%")
        df$null.90 <- sign(df$"5.0%") == sign(df$"95.0%")
        df$null.80 <- sign(df$"10.0%") == sign(df$"90.0%")
        df$null.70 <- sign(df$"15.0%") == sign(df$"85.0%")
        df$null.60 <- sign(df$"20.0%") == sign(df$"80.0%")
        df$null.50 <- sign(df$"25.0%") == sign(df$"75.0%")
        df$null.40 <- sign(df$"30.0%") == sign(df$"70.0%")
        df$null.30 <- sign(df$"35.0%") == sign(df$"65.0%")
        df$null.20 <- sign(df$"40.0%") == sign(df$"60.0%")
        df$null.10 <- sign(df$"45.0%") == sign(df$"55.0%")
        df
    }
)

if (model == "GT") {
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"
    outfiles <- list.files(
        sprintf("rds/GT/%s/", mtol),
        pattern = "ENSG*",
        full.names = TRUE
    )
    genes <- unique(gsub(".*(ENSG\\d+).*", "\\1", outfiles))
    infiles <- sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, genes)

    input_covars <- lapply(seq_along(infiles),
        function(i) {
            file <- infiles[[i]]
            inp <- readRDS(file)
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
            covars$snp <- names(inp)
            covars$gene <- genes[[i]]
            covars
        }
    )
    covar_df <- do.call(rbind, input_covars)
} else {

    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"

    # files <- list.files(dir)
    # files <- grep("refbias", files, value = TRUE)
    outfiles <- list.files(sprintf("rds/noGT/%s/", mtol), pattern = "ENSG*", full.names=TRUE)
    genes <- unique(gsub(".*(ENSG\\d+)..*", "\\1", outfiles))

    infiles <- list(
        normal_skin = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, genes),
        Psoriasis_skin = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, genes)
    )

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
    covar_dfs <- lapply(
        1:length(genes),
        function(i) {
            cat(i, "/", length(genes), "\n")
            infile_norm <- infiles[[1]][[i]]
            infile_pso <- infiles[[2]][[i]]
            inp_norm <- readRDS(infile_norm)
            inp_pso <- readRDS(infile_pso)
            covars_norm <- lapply(inp_norm, process_nogt)
            covars_pso <- lapply(inp_pso, process_nogt)
            covars_norm <- do.call(rbind, covars_norm)
            covars_pso <- do.call(rbind, covars_pso)
            covars_norm$condition <- "normal_skin"
            covars_pso$condition <- "Psoriasis_skin"
            covars <- rbind(covars_norm, covars_pso)
            covars$gene <- genes[[i]]
            covars$snp <- rownames(covars)
            covars
        }
    )
    covar_df <- do.call(rbind, covar_dfs)
}

names(dfs) <- methods
by <- if (model == "GT") {
    c("test", "gene", "snp")
} else {
    c("test", "gene", "snp", "condition")
}
# mdf <- merge(dfs[[1]], dfs[[2]], by = by, suffix = c(".vb", ".map"))
# mdf <- merge(mdf, dfs[[3]], by = by, suffix = c("", ".hmc"))
mdf <- merge(dfs[["vb"]], dfs[["sampling"]], by = by, suffix = c(".vb", ".hmc"))

mdf <- mdf %>%
    mutate(disc_vb = mean.vb - mean.hmc
    # , disc_map = mean.map - mean
)

# mtol <- sprintf("vb_%1.0e", tol)
mtol <- sprintf("vb")

cmdf <- merge(covar_df, mdf)
cmdf <- cmdf %>% mutate(discrepancy = mean.hmc - mean.vb)

x <- cmdf %>%
    arrange(-abs(disc_vb)) %>%
    top_n(50, disc_vb) %>%
    select(vb = mean.vb, hmc = mean.hmc, discrepancy = disc_vb, snp, gene)

saveRDS(x, sprintf("rds/%s_discrepancies_vb_%1.0e.rds", model, tol))

mname <- "ADVI"
method <- "vb"
r <- range(c(cmdf$mean.hmc, cmdf$mean.vb))
## not se mean but other hmc se
diag_vars <- c(
    "gene",
    "Rhat",
    "khat",
    "n_eff.hmc", "time.hmc", "n_ase",
    "se_mean.hmc", "sd.hmc",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)
for (x in diag_vars) {    
    g <- ggplot(cmdf) +
        aes_string(x, "abs(discrepancy)") +
        geom_point(size = 0.8) +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")
        ) +
        labs(x = x, y = "Discrepancy")

    ggsave(sprintf("%s/%s/diag/disc_%s_%s.png", fpath, model, x, method), width = 7, height = 7)
}

g <- ggplot(cmdf) +
    aes(time.hmc, disc_vb) +
    geom_point() +
    labs(x = "Time taken for HMC (s)", y = "Discrepancy (ADVI - HMC)")
ggsave(sprintf("%s/%s/diag/time_vb.png", fpath, model), width = 7, height = 7)

cmdf %>%
    group_by(gene) %>%
    summarise(time = sum(time.hmc), nsnps = n()) %>%
    ggplot() +
    aes(nsnps, time) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) -> g
ggsave(sprintf("%s/%s/time/nsnp.png", fpath, model), width = 7, height = 7)


cmdf <- cmdf[cmdf$n_eff.hmc > minEff & cmdf$Rhat < maxRhat, ]
## from PSIS paper, arxiv 1507.02646
# cmdf <- cmdf[cmdf$khat < 0.7, ]

gdf <- group_by(cmdf, gene)

lab_str <- "HMC: %s\nVB: %s\n"
cmdf$nullstr95 <- sprintf(lab_str,
    ifelse(cmdf$null.95.hmc, "no", "yes"),
    ifelse(cmdf$null.95.vb, "no", "yes")
)
cmdf$nullstr99 <- sprintf(lab_str,
    ifelse(cmdf$null.99.hmc, "no", "yes"),
    ifelse(cmdf$null.99.vb, "no", "yes")
)

null_levs <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "no",  "yes", "yes")
)
null_ord <- null_levs[c(1, 3, 2, 4)]
scale <- scale_colour_brewer(palette = "Paired", limits = null_levs, name = NULL)
cmdf$nullstr95 <- factor(cmdf$nullstr95, levels = null_ord)
cmdf$nullstr99 <- factor(cmdf$nullstr99, levels = null_ord)

stop()
mdfs <- cmdf
# mdfs <- cmdf[cmdf$mean.vb < 3, ]
lim <- range(c(mdfs$mean.hmc, mdfs$mean.vb))

g <- ggplot(mdfs[order(mdfs$nullstr95), ]) +
    aes(mean.hmc, mean.vb, colour = nullstr95) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = lim, y = lim) +
    labs(x = "MCMC estimate", y = "VB estimate") +
    scale
ggsave(sprintf("%s/%s/estimates/point-estimates-95.png", fpath, model), width = 7, height = 7)

g <- ggplot(mdfs[order(mdfs$nullstr99), ]) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = lim, y = lim) +
    labs(x = "MCMC estimate", y = "VB estimate") +
    scale
ggsave(sprintf("%s/%s/estimates/point-estimates-99.png", fpath, model), width = 7, height = 7)


## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t, ".hmc")
    levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
    sens_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            sens_vec(
                truth = factor(cmdf[[bc]], levels = c(TRUE, FALSE)),
                estimate = factor(cmdf[[column]], levels = c(TRUE, FALSE))
            )
        }
    )
    spec_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            spec_vec(
                truth = factor(cmdf[[bc]], levels = c(TRUE, FALSE)),
                estimate = factor(cmdf[[column]], levels = c(TRUE, FALSE))
            )
        }
    )
    time_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            sum(cmdf$time.hmc[cmdf[[column]]]) + sum(cmdf$time.vb)
        }
    )

    g <- ggplot() +
        aes(sens_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(cmdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
        ggtitle("Sensitivity vs total time for ADVI")
    ggsave(sprintf("%s/%s/time/time_vs_sens_vb_%s.png", fpath, model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(spec_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(cmdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
        ggtitle("Specificity vs total time for ADVI")
    ggsave(sprintf("%s/%s/time/time_vs_spec_vb_%s.png", fpath, model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(1 - spec_vb, sens_vb) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("VB") +
        lims(x = 0:1, y = 0:1)
    ggsave(sprintf("%s/%s/roc/roc_vb_%s.png", fpath, model, t), width = 5, height = 5)
}



gdf <- group_by(cmdf, gene)

## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t, ".hmc")
    levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
    sens_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            gdf %>%
                summarise(t = any(.data[[bc]]), e = any(.data[[column]])) %>%
                sens(
                    ## TRUE means that the HPD interval is one side of zero (sig)
                    ## FALSE means it overlaps zero (null)
                    truth = factor(t, levels = c(TRUE, FALSE)),
                    estimate = factor(e, levels = c(TRUE, FALSE))
                ) %>%
                pull(.estimate)
        }
    )
    spec_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            gdf %>%
                summarise(t = any(.data[[bc]]), e = any(.data[[column]])) %>%
                spec(
                    ## TRUE means that the HPD interval is one side of zero (sig)
                    ## FALSE means it overlaps zero (null)
                    truth = factor(t, levels = c(TRUE, FALSE)),
                    estimate = factor(e, levels = c(TRUE, FALSE))
                ) %>%
                pull(.estimate)
        }
    )
    time_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            gd <- gdf %>%
                summarise(e = any(.data[[column]]))
            sum(cmdf$time.hmc[cmdf$gene %in% gd$gene[gd$e]]) + sum(cmdf$time.vb)
        }
    )

    g <- ggplot() +
        aes(sens_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(cmdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
        ggtitle("Sensitivity vs total time for ADVI")
    ggsave(sprintf("%s/%s/time/time_vs_sens_gene_vb_%s.png", fpath, model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(spec_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(cmdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(cmdf$time.hmc), time_vb)) +
        ggtitle("Specificity vs total time for ADVI")
    ggsave(sprintf("%s/%s/time/time_vs_spec_gene_vb_%s.png", fpath, model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(1 - spec_vb, sens_vb) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("VB") +
        lims(x = 0:1, y = 0:1)
    ggsave(sprintf("%s/%s/roc/roc_vb_gene_%s.png", fpath, model, t), width = 5, height = 5)
}

g <- ggplot(cmdf) +
    aes(time.hmc, time.vb) +
    geom_pointdensity() +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_viridis() +
    theme(legend.position = "bottom") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Time (s) for HMC", y = "Time (s) for ADVI")
ggsave(sprintf("%s/%s/time/time_hmc_vs_vb.png", fpath, model), width = 5, height = 5)

cn <- lapply(dfs, colnames)
cn <- Reduce(intersect, cn)
dfs_int <- lapply(dfs, function(x) x[, cn])
df_int <- do.call(rbind, dfs_int)

g <- ggplot(df_int) +
    aes(x = time.hmc, y = method) +
    scale_x_log10(name = "Time (s)") +
    stat_pointinterval()
# ggsave("tmp.png")

g <- ggplot(df_int) +
    aes(x = time, colour = method) +
    geom_density() +
    scale_x_log10(name = "Time (s)") +
    scale_colour_brewer(palette = "Set1", name = "Method") +
    ylab("Density") +
    theme(legend.position = "bottom")
ggsave(sprintf("%s/%s/time/time_comparison.png", fpath, model), width = 5, height = 5)
