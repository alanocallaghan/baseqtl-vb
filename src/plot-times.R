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
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-3,
    type = "double"
)


args <- parser$parse_args()
tol <- args[["tolerance"]]

if (is.null(model <- args[["model"]])) {
    model <- "GT"
    model <- "noGT"
}

# "optimizing",
dfs <- lapply(c("vb", "sampling"),
    function(m) {
        cat(m, "\n")
        mtol <- if (m == "vb") sprintf("%s_%1.0e", m, tol) else m
        path <- file.path("rds", model, mtol)
        f <- list.files(path = path, pattern = "ENSG.*.rds", full.names = TRUE)
        input <- lapply(f, readRDS)
        if (model == "noGT") {
            # input <- lapply(input, 
            #     function(x) {
            #         if (!length(x)) return(NULL)
            #         c(x[[1]], x[[2]])
            #     }
            # )
            # input <- Reduce(c, input)
            df <- do.call(rbind, input)
            df <- as.data.frame(df)
            df$test <- paste(df$gene, df$snp , sep = "_")
        } else {
            # dfs <- lapply(input, function(x) {
                # d <- do.call(rbind, x)
                # d$snp <- names(x)
                # d
            # })
            df <- do.call(rbind, input)
            df <- as.data.frame(df)
            df$test <- paste(df$gene, df$snp , sep = "_")
        }
        if (m == "sampling") {
            df <- df[df$n_eff > 1000 & df$Rhat < 1.05, ]
        } else if (m == "vb") {
            ## from PSIS paper, arxiv 1507.02646
            df <- df[df$khat < 0.7, ]
        }
        
        df$method <- gsub("old/", "", m)
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


by <- if (model == "GT") c("test", "gene", "snp") else c("test", "gene", "snp", "condition")
# mdf <- merge(dfs[[1]], dfs[[2]], by = by, suffix = c(".vb", ".map"))
# mdf <- merge(mdf, dfs[[3]], by = by, suffix = c("", ".hmc"))
mdf <- merge(dfs[[1]], dfs[[2]], by = by, suffix = c(".vb", ".hmc"))

mdf <- mdf %>%
    mutate(disc_vb = mean.vb - mean.hmc
    # , disc_map = mean.map - mean
)


if (model == "GT") {
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"
    outfiles <- list.files("rds/GT/vb/", pattern = "ENSG*", full.names = TRUE)
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
    outfiles <- list.files("rds/noGT/vb/", pattern = "ENSG*", full.names=TRUE)
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
mdf <- merge(covar_df, mdf)

mdf <- mdf %>% mutate(discrepancy = mean.hmc - mean.vb)


x <- mdf %>%
    arrange(-abs(disc_vb)) %>%
    top_n(50, disc_vb) %>%
    select(vb = mean.vb, hmc = mean.hmc, discrepancy = disc_vb, snp, gene)

saveRDS(x, sprintf("rds/%s_discrepancies_vb_%1.0e.rds", model, tol))

sprintf("%s_%1.0e", m, tol)


mname <- "ADVI"
method <- "vb"
r <- range(c(mdf$mean.hmc, mdf$mean.vb))
## not se mean but other hmc se
diag_vars <- c(
    "gene",
    # "Rhat.hmc",
    "khat",
    "n_eff.hmc", "time.hmc", "n_ase",
    "se_mean.hmc", "sd.hmc",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)
for (x in diag_vars) {
    scale <- if (x == "gene") scale_colour_discrete(guide="none") else scale_colour_viridis()
    
    g <- ggplot(mdf[order(mdf[[x]]), ]) +
        aes_string("mean.hmc", "mean.vb", colour = x) +
        geom_point(size = 0.8) +
        scale +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        lims(x = r, y = r) +
        labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
    ggsave(sprintf("fig_%1.0e/%s/diag/%s_%s.png", tol, model, x, method), width = 7, height = 7)
    
    g <- ggplot(mdf) +
        aes_string(x, "abs(discrepancy)") +
        geom_point(size = 0.8) +
        # scale_colour_viridis() +
        # geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        # lims(x = r, y = r) +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")
        ) +
        labs(x = x, y = "Discrepancy")

    ggsave(sprintf("fig_%1.0e/%s/diag/disc_%s_%s.png", tol, model, x, method), width = 7, height = 7)
}

g <- ggplot(mdf) +
    aes(time.hmc, disc_vb) +
    geom_point() +
    labs(x = "Time taken for HMC (s)", y = "Discrepancy (ADVI - HMC)")
ggsave(sprintf("fig_%1.0e/%s/diag/time_vb.png", tol, model), width = 7, height = 7)

# g <- ggplot(mdf) +
#     aes(time, disc_map) +
#     geom_point() +
#     labs(x = "Time taken for HMC (s)", y = "Discrepancy (MAP - HMC)")
# ggsave(sprintf("fig/%s/diag/time_map.png", model))


# quantile(mdf$time, probs = seq(0, 1, length.out = 20))

mdf %>%
    group_by(gene) %>%
    summarise(time = sum(time.hmc), nsnps = n()) %>%
    ggplot() +
    aes(nsnps, time) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) -> g
ggsave(sprintf("fig_%1.0e/%s/time/nsnp.png", tol, model), width = 7, height = 7)

gdf <- group_by(mdf, gene)


lab_str <- "HMC: %s\nVB: %s\n"
mdf$nullstr95 <- sprintf(lab_str,
    ifelse(mdf$null.95.hmc, "no", "yes"),
    ifelse(mdf$null.95.vb, "no", "yes")
)
mdf$nullstr99 <- sprintf(lab_str,
    ifelse(mdf$null.99.hmc, "no", "yes"),
    ifelse(mdf$null.99.vb, "no", "yes")
)

null_levs <- sprintf(
    lab_str,
    c("no", "yes", "yes", "no"),
    c("no", "no",  "yes", "yes")
)
null_ord <- null_levs[c(1, 3, 2, 4)]
scale <- scale_colour_brewer(palette = "Paired", limits = null_levs, name = NULL)
mdf$nullstr95 <- factor(mdf$nullstr95, levels = null_ord)
mdf$nullstr99 <- factor(mdf$nullstr99, levels = null_ord)

g <- ggplot(mdf[order(mdf$nullstr95), ]) +
    aes(mean.hmc, mean.vb, colour = nullstr95) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "MCMC estimate", y = "VB estimate") +
    scale
ggsave(sprintf("fig_%1.0e/%s/estimates/pt-95.png", tol, model), width = 7, height = 7)
g <- ggplot(mdf[order(mdf$nullstr99), ]) +
    aes(mean.hmc, mean.vb, colour = nullstr99) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "MCMC estimate", y = "VB estimate") +
    scale
ggsave(sprintf("fig_%1.0e/%s/estimates/pt-99.png", tol, model), width = 7, height = 7)

# stop()

## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t, ".hmc")
    levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
    
    # sens_map <- sapply(
    #     levs,
    #     function(x) {
    #         column <- sprintf("null.%s.map", x)
    #         sens_vec(
    #             ## TRUE means that the HPD interval is one side of zero (sig)
    #             ## FALSE means it overlaps zero (null)
    #             truth = factor(mdf[[bc]], levels = c(TRUE, FALSE)),
    #             estimate = factor(mdf[[column]], levels = c(TRUE, FALSE))
    #         )
    #     }
    # )
    # spec_map <- sapply(
    #     levs,
    #     function(x) {
    #         column <- sprintf("null.%s.map", x)
    #         spec_vec(
    #             truth = factor(mdf[[bc]], levels = c(TRUE, FALSE)),
    #             estimate = factor(mdf[[column]], levels = c(TRUE, FALSE))
    #         )
    #     }
    # )

    sens_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            sens_vec(
                truth = factor(mdf[[bc]], levels = c(TRUE, FALSE)),
                estimate = factor(mdf[[column]], levels = c(TRUE, FALSE))
            )
        }
    )
    spec_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            spec_vec(
                truth = factor(mdf[[bc]], levels = c(TRUE, FALSE)),
                estimate = factor(mdf[[column]], levels = c(TRUE, FALSE))
            )
        }
    )

    time_vb <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.vb", x)
            sum(mdf$time.hmc[mdf[[column]]]) + sum(mdf$time.vb)
        }
    )
    # time_map <- sapply(
    #     levs,
    #     function(x) {
    #         column <- sprintf("null.%s.map", x)
    #         sum(mdf$time.hmc[mdf[[column]]]) + sum(mdf$time.map)
    #     }
    # )
    g <- ggplot() +
        aes(sens_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(mdf$time.hmc), time_vb)) +
        ggtitle("Sensitivity vs total time for ADVI")
    ggsave(sprintf("fig_%1.0e/%s/time/time_vs_sens_vb_%s.png", tol, model, t), width = 5, height = 5)

    # g <- ggplot() +
    #     aes(sens_map, time_map) +
    #     geom_line() +
    #     geom_texthline(
    #         yintercept = sum(mdf$time.hmc),
    #         label = "Total time without screening",
    #         vjust = -0.2,
    #         linetype = "dashed"
    #     ) +
    #     ylim(0, max(sum(mdf$time.hmc), time_map)) +
    #     labs(x = "Sensitivity", y = "Total time (s)") +
    #     ggtitle("Sensitivity vs total time for MAP")
    # ggsave(sprintf("fig/%s/time/time_vs_sens_map_%s.png", model, t), width = 5, height = 5)
    g <- ggplot() +
        aes(spec_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(mdf$time.hmc), time_vb)) +
        ggtitle("Specificity vs total time for ADVI")

    ggsave(sprintf("fig_%1.0e/%s/time/time_vs_spec_vb_%s.png", tol, model, t), width = 5, height = 5)

    # g <- ggplot() +
    #     aes(spec_map, time_map) +
    #     geom_line() +
    #     geom_texthline(
    #         yintercept = sum(mdf$time.hmc),
    #         label = "Total time without screening",
    #         vjust = -0.2,
    #         linetype = "dashed"
    #     ) +
    #     ylim(0, max(sum(mdf$time.hmc), time_map)) +
    #     labs(x = "Specificity", y = "Total time (s)") +
    #     ggtitle("Specificity vs total time for MAP")
    # ggsave(sprintf("fig/%s/time/time_vs_spec_map_%s.png", model, t), width = 5, height = 5)

    # g <- ggplot() +
    #     aes(1 - spec_map, sens_map) +
    #     geom_line() +
    #     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    #     labs(x = "1 - specificity", y = "Sensitivity") +
    #     ggtitle("MAP") +
    #     lims(x = 0:1, y = 0:1)
    # ggsave(sprintf("fig/%s/roc/roc_map_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(1 - spec_vb, sens_vb) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("VB") +
        lims(x = 0:1, y = 0:1)
    ggsave(sprintf("fig_%1.0e/%s/roc/roc_vb_%s.png", tol, model, t), width = 5, height = 5)
}



gdf <- group_by(mdf, gene)

## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t, ".hmc")
    levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
    
    # sens_map <- sapply(
    #     levs,
    #     function(x) {
    #         column <- sprintf("null.%s.map", x)
    #         gd <- gdf %>%
    #             summarise(t = any(.data[[bc]]), e = any(.data[[column]]))
            
    #         mdf %>% mutate(
    #                 t = gene %in% gd[["gene"]][gd$t],
    #                 e = gene %in% gd[["gene"]][gd$e]
    #             ) %>%
    #             sens(
    #                 ## TRUE means that the HPD interval is one side of zero (sig)
    #                 ## FALSE means it overlaps zero (null)
    #                 truth = factor(t, levels = c(TRUE, FALSE)),
    #                 estimate = factor(e, levels = c(TRUE, FALSE))
    #             ) %>%
    #             pull(.estimate)
    #     }
    # )

    # spec_map <- sapply(
    #     levs,
    #     function(x) {
    #         column <- sprintf("null.%s.map", x)
    #         gd <- gdf %>%
    #             summarise(t = any(.data[[bc]]), e = any(.data[[column]]))
    #         mdf %>% mutate(
    #                 t = gene %in% gd[["gene"]][gd$t],
    #                 e = gene %in% gd[["gene"]][gd$e]
    #             ) %>%
    #             spec(
    #                 ## TRUE means that the HPD interval is one side of zero (sig)
    #                 ## FALSE means it overlaps zero (null)
    #                 truth = factor(t, levels = c(TRUE, FALSE)),
    #                 estimate = factor(e, levels = c(TRUE, FALSE))
    #             ) %>%
    #             pull(.estimate)
    #     }
    # )
    
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
            sum(mdf$time.hmc[mdf$gene %in% gd$gene[gd$e]]) + sum(mdf$time.vb)
        }
    )
    # time_map <- sapply(
    #     levs,
    #     function(x) {
    #         column <- sprintf("null.%s.map", x)
    #         gd <- gdf %>%
    #             summarise(e = any(.data[[column]]))
    #         sum(mdf$time.hmc[mdf$gene %in% gd$gene[gd$e]]) + sum(mdf$time.map)
    #     }
    # )
    g <- ggplot() +
        aes(sens_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(mdf$time.hmc), time_vb)) +
        ggtitle("Sensitivity vs total time for ADVI")
    ggsave(sprintf("fig_%1.0e/%s/time/time_vs_sens_gene_vb_%s.png", tol, model, t), width = 5, height = 5)

    # g <- ggplot() +
    #     aes(sens_map, time_map) +
    #     geom_line() +
    #     geom_texthline(
    #         yintercept = sum(mdf$time.hmc),
    #         label = "Total time without screening",
    #         vjust = -0.2,
    #         linetype = "dashed"
    #     ) +
    #     ylim(0, max(sum(mdf$time.hmc), time_map)) +
    #     labs(x = "Sensitivity", y = "Total time (s)") +
    #     ggtitle("Sensitivity vs total time for MAP")
    # ggsave(sprintf("fig/%s/time/time_vs_sens_gene_map_%s.png", model, t), width = 5, height = 5)
    g <- ggplot() +
        aes(spec_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time.hmc),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(mdf$time.hmc), time_vb)) +
        ggtitle("Specificity vs total time for ADVI")

    ggsave(sprintf("fig_%1.0e/%s/time/time_vs_spec_gene_vb_%s.png", tol, model, t), width = 5, height = 5)

    # g <- ggplot() +
    #     aes(spec_map, time_map) +
    #     geom_line() +
    #     geom_texthline(
    #         yintercept = sum(mdf$time.hmc),
    #         label = "Total time without screening",
    #         vjust = -0.2,
    #         linetype = "dashed"
    #     ) +
    #     ylim(0, max(sum(mdf$time.hmc), time_map)) +
    #     labs(x = "Specificity", y = "Total time (s)") +
    #     ggtitle("Specificity vs total time for MAP")
    # ggsave(sprintf("fig/%s/time/time_vs_spec_gene_map_%s.png", model, t), width = 5, height = 5)

    # g <- ggplot() +
    #     aes(1 - spec_map, sens_map) +
    #     geom_line() +
    #     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    #     labs(x = "1 - specificity", y = "Sensitivity") +
    #     ggtitle("MAP") +
    #     lims(x = 0:1, y = 0:1)
    # ggsave(sprintf("fig/%s/roc/roc_map_gene_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(1 - spec_vb, sens_vb) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("VB") +
        lims(x = 0:1, y = 0:1)
    ggsave(sprintf("fig_%1.0e/%s/roc/roc_vb_gene_%s.png", tol, model, t), width = 5, height = 5)
}



# g <- ggplot(mdf) +
#     aes(time.vb, time.map) +
#     geom_pointdensity() +
#     scale_x_log10() +
#     scale_y_log10() +
#     scale_colour_viridis() +
#     theme(legend.position = "bottom") +
#     # geom_abline(slope = 1, intercept = 0) +
#     labs(x = "Time (s) for ADVI", y = "Time (s) for MAP")
# ggsave(sprintf("fig/%s/time/time_vb_vs_map.png", model), width = 5, height = 5)

g <- ggplot(mdf) +
    aes(time.hmc, time.vb) +
    geom_pointdensity() +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_viridis() +
    theme(legend.position = "bottom") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Time (s) for HMC", y = "Time (s) for ADVI")
ggsave(sprintf("fig_%1.0e/%s/time/time_hmc_vs_vb.png", tol, model), width = 5, height = 5)

# g <- ggplot(mdf) +
#     aes(time.hmc, time.map) +
#     geom_pointdensity() +
#     scale_x_log10() +
#     scale_y_log10() +
#     scale_colour_viridis() +
#     theme(legend.position = "bottom") +
#     # geom_abline(slope = 1, intercept = 0) +
#     labs(x = "Time (s) for HMC", y = "Time (s) for MAP")
# ggsave(sprintf("fig/%s/time/time_hmc_vs_map.png", model), width = 5, height = 5)



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
ggsave(sprintf("fig_%1.0e/%s/time/time_comparison.png", tol, model), width = 5, height = 5)

