library("ggplot2")
library("ggdist")
library("ggpointdensity")
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
args <- parser$parse_args()

if (is.null(model <- args[["model"]])) {
    model <- "noGT"
}

dfs <- lapply(c("vb", "optimizing", "sampling"),
    function(m) {
        print(m)
        f <- list.files(path = file.path("rds", model, m), pattern = "ENSG.*.rds", full.names = TRUE)
        input <- lapply(f, readRDS)
        if (model == "noGT") {
            input <- lapply(input, 
                function(x) {
                    if (!length(x)) return(NULL)
                    c(x[[1]], x[[2]])
                }
            )
            input <- Reduce(c, input)
            df <- do.call(rbind, input)
            df <- as.data.frame(df)
            df$test <- paste(df$gene, df$snp , sep = "_")
        } else {
            df <- do.call(rbind, lapply(input, function(x) do.call(rbind, x)))
            df <- as.data.frame(df)
            df$test <- paste(df$gene, rownames(df) , sep = "_")
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
cn <- lapply(dfs, colnames)
cn <- Reduce(intersect, cn)
dfs <- lapply(dfs, function(x) x[, cn])

mdf <- merge(dfs[[1]], dfs[[2]], by = "test", suffix = c(".vb", ".map"))
mdf <- merge(mdf, dfs[[3]], by = "test", suffix = c("", ".hmc"))

# quantile(mdf$time, probs = seq(0, 1, length.out = 20))

mdf %>%
    group_by(gene) %>%
    summarise(time = sum(time), nsnps = n()) %>%
    ggplot() +
    aes(nsnps, time) +
    geom_point() +
    geom_smooth(method = "lm") -> g
ggsave(sprintf("fig/%s/time_nsnp.png", model))


# sum(mdf$time)
# sum(mdf$time[mdf$null.99.vb]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.95.vb]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.90.vb]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.80.vb]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.70.vb]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.60.vb]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.50.vb]) + sum(mdf$time.vb)


# sum(mdf$time[mdf$null.99.map]) + sum(mdf$time.map)
# sum(mdf$time[mdf$null.95.map]) + sum(mdf$time.map)
# sum(mdf$time[mdf$null.90.map]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.80.map]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.70.map]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.60.map]) + sum(mdf$time.vb)
# sum(mdf$time[mdf$null.50.map]) + sum(mdf$time.map)


gdf <- group_by(mdf, gene)

## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t)
    levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
    
    sens_map <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.map", x)
            sens_vec(
                ## TRUE means that the HPD interval is one side of zero (sig)
                ## FALSE means it overlaps zero (null)
                truth = factor(mdf[[bc]], levels = c(TRUE, FALSE)),
                estimate = factor(mdf[[column]], levels = c(TRUE, FALSE))
            )
        }
    )
    spec_map <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.map", x)
            spec_vec(
                truth = factor(mdf[[bc]], levels = c(TRUE, FALSE)),
                estimate = factor(mdf[[column]], levels = c(TRUE, FALSE))
            )
        }
    )

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
            sum(mdf$time[mdf[[column]]]) + sum(mdf$time.vb)
        }
    )
    time_map <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.map", x)
            sum(mdf$time[mdf[[column]]]) + sum(mdf$time.map)
        }
    )
    g <- ggplot() +
        aes(sens_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(mdf$time), time_vb)) +
        ggtitle("Sensitivity vs total time for ADVI")
    ggsave(sprintf("fig/%s/time_vs_sens_vb_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(sens_map, time_map) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time),
            label = "Total time without screening",
            vjust = -0.2,
            linetype = "dashed"
        ) +
        ylim(0, max(sum(mdf$time), time_map)) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ggtitle("Sensitivity vs total time for MAP")
    ggsave(sprintf("fig/%s/time_vs_sens_map_%s.png", model, t), width = 5, height = 5)
    g <- ggplot() +
        aes(spec_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(mdf$time), time_vb)) +
        ggtitle("Specificity vs total time for ADVI")

    ggsave(sprintf("fig/%s/time_vs_spec_vb_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(spec_map, time_map) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time),
            label = "Total time without screening",
            vjust = -0.2,
            linetype = "dashed"
        ) +
        ylim(0, max(sum(mdf$time), time_map)) +
        labs(x = "Specificity", y = "Total time (s)") +
        ggtitle("Specificity vs total time for MAP")
    ggsave(sprintf("fig/%s/time_vs_spec_map_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(1 - spec_map, sens_map) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("MAP") +
        lims(x = 0:1, y = 0:1)
    ggsave(sprintf("fig/%s/roc_map_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(1 - spec_vb, sens_vb) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("VB") +
        lims(x = 0:1, y = 0:1)
    ggsave(sprintf("fig/%s/roc_vb_%s.png", model, t), width = 5, height = 5)
}



gdf <- group_by(mdf, gene)

## by snp
for (t in c(99, 95)) {
    bc <- paste0("null.", t)
    levs <- c(99, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10)
    
    sens_map <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.map", x)
            gd <- gdf %>%
                summarise(t = any(.data[[bc]]), e = any(.data[[column]]))
            
            mdf %>% mutate(
                    t = gene %in% gd[["gene"]][gd$t],
                    e = gene %in% gd[["gene"]][gd$e]
                ) %>%
                sens(
                    ## TRUE means that the HPD interval is one side of zero (sig)
                    ## FALSE means it overlaps zero (null)
                    truth = factor(t, levels = c(TRUE, FALSE)),
                    estimate = factor(e, levels = c(TRUE, FALSE))
                ) %>%
                pull(.estimate)
        }
    )

    spec_map <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.map", x)
            gd <- gdf %>%
                summarise(t = any(.data[[bc]]), e = any(.data[[column]]))
            mdf %>% mutate(
                    t = gene %in% gd[["gene"]][gd$t],
                    e = gene %in% gd[["gene"]][gd$e]
                ) %>%
                spec(
                    ## TRUE means that the HPD interval is one side of zero (sig)
                    ## FALSE means it overlaps zero (null)
                    truth = factor(t, levels = c(TRUE, FALSE)),
                    estimate = factor(e, levels = c(TRUE, FALSE))
                ) %>%
                pull(.estimate)
        }
    )
    
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
            sum(mdf$time[mdf$gene %in% gd$gene[gd$e]]) + sum(mdf$time.vb)
        }
    )
    time_map <- sapply(
        levs,
        function(x) {
            column <- sprintf("null.%s.map", x)
            gd <- gdf %>%
                summarise(e = any(.data[[column]]))
            sum(mdf$time[mdf$gene %in% gd$gene[gd$e]]) + sum(mdf$time.map)
        }
    )
    g <- ggplot() +
        aes(sens_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ylim(0, max(sum(mdf$time), time_vb)) +
        ggtitle("Sensitivity vs total time for ADVI")
    ggsave(sprintf("fig/%s/time_vs_sens_gene_vb_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(sens_map, time_map) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time),
            label = "Total time without screening",
            vjust = -0.2,
            linetype = "dashed"
        ) +
        ylim(0, max(sum(mdf$time), time_map)) +
        labs(x = "Sensitivity", y = "Total time (s)") +
        ggtitle("Sensitivity vs total time for MAP")
    ggsave(sprintf("fig/%s/time_vs_sens_gene_map_%s.png", model, t), width = 5, height = 5)
    g <- ggplot() +
        aes(spec_vb, time_vb) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time),
            label = "Total time without screening", 
            vjust = -0.2,
            linetype = "dashed"
        ) +
        labs(x = "Specificity", y = "Total time (s)") +
        ylim(0, max(sum(mdf$time), time_vb)) +
        ggtitle("Specificity vs total time for ADVI")

    ggsave(sprintf("fig/%s/time_vs_spec_gene_vb_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(spec_map, time_map) +
        geom_line() +
        geom_texthline(
            yintercept = sum(mdf$time),
            label = "Total time without screening",
            vjust = -0.2,
            linetype = "dashed"
        ) +
        ylim(0, max(sum(mdf$time), time_map)) +
        labs(x = "Specificity", y = "Total time (s)") +
        ggtitle("Specificity vs total time for MAP")
    ggsave(sprintf("fig/%s/time_vs_spec_gene_map_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(1 - spec_map, sens_map) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("MAP") +
        lims(x = 0:1, y = 0:1)
    ggsave(sprintf("fig/%s/roc_map_gene_%s.png", model, t), width = 5, height = 5)

    g <- ggplot() +
        aes(1 - spec_vb, sens_vb) +
        geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(x = "1 - specificity", y = "Sensitivity") +
        ggtitle("VB") +
        lims(x = 0:1, y = 0:1)
    ggsave(sprintf("fig/%s/roc_vb_gene_%s.png", model, t), width = 5, height = 5)
}



g <- ggplot(mdf) +
    aes(time.vb, time.map) +
    geom_pointdensity() +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_viridis() +
    theme(legend.position = "bottom") +
    # geom_abline(slope = 1, intercept = 0) +
    labs(x = "Time (s) for ADVI", y = "Time (s) for MAP")
ggsave(sprintf("fig/%s/time_vb_vs_map.png", model), width = 5, height = 5)

g <- ggplot(mdf) +
    aes(time, time.vb) +
    geom_pointdensity() +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_viridis() +
    theme(legend.position = "bottom") +
    # geom_abline(slope = 1, intercept = 0) +
    labs(x = "Time (s) for HMC", y = "Time (s) for ADVI")
ggsave(sprintf("fig/%s/time_hmc_vs_vb.png", model), width = 5, height = 5)

g <- ggplot(mdf) +
    aes(time, time.map) +
    geom_pointdensity() +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_viridis() +
    theme(legend.position = "bottom") +
    # geom_abline(slope = 1, intercept = 0) +
    labs(x = "Time (s) for HMC", y = "Time (s) for MAP")
ggsave(sprintf("fig/%s/time_hmc_vs_map.png", model), width = 5, height = 5)





df <- do.call(rbind, dfs)

g <- ggplot(df) +
    aes(x = time, y = method) +
    scale_x_log10(name = "Time (s)") +
    stat_pointinterval()

# ggsave("tmp.png")


g <- ggplot(df) +
    aes(x = time, colour = method) +
    geom_density() +
    scale_x_log10(name = "Time (s)") +
    scale_colour_brewer(palette = "Set1", name = "Method") +
    ylab("Density") +
    theme(legend.position = "bottom")
ggsave(sprintf("fig/%s/time-comparison.png", model), width = 5, height = 5)

