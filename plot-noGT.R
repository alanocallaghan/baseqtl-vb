library("argparse")
library("ggplot2")
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
    method <- "optimizing"
}
mname <- switch(method,
    "vb" = "ADVI",
    "optimizing" = "MAP",
    "sampling" = "HMC"
)

dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"

files <- list.files(dir)
files <- grep("refbias", files, value = TRUE)
genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", files))


all_stan_res <- lapply(
    list.files(sprintf("rds/noGT/%s/", method), pattern = "ENSG*", full.names=TRUE),
    readRDS
)
all_stan_res <- lapply(all_stan_res,
    function(x) {
        if (!length(x)) return(NULL)
        c(x[[1]], x[[2]])
    }
)
all_stan_res <- Reduce(c, all_stan_res)
all_stan_res <- do.call(rbind, all_stan_res)
approx_res_df <- as.data.frame(all_stan_res)
approx_res_df$test <- paste(approx_res_df$gene, approx_res_df$snp , sep = "_")

summary_files <- list.files(dir, pattern = "stan.summary", full.names = TRUE)
sdata <- data.frame(
    tissue = gsub(".*(Psoriasis_skin|normal_skin).*", "\\1", summary_files)
)
ll <- lapply(seq_along(summary_files),
    function(i) {
        x <- read.table(summary_files[[i]], header = TRUE)
        x$tissue <- sdata[i, ]
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

g <- ggplot() +
    aes(approx_res_df$mean, sample_res_df$log2_aFC_mean) +
    geom_pointdensity() +
    scale_colour_viridis() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/noGT/%s_mcmc_estimates_all.png", method), width = 7, height = 7)


mdf <- mdf[order(mdf$Null.99), ]

g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.99) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/noGT/%s_mcmc_estimates_all_categorical_99.png", method), width = 7, height = 7)



mdf <- mdf[order(mdf$Null.95), ]
g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.95) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/noGT/%s_mcmc_estimates_all_categorical_95.png", method), width = 7, height = 7)

mdf <- mdf[order(mdf$Null.95.50), ]
g <- ggplot(mdf) +
    aes(mean, log2_aFC_mean, colour = Null.95.50) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = sprintf("%s estimate", mname), y = "MCMC estimate")
ggsave(file = sprintf("fig/noGT/%s_mcmc_estimates_all_categorical_50.png", method), width = 7, height = 7)


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
ggsave(sprintf("fig/noGT/time_dist_all_%s.png", method), width = 6, height = 6)



# sample_data <- readRDS(
#     "/rds/project/cew54/rds-cew54-wallace-share/Projects/baseqtl/data/btrecase.GT.nGT.all.rds"
# )
# ss <- as.data.frame(sample_data$RNA)
# ss$test <- paste(ss$Gene_id, ss$tag , sep = "_")
# ss$null.50 <- sign(ss$"log2_aFC_75%") == sign(ss$"log2_aFC_25%")

# ss <- ss[ss$test %in% approx_res_df$test, ]
# im <- match(ss$test, approx_res_df$test)
# approx_res_df_sub <- approx_res_df[im, ]
# approx_res_df_sub$null.50 <- sign(approx_res_df_sub$"75.0%") == sign(approx_res_df_sub$"25.0%")
# stopifnot(all(ss$test == approx_res_df_sub$test))
# r2 <- range(c(ss$log2_aFC_mean, approx_res_df_sub$mean))
# ## todo


# g <- ggplot() +
#     aes(approx_res_df_sub$mean, ss$log2_aFC_mean) +
#     geom_pointdensity() +
#     scale_colour_viridis() +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     lims(x = r2, y = r2) +
#     labs(x = sprintf("%s estimate", method), y = "MCMC estimate")
# ggsave(file = sprintf("fig/noGT/%s_mcmc_estimates.png", method), width = 7, height = 7)

# mdf2 <- merge(approx_res_df_sub, ss, by = "test", suffix = c(".VB", ".HMC"))
# mdf2$null.99.VB <- ifelse(mdf2$null.99.VB, "no", "yes")
# mdf2$null.95.VB <- ifelse(mdf2$null.95.VB, "no", "yes")
# mdf2$null.50.VB <- ifelse(mdf2$null.50.VB, "no", "yes")
# mdf2$null.50.HMC <- ifelse(mdf2$null.50.HMC, "no", "yes")


# mdf2$Null.99 <- sprintf(lab_str, mdf2$null.99.VB, mdf2$null.99.HMC)
# mdf2$Null.95 <- sprintf(lab_str, mdf2$null.95.VB, mdf2$null.95.HMC)
# mdf2$Null.50 <- sprintf(lab_str, mdf2$null.50.VB, mdf2$null.50.HMC)
# mdf2$Null.95.50 <- sprintf(lab_str, mdf2$null.50.VB, mdf2$null.95.HMC)

# mdf2$Null.99 <- factor(mdf2$Null.99, levels = flev)
# mdf2$Null.95 <- factor(mdf2$Null.95, levels = flev)
# mdf2$Null.50 <- factor(mdf2$Null.50, levels = flev)
# mdf2$Null.95.50 <- factor(mdf2$Null.95.50, levels = flev)


# mdf2 <- mdf2[order(mdf2$Null.99), ]
# g <- ggplot(mdf2) +
#     aes(mean, log2_aFC_mean, colour = Null.99) +
#     geom_point() +
#     scale_colour_brewer(palette = "Paired", limits = levs) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     lims(x = r2, y = r2) +
#     labs(x = sprintf("%s estimate", method), y = "MCMC estimate")

# # ggMarginal(g)
# ggsave(file = sprintf("fig/noGT/%s_mcmc_estimates_categorical_99.png", method), width = 7, height = 7)

# mdf2 <- mdf2[order(mdf2$Null.95), ]
# g <- ggplot(mdf2) +
#     aes(mean, log2_aFC_mean, colour = Null.95) +
#     geom_point() +
#     scale_colour_brewer(palette = "Paired", limits = levs) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     lims(x = r2, y = r2) +
#     labs(x = sprintf("%s estimate", method), y = "MCMC estimate")

# # ggMarginal(g)
# ggsave(file = sprintf("fig/noGT/%s_mcmc_estimates_categorical_95.png", method), width = 7, height = 7)


# mdf2 <- mdf2[order(mdf2$Null.95.50), ]
# g <- ggplot(mdf2) +
#     aes(mean, log2_aFC_mean, colour = Null.95.50) +
#     geom_point() +
#     scale_colour_brewer(palette = "Paired", limits = levs) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     lims(x = r2, y = r2) +
#     labs(x = sprintf("%s estimate", method), y = "MCMC estimate")

# # ggMarginal(g)
# ggsave(file = sprintf("fig/noGT/%s_mcmc_estimates_categorical_50.png", method), width = 7, height = 7)


# summ2 <- mdf2 %>%
#     group_by(gene) %>%
#     summarise(
#         anynull99vb = any(null.99.VB == "no"),
#         anynull99hmc = any(null.99.HMC == "no"),
#         anynull95vb = any(null.95.VB == "no"),
#         anynull95hmc = any(null.95.HMC == "no"),
#         anynull50vb = any(null.50.VB == "no"),
#         anynull50hmc = any(null.50.HMC == "no")
#     )
# summ2[-1] <- lapply(summ2[-1], factor)

# stop()

# ggplot() +
#     aes(mdf2$time) +
#     geom_histogram(
#         bins = nclass.FD(mdf2$time),
#         colour = "grey60",
#         fill = "grey90",
#         boundary = 0
#     ) +
#     geom_vline(xintercept = median(mdf2$time), linetype = "dashed") +
#     labs(x = "Time (s)", y = "Frequency")
# ggsave(sprintf("fig/noGT/time_dist_all_%s.png", method), width = 6, height = 6)
