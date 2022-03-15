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
parser$add_argument(
    "-n", "--nogt",
    action = "store_true",
    default = TRUE
)


args <- parser$parse_args()
method <- args[["method"]]
tol <- args[["tolerance"]]
model <- if (args[["nogt"]]) "noGT" else "GT"

mname <- switch(method,
    "vb" = "ADVI",
    "optimizing" = "MAP",
    "sampling" = "HMC"
)

mtol <- if (method == "vb") sprintf("vb_%1.0e", tol) else method

maxRhat <- 1.1
minEff <- 500


if (model == "GT") {
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
            cn <- setdiff(
                colnames(out),
                c("n_eff", "Rhat", "null.99", "gene", "time", "snp")
            )
            out[, cn] <- out[, cn] / log(2)
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
    by <- c("gene", "snp", "test")
} else {
    
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"

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
            df_all <- merge(covars, out)
            df_all
        },
        mc.cores = 8
    )

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
    by <- c("gene", "snp", "test", "condition")
}


sample_res_df <- sample_res_df[sample_res_df$n_eff > minEff & sample_res_df$Rhat < maxRhat, ]
if (method == "sampling") {
    approx_res_df <- approx_res_df[approx_res_df$n_eff > minEff & approx_res_df$Rhat < maxRhat, ]
}

approx_res_df$test <- paste(approx_res_df$gene, approx_res_df$snp, sep = "_")
sample_res_df$test <- paste(sample_res_df$Gene_id, sample_res_df$tag , sep = "_")
sample_res_df$gene <- sample_res_df$Gene_id
sample_res_df$snp <- sample_res_df$tag
sample_res_df$Gene_id <- NULL
sample_res_df$tag <- NULL

approx_res_df$null.95 <- sign(approx_res_df$"2.5%") == sign(approx_res_df$"97.5%")
approx_res_df$null.50 <- sign(approx_res_df$"25.0%") == sign(approx_res_df$"75.0%")
sample_res_df$null.95 <- sign(sample_res_df$"log2_aFC_2.5.") == sign(sample_res_df$"log2_aFC_97.5.")
sample_res_df$null.50 <- sign(sample_res_df$"log2_aFC_25.") == sign(sample_res_df$"log2_aFC_75.")

sample_res_df <- sample_res_df[sample_res_df$test %in% approx_res_df$test, ]
im <- match(sample_res_df$test, approx_res_df$test)
approx_res_df <- approx_res_df[im, ]
stopifnot(all(sample_res_df$test == approx_res_df$test))

r <- range(c(sample_res_df$log2_aFC_mean, approx_res_df$mean))

mdf <- merge(approx_res_df, sample_res_df, by = by, suffix = c(".VB", ".HMC"))
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

## not se mean but other hmc se
diag_vars <- c(
    "gene",
    # "Rhat",
    # "n_eff",
    "time", "n_ase",
    "log2_aFC_se_mean", "log2_aFC_sd",
    "mean_count", "sd_count", "n_wt", "n_het", "n_hom"
)

mdf2 <- mdf[mdf$Rhat.VB < maxRhat & mdf$Rhat.HMC < maxRhat & mdf$n_eff.VB > minEff & mdf$n_eff.HMC > minEff, ]

g <- ggplot(mdf2) +
    aes_string("log2_aFC_mean", "mean", colour = "log2_aFC_se_mean") +
    geom_point(size = 0.8) +
    scale_colour_viridis(name = "SE (x)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "Old estimate", y = "New estimate")
ggsave(sprintf("fig/%s/estimates/se1.png", model), width = 7, height = 7)
g <- ggplot(mdf2) +
    aes_string("log2_aFC_mean", "mean", colour = "se_mean") +
    geom_point(size = 0.8) +
    scale_colour_viridis(name = "SE (y)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "Old estimate", y = "New estimate")
ggsave(sprintf("fig/%s/estimates/se2.png", model), width = 7, height = 7)

for (x in diag_vars) {
    scale <- if (x == "gene") scale_colour_discrete(guide="none") else scale_colour_viridis()

    g <- ggplot(mdf) +
        aes_string(x, "abs(discrepancy)") +
        geom_point(size = 0.8) +
        # geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
        labs(x = x, y = "Discrepancy")

    ggsave(sprintf("fig/%s/diag/disc_%s_%s.png", model, x, mtol), width = 7, height = 7)
}

mdf <- mdf[order(mdf$Null.99), ]

g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.99) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/%s/estimates/%s_mcmc_all_categorical_99.png", model, mtol), width = 7, height = 7)

mdf <- mdf[order(mdf$Null.95), ]
g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.95) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/%s/estimates/%s_mcmc_all_categorical_95.png", model, mtol), width = 7, height = 7)

mdf <- mdf[order(mdf$Null.95.50), ]
g <- ggplot(mdf) +
    aes(log2_aFC_mean, mean, colour = Null.95.50) +
    geom_point() +
    scale_colour_brewer(palette = "Paired", limits = levs) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    lims(x = r, y = r) +
    labs(x = "MCMC estimate", y = sprintf("%s estimate", mname))
ggsave(file = sprintf("fig/%s/estimates/%s_mcmc_all_categorical_50.png", model, mtol), width = 7, height = 7)
