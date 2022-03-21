library("argparse")
library("baseqtl")
library("dplyr")

parser <- ArgumentParser()
parser$add_argument(
    "-i", "--inference",
    default = "sampling",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-2,
    type = "double"
)
parser$add_argument(
    "-m", "--model",
    default = "GT",
    type = "character"
)

args <- parser$parse_args()
method <- args[["inference"]]
tol <- args[["tolerance"]]
model <- args[["model"]]

mname <- switch(method,
    "vb" = "ADVI",
    "optimizing" = "MAP",
    "sampling" = "HMC"
)

# mtol <- if (method == "vb") sprintf("vb_%1.0e", tol) else method
mtol <- method

sfile <- sprintf("rds/%s/sfile.rds", model)
combfile <- sprintf("rds/%s/%s_combined.rds", model, method)

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

    sfile <- sprintf("rds/%s/sfile.rds", model)
    combfile <- sprintf("rds/%s/%s_combined.rds", model, method)

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
        },
        mc.cores = 8
    )

    approx_res_df <- do.call(rbind, dfs)
    approx_res_df$test <- paste(approx_res_df$gene, approx_res_df$snp , sep = "_")
    sample_res_df <- do.call(
        rbind,
        lapply(
            genes,
            function(gene) {
                read.table(sprintf("%s/rbias.%s.stan.summary.txt", dir, gene), header=TRUE)
            }
        )
    )
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

}

saveRDS(approx_res_df, combfile)
saveRDS(sample_res_df, sfile)
