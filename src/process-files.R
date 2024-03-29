library("argparse")
library("baseqtl")
library("dplyr")
library("BiocParallel")

parser <- ArgumentParser()
parser$add_argument(
    "-i", "--inference",
    default = "vb",
    type = "character"
)
parser$add_argument(
    "-t", "--tolerance",
    default = 1e-2,
    type = "double"
)
parser$add_argument(
    "-m", "--model",
    default = "noGT",
    type = "character"
)

args <- parser$parse_args()
method <- args[["inference"]]
tol <- args[["tolerance"]]
model <- args[["model"]]
source("src/functions.R")
options(mc.cores = 8)

mtol <- mtol(method, tol)

sfile <- sprintf("rds/%s/sfile.rds", model)
combfile <- sprintf("rds/%s/%s_combined.rds", model, mtol)

if (model == "GT") {
    dir <- paste0(
        "/home/abo27/rds/rds-mrc-bsu/",
        "ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"
    )

    files <- list.files(dir)
    files <- grep("rbias", files, value = TRUE)
    stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
    genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))


    outfiles <- list.files(
        sprintf("rds/GT/%s/", mtol), pattern = "ENSG.*.rds", full.names = TRUE
    )
    genes <- unique(gsub(".*(ENSG\\d+).*", "\\1", outfiles))
    infiles <- sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, genes)

    dfs <- bplapply(
        seq_len(length(genes)),
        function(i) {
            cat(i, "/", length(genes), "\n")
            gene <- genes[[i]]
            infile <- infiles[[i]]
            outfiles <- list.files(
                sprintf("rds/GT/%s/", mtol),
                pattern = paste0(gene, ".*"),
                full.names = TRUE
            )

            outs <- lapply(outfiles, readRDS)
            outs <- lapply(outs, function(x) do.call(rbind, x))
            out <- do.call(rbind, outs)
            out <- as.data.frame(out)
            if (method == "vb") {
                out$niter <- sapply(outs, function(x) max(attr(x, "elbo")[[1]]))
            }
            cn <- setdiff(
                colnames(out),
                c(
                    "n_eff", "Rhat", "niter", "converged", "null.99",
                    "gene", "time", "snp", "variable", "seed", "PEP",
                    "rbias"
                )
            )
            out[, cn] <- out[, cn] / log(2)
            if (method == "vb") {
                out <- cbind(out,
                    max_iter = sapply(outs, function(x) nrow(attr(x, "elbo")))
                )
            }
            inp <- readRDS(infile)
            if (!length(out)) {
                return(list())
            }
            snps <- out$snp
            covars <- lapply(
                unique(snps),
                function(snp) {
                    if (!snp %in% names(inp)) {
                        return(NULL)
                    }
                    x <- inp[[snp]]
                    inp1 <- in.neg.beta.prob.eff2(x)
                    data.frame(
                        gene = gene,
                        snp = snp,
                        n_tot = inp1$N,
                        n_ase = inp1$A,
                        mean_count = mean(log1p(inp1$Y)),
                        sd_count = sd(log1p(inp1$Y)),
                        n_wt = sum(inp1$g == 0),
                        n_het = sum(abs(inp1$g) == 1),
                        p_het = sum(abs(inp1$g) == 1) / inp1$N,
                        n_hom = sum(abs(inp1$g) == 2)
                    )
                }
            )
            covars <- do.call(rbind, covars)
            df <- out
            df <- merge(covars, df)
            df
        }
    )

    df_out <- do.call(rbind, dfs)
    df_out$test <- paste(df_out$gene, df_out$snp, sep = "_")
} else {
    dir <- paste0(
        "/home/abo27/rds/rds-mrc-bsu",
        "/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"
    )

    outfiles <- list.files(
        sprintf("rds/noGT/%s/", mtol), pattern = "ENSG*", full.names = TRUE
    )
    genes <- unique(gsub(".*(ENSG\\d+).*", "\\1", outfiles))

    infiles <- list(
        normal_skin = sprintf(
            "%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, genes
        ),
        Psoriasis_skin = sprintf(
            "%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, genes
        )
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
            p_het = sum(abs(inp1$gase) == 1) / inp1$N,
            n_hom = sum(abs(inp1$gase) == 2)
        )
    }
    dfs <- bplapply(
        seq_len(length(genes)),
        function(i) {
            cat(i, "/", length(genes), "\n")
            infile_norm <- infiles[[1]][[i]]
            infile_pso <- infiles[[2]][[i]]

            outfiles <- list.files(
                sprintf("rds/noGT/%s/", mtol),
                pattern = sprintf("%s_.*.rds", genes[[i]]),
                full.names = TRUE
            )
            outs <- lapply(outfiles, readRDS)
            outs <- lapply(outs, function(x) {
                out <- as.data.frame(do.call(rbind, x))
                if (method == "vb") {
                    out$niter <- sapply(x,
                        function(y) max(attr(y, "elbo")[[1]])
                    )
                }
                out
            })
            out <- do.call(rbind, outs)
            ## should figure out why the non-rbias part is so horrible
            # out <- out[out$rbias, ]

            cn <- setdiff(
                colnames(out),
                c(
                    "n_eff", "Rhat", "null.99", "gene", "time", "snp",
                    "condition", "variable", "seed", "PEP", "rbias"
                )
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
            covars_norm$snp <- rownames(covars_norm)
            covars_pso$snp <- rownames(covars_pso)
            covars_norm$condition <- "normal_skin"
            covars_pso$condition <- "Psoriasis_skin"
            covars <- rbind(covars_norm, covars_pso)
            covars$gene <- genes[[i]]
            df_all <- merge(covars, out, by = c("snp", "gene", "condition"))
            df_all
        }
    )
    df_out <- do.call(rbind, dfs)
    df_out$test <- paste(df_out$gene, df_out$snp, sep = "_")
}

cat("nrow:", nrow(df_out), "\n")
print(combfile)
saveRDS(df_out, combfile)
